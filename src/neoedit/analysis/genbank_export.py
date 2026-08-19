"""Write GenBank records that carry BOTH mitochondrially-translated genes and
non-canonical ORFs translated by cytoplasmic ribosomes.

Mitochondrial-derived peptides (MDPs) — humanin and SHLP1-6 inside MT-RNR2, MOTS-c
inside MT-RNR1, and others — are encoded in mtDNA but translated on cytoplasmic
ribosomes, so they use the *standard* genetic code (transl_table=1) while the
canonical mitochondrial genes around them use the organism's mitochondrial code
(transl_table=2 in vertebrates). NCBI records normally annotate only the latter.

`build_record` merges user-found ORFs into an existing record (or a bare sequence),
tagging each new CDS with its own transl_table plus notes describing where it is
translated, so any GenBank-aware viewer shows both classes in context.
"""
from __future__ import annotations

import io
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .orf_finder import ORF

# Names of described mitochondrial-derived peptides, offered as naming suggestions.
# NOTE: this is a name/host list only. An ORF is NOT called one of these unless the
# user supplies reference peptide sequences and the ORF is similar to one of them
# (see `name_by_similarity`); length alone is not evidence of homology.
KNOWN_MDPS = [
    ("humanin", "MT-RNR2", "cytoprotective MDP encoded within the 16S rRNA gene"),
    ("MOTS-c", "MT-RNR1", "mitochondrial ORF of the 12S rRNA type-c; metabolic regulator"),
    ("SHLP1", "MT-RNR2", "small humanin-like peptide 1"),
    ("SHLP2", "MT-RNR2", "small humanin-like peptide 2"),
    ("SHLP3", "MT-RNR2", "small humanin-like peptide 3"),
    ("SHLP4", "MT-RNR2", "small humanin-like peptide 4"),
    ("SHLP5", "MT-RNR2", "small humanin-like peptide 5"),
    ("SHLP6", "MT-RNR2", "small humanin-like peptide 6"),
    ("SHMOOSE", "MT-ND4", "MDP encoded within ND4"),
]
KNOWN_MDP_NAMES = [n for n, _h, _d in KNOWN_MDPS]

CYTOPLASMIC_NOTE = ("non-canonical ORF within the mitochondrial genome; translated by cytoplasmic "
                    "ribosomes using the standard genetic code (transl_table=1), unlike the "
                    "surrounding mitochondrially-translated genes")


@dataclass
class ORFAnnotation:
    """One ORF to write, with the qualifiers a reader needs to interpret it."""
    orf: ORF
    name: str = ""                     # /gene
    product: str = ""                  # /product
    table: int = 1                     # /transl_table for THIS feature
    cytoplasmic: bool = True           # adds the explanatory note
    feature_type: str = "CDS"
    note: str = ""
    extra: dict = field(default_factory=dict)

    def qualifiers(self, translation: str) -> dict:
        q: dict[str, list[str]] = {}
        if self.name:
            q["gene"] = [self.name]
        q["product"] = [self.product or self.name or "hypothetical peptide"]
        q["codon_start"] = ["1"]
        q["transl_table"] = [str(self.table)]
        q["translation"] = [translation]
        notes = []
        if self.cytoplasmic:
            notes.append(CYTOPLASMIC_NOTE)
        if self.orf.partial5 or self.orf.partial3:
            notes.append("partial ORF at a sequence end")
        if self.note:
            notes.append(self.note)
        notes.append(f"identified with NeoEdit ORF finder (genetic code {self.table}, "
                     f"start codon {self.orf.start_codon}, {self.orf.length_aa} aa)")
        q["note"] = notes
        for k, v in self.extra.items():
            q[k] = [str(v)]
        return q


def translate_orf(nt: str, table: int) -> str:
    s = nt.upper().replace("U", "T")
    s = s[: len(s) - (len(s) % 3)]
    if not s:
        return ""
    aa = str(Seq(s).translate(table=table))
    return aa[:-1] if aa.endswith("*") else aa


def orf_to_feature(a: ORFAnnotation, offset: int = 0) -> SeqFeature:
    o = a.orf
    loc = SimpleLocation(o.start + offset, o.end + offset, strand=1 if o.strand > 0 else -1)
    tr = translate_orf(o.nt, a.table)
    return SeqFeature(loc, type=a.feature_type, qualifiers=a.qualifiers(tr))


def load_source_record(path: str, index: int = 0) -> SeqRecord:
    with open(path, "r", errors="replace") as fh:
        recs = list(SeqIO.parse(fh, "genbank"))
    if not recs:
        raise ValueError(f"No GenBank records in {path}")
    return recs[index]


def build_record(sequence: str, annotations: Sequence[ORFAnnotation],
                 source: Optional[SeqRecord] = None,
                 record_id: str = "", description: str = "",
                 organism: str = "", mol_type: str = "genomic DNA",
                 topology: str = "circular", default_table: int = 2,
                 keep_source_features: bool = True) -> SeqRecord:
    """Merge `annotations` into `source` (or create a fresh record).

    Existing features are kept verbatim, so canonical mitochondrial genes retain
    their own transl_table; the added ORFs carry theirs.
    """
    seq = Seq("".join(c for c in sequence if c not in "-.~").upper())
    if source is not None:
        rec = SeqRecord(seq if len(seq) == len(source.seq) else source.seq,
                        id=record_id or source.id, name=(record_id or source.name or "record")[:16],
                        description=description or source.description)
        rec.annotations = dict(source.annotations)
        if keep_source_features:
            rec.features = list(source.features)
    else:
        rec = SeqRecord(seq, id=record_id or "NEOEDIT_1", name=(record_id or "NEOEDIT_1")[:16],
                        description=description or "sequence annotated with NeoEdit")
        rec.annotations = {}
        rec.features = [SeqFeature(SimpleLocation(0, len(seq), strand=1), type="source",
                                   qualifiers={"organism": [organism or "unspecified"],
                                               "mol_type": [mol_type],
                                               **({"organelle": ["mitochondrion"]} if default_table == 2 else {})})]
    rec.annotations.setdefault("molecule_type", "DNA")
    rec.annotations.setdefault("topology", topology)
    if organism:
        rec.annotations["organism"] = organism
    rec.annotations.setdefault("data_file_division", "INV" if default_table == 2 else "UNK")
    existing = {(int(f.location.start), int(f.location.end), f.type,
                 tuple(f.qualifiers.get("gene", f.qualifiers.get("product", []))))
                for f in rec.features}
    added = 0
    for a in annotations:
        feat = orf_to_feature(a)
        key = (int(feat.location.start), int(feat.location.end), feat.type,
               tuple(feat.qualifiers.get("gene", feat.qualifiers.get("product", []))))
        if key in existing:
            continue
        rec.features.append(feat)
        existing.add(key)
        added += 1
    rec.features.sort(key=lambda f: (int(f.location.start), -int(f.location.end)))
    rec.annotations.setdefault("comment", "")
    note = (f"{added} non-canonical ORF(s) added with NeoEdit. Features carry per-feature "
            f"/transl_table qualifiers: mitochondrially-translated genes use the organism's "
            f"mitochondrial code, while cytoplasmically-translated ORFs (e.g. mitochondrial-derived "
            f"peptides) use the standard code.")
    rec.annotations["comment"] = (rec.annotations["comment"] + "\n" + note).strip() if rec.annotations["comment"] else note
    return rec


def write_genbank(rec: SeqRecord, path: str):
    with open(path, "w") as fh:
        SeqIO.write(rec, fh, "genbank")


def dumps_genbank(rec: SeqRecord) -> str:
    out = io.StringIO()
    SeqIO.write(rec, out, "genbank")
    return out.getvalue()


def suggest_mdp_name(o: ORF, host_gene: str = "") -> tuple[str, str]:
    """A neutral, descriptive name for an ORF. Deliberately does NOT guess published
    MDP identities from length — use `name_by_similarity` with reference peptides for that."""
    base = f"ORF{o.length_aa}aa"
    if host_gene:
        return f"{base}_in_{host_gene}", f"non-canonical ORF within {host_gene}"
    return base, "non-canonical ORF"


def load_reference_peptides(path: str) -> dict[str, str]:
    """Reference peptides (FASTA) used to name ORFs by similarity, e.g. human humanin,
    MOTS-c, SHLP1-6 sequences from the literature."""
    from Bio import SeqIO as _S
    out = {}
    with open(path, "r", errors="replace") as fh:
        for rec in _S.parse(fh, "fasta"):
            out[rec.id] = str(rec.seq).upper().replace("*", "")
    return out


def peptide_identity(a: str, b: str) -> float:
    """Global-alignment identity of two peptides (0-1)."""
    from Bio.Align import PairwiseAligner
    from Bio.Align import substitution_matrices
    if not a or not b:
        return 0.0
    al = PairwiseAligner()
    al.mode = "global"
    try:
        al.substitution_matrix = substitution_matrices.load("BLOSUM62")
        al.open_gap_score, al.extend_gap_score = -11, -1
    except Exception:
        al.match_score, al.mismatch_score = 1, -1
        al.open_gap_score, al.extend_gap_score = -2, -0.5
    try:
        aln = al.align(a, b)[0]
    except Exception:
        return 0.0
    ta, tb = aln[0], aln[1]
    same = sum(1 for x, y in zip(ta, tb) if x == y and x != "-")
    return same / max(len(a), len(b))


def name_by_similarity(orfs: Sequence[ORF], references: dict[str, str],
                       min_identity: float = 0.5) -> dict[int, tuple[str, float, str]]:
    """Match ORF peptides against reference peptides.

    Returns {index_in_orfs: (reference_name, identity, note)} for ORFs whose best
    match reaches `min_identity`. Evidence-based alternative to guessing by length.
    """
    out = {}
    for i, o in enumerate(orfs):
        pep = o.aa.rstrip("*")
        best_name, best_id = "", 0.0
        for name, ref in references.items():
            ident = peptide_identity(pep, ref)
            if ident > best_id:
                best_name, best_id = name, ident
        if best_name and best_id >= min_identity:
            out[i] = (best_name, best_id,
                      f"similar to {best_name} ({best_id:.0%} aa identity over {len(pep)} aa)")
    return out
