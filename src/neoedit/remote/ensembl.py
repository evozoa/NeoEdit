"""Ensembl REST client pinned to a release (default 116).

Ensembl keeps one REST server per release (`e115.rest.ensembl.org`, …) next to the
live `rest.ensembl.org`. `EnsemblClient.resolve()` tries the archive host for the wanted
release first and falls back to the live server, reporting the release it actually got,
so the caller can warn when the data is not from the release the user asked for.

Genomic imports are turned into a GenBank record with gene / mRNA|tRNA|rRNA|ncRNA / CDS
features built from the `/overlap` endpoint, so NeoEdit's gene view and the amino-acid
line work on Ensembl data exactly as on an NCBI record."""
from __future__ import annotations

import json
import re
import urllib.parse
from dataclasses import dataclass, field
from datetime import date

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, BeforePosition, AfterPosition, ExactPosition

from .http import http_get, RemoteError, write_download

DEFAULT_RELEASE = 116
LIVE_SERVER = "https://rest.ensembl.org"
ARCHIVE_SERVER = "https://e{release}.rest.ensembl.org"
GRCH37_SERVER = "https://grch37.rest.ensembl.org"
DIVISIONS = [("Vertebrates", "EnsemblVertebrates"), ("Plants", "EnsemblPlants"), ("Metazoa", "EnsemblMetazoa"),
             ("Fungi", "EnsemblFungi"), ("Protists", "EnsemblProtists")]
SEQ_TYPES = [("Genomic (with gene models)", "genomic"), ("cDNA (spliced transcript)", "cdna"),
             ("CDS (coding sequence)", "cds"), ("Protein", "protein")]
MAX_REGION = 5_000_000          # /overlap caps at 5 Mb; /sequence/region at 10 Mb
_MITO = {"mt", "chrm", "chrmt", "mito", "mitochondrion", "mitochondrion_genome", "m"}
_ID_RE = re.compile(r"^(ENS[A-Z]*[GTPE]\d{6,}|[A-Z]{2,}\w*[GTP]\d{5,})(\.\d+)?$", re.I)   # ENSG…, ENSDARG…, FBgn (loosely)


@dataclass
class Species:
    name: str
    display_name: str
    common_name: str = ""
    assembly: str = ""
    division: str = ""

    def label(self) -> str:
        cn = f" ({self.common_name})" if self.common_name and self.common_name.lower() != self.display_name.lower() else ""
        return f"{self.display_name}{cn} — {self.name} [{self.assembly}]"


@dataclass
class Lookup:
    id: str
    object_type: str
    species: str
    seq_region_name: str
    start: int
    end: int
    strand: int
    display_name: str = ""
    biotype: str = ""
    description: str = ""
    assembly_name: str = ""
    parent: str = ""
    raw: dict = field(default_factory=dict)

    def region(self) -> str:
        return f"{self.seq_region_name}:{self.start}-{self.end}:{'+' if self.strand >= 0 else '-'}"


def parse_region(text: str) -> tuple[str, int, int, int]:
    """'17:43,044,295-43,170,245', 'chr17:1..100:-1', 'MT:1-16569:-' -> (chrom, start, end, strand); 1-based inclusive."""
    t = (text or "").strip().replace(",", "").replace(" ", "")
    m = re.match(r"^([\w.-]+):(\d+)(?:\.\.|-)(\d+)(?::(-1|1|\+|-))?$", t)
    if not m:
        raise RemoteError(f"Cannot parse region '{text}'. Use chromosome:start-end[:strand], e.g. 17:43044295-43170245:-1")
    chrom, s, e, st = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    if s > e:
        s, e = e, s
    if s < 1:
        raise RemoteError("Region start must be ≥ 1 (coordinates are 1-based).")
    strand = -1 if st in ("-1", "-") else 1
    return chrom, s, e, strand


def normalize_species(s: str) -> str:
    return re.sub(r"\s+", "_", (s or "").strip()).lower()


def looks_like_id(q: str) -> bool:
    return bool(_ID_RE.match((q or "").strip()))


class EnsemblClient:
    def __init__(self, release: int | None = DEFAULT_RELEASE, server: str | None = None, fetch=http_get):
        self.release_wanted = release
        self.server = server.rstrip("/") if server else None
        self.release = None                  # release actually served (after resolve)
        self._fetch = fetch
        self._resolved = False

    # ------------------------------------------------------------ plumbing
    def _raw(self, server: str, path: str, content_type: str, params: dict | None = None,
             timeout: float = 60.0, retries: int = 3) -> bytes:
        q = dict(params or {})
        q["content-type"] = content_type
        url = server + path + "?" + urllib.parse.urlencode(q, doseq=True, safe=":,")
        return self._fetch(url, headers={"Accept": content_type}, timeout=timeout, retries=retries)

    def _software_release(self, server: str) -> int | None:
        try:
            d = json.loads(self._raw(server, "/info/software", "application/json", timeout=15, retries=0))
            return int(d.get("release"))
        except Exception:
            return None

    def resolve(self) -> tuple[str, int | None]:
        """Pick the server: explicit > archive for the wanted release > live. Returns (server, release)."""
        if self._resolved:
            return self.server, self.release
        cands = []
        if self.server:
            cands.append(self.server)
        else:
            if self.release_wanted:
                cands.append(ARCHIVE_SERVER.format(release=self.release_wanted))
            cands.append(LIVE_SERVER)
        chosen, rel = None, None
        for srv in cands:
            r = self._software_release(srv)
            if r is None:
                continue
            chosen, rel = srv, r
            if not self.release_wanted or r == self.release_wanted:
                break
            if srv == self.server:
                break
        if chosen is None:
            raise RemoteError("Cannot reach the Ensembl REST servers (" + ", ".join(cands) + ").")
        self.server, self.release, self._resolved = chosen, rel, True
        return chosen, rel

    def release_note(self) -> str:
        srv, rel = self.resolve()
        host = urllib.parse.urlsplit(srv).netloc
        if self.release_wanted and rel != self.release_wanted:
            return f"Ensembl release {rel} ({host}) — release {self.release_wanted} is not available on this server"
        return f"Ensembl release {rel} ({host})"

    def get_json(self, path: str, **params):
        srv, _ = self.resolve()
        body = self._raw(srv, path, "application/json", params)
        try:
            return json.loads(body)
        except ValueError:
            raise RemoteError(f"Ensembl returned non-JSON for {path}: {body[:120]!r}") from None

    def get_text(self, path: str, content_type: str = "text/x-fasta", **params) -> str:
        srv, _ = self.resolve()
        return self._raw(srv, path, content_type, params).decode("utf-8", "replace")

    # ------------------------------------------------------------ API
    def species(self, division: str = "EnsemblVertebrates") -> list[Species]:
        d = self.get_json("/info/species", division=division)
        out = []
        for s in d.get("species", []):
            out.append(Species(s.get("name", ""), s.get("display_name") or s.get("name", ""),
                               s.get("common_name") or "", s.get("assembly") or "", s.get("division") or division))
        out.sort(key=lambda x: x.display_name.lower())
        return out

    def lookup(self, query: str, species: str | None = None) -> Lookup:
        q = (query or "").strip()
        if not q:
            raise RemoteError("Enter a gene symbol or an Ensembl stable ID.")
        if looks_like_id(q):
            d = self.get_json(f"/lookup/id/{q.split('.')[0]}", expand=0)
        else:
            sp = normalize_species(species or "")
            if not sp:
                raise RemoteError("A species is needed to look up a gene symbol.")
            d = self.get_json(f"/lookup/symbol/{sp}/{urllib.parse.quote(q)}", expand=0)
        if isinstance(d, dict) and "seq_region_name" not in d and d.get("Parent") and species is not None:
            # a Translation (ENSP…) has no location of its own: use its transcript
            parent = self.lookup(d["Parent"], species)
            parent.raw = {"translation": d, **parent.raw}
            return parent
        if not isinstance(d, dict) or "seq_region_name" not in d:
            raise RemoteError(f"No location found for '{q}'.")
        return Lookup(id=d.get("id", q), object_type=d.get("object_type", ""), species=d.get("species", species or ""),
                      seq_region_name=str(d["seq_region_name"]), start=int(d["start"]), end=int(d["end"]),
                      strand=int(d.get("strand", 1)), display_name=d.get("display_name", ""),
                      biotype=d.get("biotype", ""), description=d.get("description") or "",
                      assembly_name=d.get("assembly_name", ""), parent=d.get("Parent", ""), raw=d)

    def sequence_by_id(self, stable_id: str, seq_type: str = "cds", mask: str | None = None) -> str:
        """FASTA text for an ID. Gene IDs with cds/cdna/protein return every transcript."""
        params = {"type": seq_type, "multiple_sequences": 1}
        if mask and seq_type == "genomic":
            params["mask"] = mask
        return self.get_text(f"/sequence/id/{stable_id}", **params)

    def sequence_region(self, species: str, chrom: str, start: int, end: int, strand: int = 1,
                        mask: str | None = None, coord_system_version: str | None = None) -> str:
        if end - start + 1 > 10_000_000:
            raise RemoteError("Ensembl serves at most 10 Mb of sequence per request.")
        params = {}
        if mask:
            params["mask"] = mask
        if coord_system_version:
            params["coord_system_version"] = coord_system_version
        sp = normalize_species(species)
        return self.get_text(f"/sequence/region/{sp}/{chrom}:{start}..{end}:{strand}", **params)

    def overlap_region(self, species: str, chrom: str, start: int, end: int,
                       features=("gene", "transcript", "exon", "cds")) -> list[dict]:
        if end - start + 1 > MAX_REGION:
            raise RemoteError(f"Gene models can be fetched for at most {MAX_REGION // 1_000_000} Mb at a time.")
        sp = normalize_species(species)
        d = self.get_json(f"/overlap/region/{sp}/{chrom}:{start}..{end}", feature=list(features))
        if isinstance(d, dict) and d.get("error"):
            raise RemoteError(f"Ensembl: {d['error']}")
        return list(d)

    # ------------------------------------------------------------ high level
    def fetch_genomic(self, species: str, chrom: str, start: int, end: int, strand: int = 1,
                      annotate: bool = True, mask: str | None = None, name: str | None = None,
                      description: str = "", gene: Lookup | None = None) -> SeqRecord:
        """Region (1-based inclusive) as an annotated SeqRecord oriented to `strand`."""
        fasta = self.sequence_region(species, chrom, start, end, 1, mask)
        seq = "".join(ln.strip() for ln in fasta.splitlines() if not ln.startswith(">"))
        if not seq:
            raise RemoteError(f"Ensembl returned no sequence for {chrom}:{start}-{end} ({species}).")
        feats = self.overlap_region(species, chrom, start, end) if annotate else []
        asm = ""
        m = re.search(r"chromosome:([^:]+):", fasta)
        if m:
            asm = m.group(1)
        _srv, rel = self.resolve()
        return build_record(seq, chrom, start, end, feats, species=species, assembly=asm, strand=strand,
                            name=name, description=description, release=rel, gene=gene, masked=bool(mask))

    def fetch_sequences(self, stable_id: str, seq_type: str, label: str = "", species: str = "") -> list[SeqRecord]:
        """cDNA / CDS / protein FASTA for an ID (all transcripts of a gene) as SeqRecords."""
        txt = self.sequence_by_id(stable_id, seq_type)
        recs = _parse_fasta(txt)
        if not recs:
            raise RemoteError(f"Ensembl returned no {seq_type} sequence for {stable_id}.")
        _srv, rel = self.resolve()
        for r in recs:
            bits = [label, seq_type, species.replace("_", " "), f"Ensembl {rel}", stable_id]
            r.description = " ".join(b for b in bits if b)
            r.annotations["molecule_type"] = "protein" if seq_type == "protein" else "DNA"
        return recs


# ---------------------------------------------------------------- record building (pure)
def _parse_fasta(text: str) -> list[SeqRecord]:
    recs, name, desc, buf = [], None, "", []
    for ln in text.splitlines():
        if ln.startswith(">"):
            if name is not None:
                recs.append(SeqRecord(Seq("".join(buf)), id=name, name=name[:16], description=desc))
            hdr = ln[1:].strip().split(None, 1)
            name, desc, buf = (hdr[0] if hdr else "seq"), (hdr[1] if len(hdr) > 1 else ""), []
        else:
            buf.append(ln.strip())
    if name is not None:
        recs.append(SeqRecord(Seq("".join(buf)), id=name, name=name[:16], description=desc))
    return recs


_TRANSCRIPT_TYPE = {"protein_coding": "mRNA", "Mt_tRNA": "tRNA", "tRNA": "tRNA", "Mt_rRNA": "rRNA", "rRNA": "rRNA"}
_NCRNA_CLASS = {"lncRNA": "lncRNA", "miRNA": "miRNA", "snoRNA": "snoRNA", "snRNA": "snRNA", "scaRNA": "scaRNA",
                "misc_RNA": "other", "ribozyme": "ribozyme", "vault_RNA": "vault_RNA", "sRNA": "other",
                "scRNA": "scRNA", "Y_RNA": "Y_RNA", "piRNA": "piRNA", "siRNA": "siRNA"}


def _loc(s0: int, e0: int, strand: int, L: int, flip: bool, clip5: bool, clip3: bool) -> SimpleLocation:
    """0-based half-open piece -> SimpleLocation, mirrored when `flip` (minus-strand output)."""
    if flip:
        s0, e0, strand = L - e0, L - s0, -strand
        clip5, clip3 = clip5, clip3          # biological 5'/3' unchanged by mirroring
    # which end is biological 5'? plus: start; minus: end
    if strand >= 0:
        sp = BeforePosition(s0) if clip5 else ExactPosition(s0)
        ep = AfterPosition(e0) if clip3 else ExactPosition(e0)
    else:
        sp = BeforePosition(s0) if clip3 else ExactPosition(s0)
        ep = AfterPosition(e0) if clip5 else ExactPosition(e0)
    return SimpleLocation(sp, ep, strand)


def _pieces(items: list[dict], rstart: int, rend: int, L: int, flip: bool) -> tuple[list[SimpleLocation], int]:
    """Clip Ensembl features (1-based inclusive) to the region, order biologically, return
    (locations, codon_start) where codon_start uses the GFF phase of the first surviving piece."""
    if not items:
        return [], 1
    strand = int(items[0].get("strand", 1))
    items = sorted(items, key=lambda d: int(d["start"]), reverse=(strand < 0))   # biological order
    kept = []                                   # (cs, ce, clip5, clip3)
    dropped_before = dropped_after = False
    for it in items:
        s, e = int(it["start"]), int(it["end"])
        cs, ce = max(s, rstart), min(e, rend)
        if ce < cs:                             # whole piece outside the region
            if kept:
                dropped_after = True
            else:
                dropped_before = True
            continue
        clip_lo, clip_hi = cs > s, ce < e
        clip5 = clip_lo if strand >= 0 else clip_hi
        clip3 = clip_hi if strand >= 0 else clip_lo
        kept.append([it, cs, ce, clip5, clip3])
    if not kept:
        return [], 1
    if dropped_before:
        kept[0][3] = True                       # a 5' piece is missing: partial at the 5' end
    if dropped_after:
        kept[-1][4] = True
    it, cs, ce, _c5, _c3 = kept[0]
    s, e = int(it["start"]), int(it["end"])
    d5 = (cs - s) if strand >= 0 else (e - ce)
    try:
        phase = int(it.get("phase", 0) or 0)
    except (TypeError, ValueError):
        phase = 0
    codon_start = ((phase - d5) % 3) + 1 if d5 or phase else 1
    locs = [_loc(cs - rstart, ce - rstart + 1, strand, L, flip, c5, c3) for _it, cs, ce, c5, c3 in kept]
    return locs, codon_start


def _compound(locs: list[SimpleLocation]):
    return locs[0] if len(locs) == 1 else CompoundLocation(locs)


def build_record(seq: str, chrom: str, start: int, end: int, feats: list[dict], *, species: str = "",
                 assembly: str = "", strand: int = 1, name: str | None = None, description: str = "",
                 release: int | None = None, gene: Lookup | None = None, masked: bool = False) -> SeqRecord:
    """Assemble a GenBank-style SeqRecord from a + strand region sequence and /overlap features.

    `strand=-1` reverse-complements the sequence and mirrors every feature, so a minus-strand
    gene reads left to right. Features straddling the region edge are clipped and marked
    partial (< / >), with /codon_start derived from the GFF phase."""
    L = len(seq)
    flip = strand < 0
    if flip:
        seq = str(Seq(seq).reverse_complement())
    sp_disp = species.replace("_", " ").strip()
    sp_disp = sp_disp[:1].upper() + sp_disp[1:]
    region = f"{chrom}:{start}-{end}({'-' if flip else '+'})"
    if gene is not None and gene.display_name:
        rid = gene.display_name
    elif gene is not None:
        rid = gene.id
    else:
        rid = f"{chrom}_{start}-{end}"
    rid = name or rid
    desc_bits = [description] if description else []
    if gene is not None:
        desc_bits.append(f"{gene.id} {gene.biotype}".strip())
    desc_bits.append(f"{sp_disp} {assembly}".strip())
    desc_bits.append(region)
    if release:
        desc_bits.append(f"Ensembl {release}")
    if masked:
        desc_bits.append("soft-masked")
    rec = SeqRecord(Seq(seq), id=rid, name=re.sub(r"[^\w.-]", "_", rid)[:16], description=" ".join(desc_bits))
    rec.annotations.update({"molecule_type": "DNA", "topology": "linear", "organism": sp_disp or ".",
                            "source": f"{sp_disp} {assembly}".strip() or ".",
                            "date": date.today().strftime("%d-%b-%Y").upper(),
                            "data_file_division": "UNK",
                            "comment": f"Imported by NeoEdit from Ensembl{' release ' + str(release) if release else ''}: "
                                       f"{sp_disp} {assembly} {region}"})
    mito = chrom.lower() in _MITO
    # source feature
    src = SeqFeature(SimpleLocation(0, L, 1), type="source",
                     qualifiers={"organism": [sp_disp or "unknown"], "mol_type": ["genomic DNA"],
                                 "chromosome": [chrom], "note": [f"{assembly} {region}".strip()]})
    if mito:
        src.qualifiers["organelle"] = ["mitochondrion"]
    rec.features.append(src)

    genes = [f for f in feats if f.get("feature_type") == "gene"]
    transcripts = [f for f in feats if f.get("feature_type") == "transcript"]
    exons = [f for f in feats if f.get("feature_type") == "exon"]
    cdss = [f for f in feats if f.get("feature_type") == "cds"]
    ex_by_t: dict[str, list[dict]] = {}
    for e in exons:
        ex_by_t.setdefault(e.get("Parent", ""), []).append(e)
    cds_by_t: dict[str, list[dict]] = {}
    for c in cdss:
        cds_by_t.setdefault(c.get("Parent", ""), []).append(c)
    gname = {g["id"]: (g.get("external_name") or g["id"]) for g in genes}
    tr_by_g: dict[str, list[dict]] = {}
    for t in transcripts:
        tr_by_g.setdefault(t.get("Parent", ""), []).append(t)

    ordered: list[SeqFeature] = []
    vert_mito_table = "2" if mito else None
    for g in sorted(genes, key=lambda d: int(d["start"])):
        locs, _ = _pieces([g], start, end, L, flip)
        if not locs:
            continue
        gn = gname[g["id"]]
        q = {"gene": [gn], "db_xref": [f"Ensembl:{g['id']}"]}
        notes = [b for b in (g.get("biotype", ""), g.get("description") or "") if b]
        if notes:
            q["note"] = ["; ".join(notes)]
        ordered.append(SeqFeature(_compound(locs), type="gene", qualifiers=q))
        for t in sorted(tr_by_g.get(g["id"], []), key=lambda d: (0 if d.get("is_canonical") else 1, int(d["start"]))):
            bt = t.get("biotype", "")
            ttype = _TRANSCRIPT_TYPE.get(bt, "ncRNA")
            ex = ex_by_t.get(t["id"]) or [t]
            tl, _ = _pieces(ex, start, end, L, flip)
            if not tl:
                continue
            tq = {"gene": [gn], "transcript_id": [f"{t['id']}.{t['version']}" if t.get("version") else t["id"]],
                  "db_xref": [f"Ensembl:{t['id']}"]}
            tn = [bt] if bt else []
            if t.get("is_canonical"):
                tn.append("Ensembl canonical")
            if tn:
                tq["note"] = ["; ".join(tn)]
            if ttype == "ncRNA":
                tq["ncRNA_class"] = [_NCRNA_CLASS.get(bt, "other")]
            if ttype in ("tRNA", "rRNA", "ncRNA"):
                tq["product"] = [t.get("external_name") or gn]
            ordered.append(SeqFeature(_compound(tl), type=ttype, qualifiers=tq))
            cds = cds_by_t.get(t["id"])
            if cds:
                cl, codon_start = _pieces(cds, start, end, L, flip)
                if cl:
                    cq = {"gene": [gn], "codon_start": [codon_start],
                          "transcript_id": tq["transcript_id"]}
                    pid = cds[0].get("protein_id")
                    if pid:
                        cq["protein_id"] = [pid]
                    if vert_mito_table:
                        cq["transl_table"] = [vert_mito_table]
                    if t.get("external_name"):
                        cq["product"] = [t["external_name"]]
                    ordered.append(SeqFeature(_compound(cl), type="CDS", qualifiers=cq))
    # transcripts without a gene in the set (rare: gene starts outside the 5 Mb overlap limit) are skipped
    ordered.sort(key=lambda f: (int(f.location.start), 0 if f.type == "gene" else 1))
    rec.features.extend(ordered)
    return rec


def save_record(rec: SeqRecord, out_dir: str, fmt: str = "genbank") -> str:
    from Bio import SeqIO
    import io as _io
    buf = _io.StringIO()
    SeqIO.write(rec, buf, fmt)
    ext = ".gb" if fmt == "genbank" else ".fasta"
    return write_download(buf.getvalue(), out_dir, rec.id, ext)


def save_records(recs: list[SeqRecord], out_dir: str, stem: str) -> str:
    from Bio import SeqIO
    import io as _io
    buf = _io.StringIO()
    SeqIO.write(recs, buf, "fasta")
    return write_download(buf.getvalue(), out_dir, stem, ".fasta")
