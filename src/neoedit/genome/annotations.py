"""Gene-model annotations (GFF3 / GTF / BED) and synteny blocks (PAF), with interval queries."""
from __future__ import annotations

import bisect
import gzip
import re
from dataclasses import dataclass, field
from typing import Iterable, Iterator


@dataclass
class Transcript:
    id: str
    name: str
    start: int            # 0-based half-open, genomic
    end: int
    strand: int
    biotype: str = ""
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)

    def utrs(self) -> list[tuple[int, int]]:
        """Exon segments not covered by CDS."""
        if not self.cds:
            return []
        out = []
        for es, ee in self.exons:
            segs = [(es, ee)]
            for cs, ce in self.cds:
                nsegs = []
                for s, e in segs:
                    if ce <= s or cs >= e:
                        nsegs.append((s, e))
                    else:
                        if s < cs:
                            nsegs.append((s, cs))
                        if ce < e:
                            nsegs.append((ce, e))
                segs = nsegs
            out += segs
        return out


@dataclass
class Gene:
    id: str
    name: str
    seqid: str
    start: int
    end: int
    strand: int
    biotype: str = ""
    transcripts: list[Transcript] = field(default_factory=list)
    attrs: dict = field(default_factory=dict)

    def __len__(self):
        return self.end - self.start

    @property
    def cytoplasmic(self) -> bool:
        """CDS translated with the standard code inside an organelle genome, i.e. a
        mitochondrial-derived peptide or similar non-canonical ORF."""
        return self.attrs.get("transl_table") == "1"

    @property
    def low_confidence(self) -> bool:
        """Liftoff flags: partial mapping or low identity."""
        a = self.attrs
        return a.get("partial_mapping") == "True" or a.get("low_identity") == "True"


@dataclass
class SyntenyBlock:
    qname: str; qstart: int; qend: int; strand: int
    tname: str; tstart: int; tend: int
    matches: int; alen: int; mapq: int

    @property
    def identity(self):
        return self.matches / self.alen if self.alen else 0.0


def _open(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, "r", errors="replace")


def _attrs_gff3(s: str) -> dict:
    out = {}
    for part in s.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip().replace("%2C", ",").replace("%3B", ";").replace("%3D", "=")
    return out


def _attrs_gtf(s: str) -> dict:
    out = {}
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', s):
        out.setdefault(m.group(1), m.group(2))
    return out


class Annotation:
    """All genes of one or more seqids, with fast overlap queries per seqid."""

    def __init__(self):
        self.genes_by_seq: dict[str, list[Gene]] = {}
        self._starts: dict[str, list[int]] = {}
        self._maxlen: dict[str, int] = {}
        self._wrapping: dict[str, list[Gene]] = {}    # genes whose end > sequence length (cross the origin)
        self.by_name: dict[str, Gene] = {}
        self.source: str = ""
        self.lengths: dict[str, int] = {}             # seqid -> molecule length, when known
        self.circular: dict[str, bool] = {}           # seqid -> declared circular topology

    # ------------------------------------------------------------ topology
    def set_topology(self, seqid: str, length: int | None = None, circular: bool | None = None):
        """Declare a molecule's length and/or circularity (re-finalizes so wrap-around genes index correctly)."""
        if length:
            self.lengths[seqid] = int(length)
        if circular is not None:
            self.circular[seqid] = bool(circular)
        self.finalize()

    def is_circular(self, seqid: str) -> bool:
        return bool(self.circular.get(seqid))

    def length_of(self, seqid: str) -> int | None:
        return self.lengths.get(seqid)

    # ------------------------------------------------------------ building
    def add_gene(self, g: Gene):
        self.genes_by_seq.setdefault(g.seqid, []).append(g)

    def finalize(self):
        for sid, genes in self.genes_by_seq.items():
            L = self.lengths.get(sid)
            if L and self.circular.get(sid):
                for g in genes:
                    unwrap_gene(g, L)
            genes.sort(key=lambda g: (g.start, g.end))
            self._starts[sid] = [g.start for g in genes]
            self._maxlen[sid] = max((len(g) for g in genes), default=0)
            self._wrapping[sid] = [g for g in genes if L and g.end > L]
            for g in genes:
                for t in g.transcripts:
                    t.exons.sort(); t.cds.sort()
                g.transcripts.sort(key=lambda t: (t.start, t.end))
                self.by_name[g.name.lower()] = g
                self.by_name.setdefault(g.id.lower(), g)

    # ------------------------------------------------------------ queries
    def seqids(self) -> list[str]:
        return list(self.genes_by_seq)

    def overlapping(self, seqid: str, start: int, end: int) -> list[Gene]:
        genes = self.genes_by_seq.get(seqid, [])
        if not genes:
            return []
        starts = self._starts[seqid]
        lo = bisect.bisect_left(starts, start - self._maxlen[seqid])
        hi = bisect.bisect_right(starts, end)
        out = [g for g in genes[lo:hi] if g.end > start and g.start < end]
        L = self.lengths.get(seqid)
        if L:
            # genes crossing the origin also occupy [0, end - L)
            for g in self._wrapping.get(seqid, ()):
                if g.end - L > start and g not in out:
                    out.append(g)
            out.sort(key=lambda g: (g.start, g.end))
        return out

    def find(self, text: str) -> list[Gene]:
        t = text.lower().strip()
        if not t:
            return []
        if t in self.by_name:
            return [self.by_name[t]]
        return [g for genes in self.genes_by_seq.values() for g in genes
                if t in g.name.lower() or t in g.id.lower() or any(t in v.lower() for v in g.attrs.values())][:200]

    def count(self) -> int:
        return sum(len(v) for v in self.genes_by_seq.values())


def split_span(start: int, end: int, length: int | None) -> list[tuple[int, int]]:
    """[start,end) in unwrapped coordinates -> the 1 or 2 pieces that lie on the molecule [0,length)."""
    if not length or end <= length:
        return [(start, end)]
    pieces = []
    if start < length:
        pieces.append((start, length))
    pieces.append((0 if start < length else start - length, end - length))
    return pieces


def fmt_span(start: int, end: int, length: int | None = None) -> str:
    """1-based display of [start,end); a span past the origin reads '16,024-576 (across origin)'."""
    if length and end > length:
        return f"{start + 1:,}-{end - length:,} (across origin)"
    return f"{start + 1:,}-{end:,}"


def _unwrap_parts(parts: list[tuple[int, int]], length: int) -> tuple[list[tuple[int, int]], bool]:
    """Pieces of one feature in transcription order (as GenBank/GFF list them) -> unwrapped
    coordinates: once the coordinates step backwards past the origin, `length` is added to
    every following piece. Returns (pieces, wrapped?)."""
    out, shift, wrapped = [], 0, False
    prev_end = None
    for s, e in parts:
        if prev_end is not None and s + shift < prev_end:
            shift += length; wrapped = True
        out.append((s + shift, e + shift))
        prev_end = e + shift
    return out, wrapped


def unwrap_gene(g: Gene, length: int):
    """GFF/BED-style genes on a circular molecule: a transcript whose exons jump back to the
    start of the sequence is re-expressed in unwrapped coordinates (end > length).
    Idempotent; exact when the feature was built from ordered parts (see load_genbank)."""
    if not g.transcripts or g.attrs.get("wraps_origin") == "true":
        return
    changed = False
    for t in g.transcripts:
        ex = sorted(t.exons)
        if len(ex) < 2 or (ex[-1][1] - ex[0][0]) <= 0.5 * length:
            continue
        # biggest gap between consecutive exons is the fake "intron" across the origin
        gaps = [(ex[i + 1][0] - ex[i][1], i) for i in range(len(ex) - 1)]
        gap, i = max(gaps)
        if gap < 0.5 * length:
            continue
        head, tail = ex[i + 1:], [(a + length, b + length) for a, b in ex[:i + 1]]
        t.exons = head + tail
        if t.cds:
            cds = sorted(t.cds)
            t.cds = [(a, b) for a, b in cds if a >= head[0][0]] + [(a + length, b + length) for a, b in cds if a < head[0][0]]
        t.start, t.end = t.exons[0][0], t.exons[-1][1]
        changed = True
    if changed:
        g.start = min(t.start for t in g.transcripts)
        g.end = max(t.end for t in g.transcripts)
        g.attrs["wraps_origin"] = "true"


def load_gff(path: str, only_seqid: str | None = None) -> Annotation:
    """Load GFF3 or GTF (auto-detected). Handles gene→mRNA→exon/CDS hierarchies and
    Ensembl/NCBI conventions; features without a gene parent become single-transcript genes."""
    ann = Annotation()
    ann.source = path
    genes: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}
    tx_parent: dict[str, str] = {}
    orphan_n = 0
    is_gtf = None
    pending_children: list[tuple[str, str, int, int]] = []   # (kind, parent_id, start, end)
    gene_parts: dict[str, list[tuple[int, int]]] = {}        # repeated gene IDs = pieces of one feature (GFF3 discontinuous)
    tx_parts: dict[str, list[tuple[int, int]]] = {}
    with _open(path) as fh:
        for ln in fh:
            if ln.startswith("##sequence-region"):
                bits = ln.split()
                if len(bits) >= 4 and bits[3].isdigit():
                    ann.lengths[bits[1]] = int(bits[3])
                continue
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, ftype, s, e, _score, strand, _phase, attr = parts[:9]
            if only_seqid and seqid != only_seqid:
                continue
            if ftype.lower() in ("region", "chromosome", "contig", "sequence_feature") and "Is_circular=true" in attr:
                ann.circular[seqid] = True
                ann.lengths.setdefault(seqid, int(e))
            if is_gtf is None:
                is_gtf = ("=" not in attr.split(";")[0]) and ('"' in attr)
            a = _attrs_gtf(attr) if is_gtf else _attrs_gff3(attr)
            start, end = int(s) - 1, int(e)
            st = -1 if strand == "-" else 1
            ftl = ftype.lower()
            if is_gtf:
                gid = a.get("gene_id", f"gene{start}")
                tid = a.get("transcript_id")
                gname = a.get("gene_name", a.get("gene", gid))
                if ftl == "gene" or gid not in genes:
                    g = genes.get(gid)
                    if g is None:
                        g = Gene(gid, gname, seqid, start, end, st, a.get("gene_biotype", a.get("gene_type", "")), attrs={})
                        genes[gid] = g
                        ann.add_gene(g)
                    if ftl == "gene":
                        g.start, g.end = start, end
                        continue
                if tid:
                    t = transcripts.get(tid)
                    if t is None:
                        t = Transcript(tid, a.get("transcript_name", tid), start, end, st,
                                       a.get("transcript_biotype", a.get("transcript_type", "")))
                        transcripts[tid] = t
                        genes[gid].transcripts.append(t)
                    if ftl == "exon":
                        t.exons.append((start, end))
                    elif ftl == "cds":
                        t.cds.append((start, end))
                    elif ftl in ("transcript", "mrna"):
                        t.start, t.end = start, end
                    g = genes[gid]
                    g.start, g.end = min(g.start, start), max(g.end, end)
                continue
            # ---- GFF3
            fid = a.get("ID", "")
            parent = a.get("Parent", "")
            if ftl in ("gene", "pseudogene", "ncrna_gene", "transposable_element_gene"):
                if fid and fid in genes:
                    # another piece of the same gene (e.g. an origin-spanning gene written as two rows)
                    g = genes[fid]
                    g.start, g.end = min(g.start, start), max(g.end, end)
                    gene_parts[fid].append((start, end))
                    continue
                name = a.get("Name", a.get("gene", a.get("gene_name", fid)))
                if fid:
                    gene_parts[fid] = [(start, end)]
                g = Gene(fid or name, name, seqid, start, end, st,
                         a.get("gene_biotype", a.get("biotype", a.get("gene_type", ftl))),
                         attrs={k: v for k, v in a.items() if k in ("description", "product", "gene", "Note", "Dbxref",
                                                                   "coverage", "sequence_ID", "partial_mapping", "low_identity", "extra_copy_number")})
                genes[g.id] = g
                ann.add_gene(g)
            elif ftl in ("mrna", "transcript", "ncrna", "lnc_rna", "rrna", "trna", "snorna", "snrna", "mirna",
                         "primary_transcript", "pseudogenic_transcript", "guide_rna", "scrna", "v_gene_segment",
                         "c_gene_segment", "d_gene_segment", "j_gene_segment", "lncrna", "antisense_rna", "rnase_p_rna", "srp_rna", "telomerase_rna", "y_rna", "vault_rna"):
                pid = parent.split(",")[0]
                if fid and fid in transcripts:
                    t = transcripts[fid]
                    t.start, t.end = min(t.start, start), max(t.end, end)
                    tx_parts[fid].append((start, end))
                    continue
                t = Transcript(fid or f"tx{start}", a.get("Name", a.get("transcript_id", fid)), start, end, st,
                               a.get("transcript_biotype", a.get("biotype", ftl)))
                if not fid:
                    fid = t.id
                tx_parts[fid] = [(start, end)]
                transcripts[fid] = t
                tx_parent[fid] = pid
                if pid in genes:
                    genes[pid].transcripts.append(t)
                else:
                    # transcript without a gene record: synthesise one
                    g = genes.get(pid) if pid else None
                    if g is None:
                        g = Gene(pid or fid, a.get("gene", a.get("Name", fid)), seqid, start, end, st, ftl)
                        genes[g.id] = g; ann.add_gene(g)
                    g.transcripts.append(t)
            elif ftl in ("exon", "cds", "five_prime_utr", "three_prime_utr"):
                for pid in parent.split(","):
                    t = transcripts.get(pid)
                    if t is None:
                        if pid in genes:
                            # exon directly under gene (e.g. AUGUSTUS-style): make an implicit transcript
                            g = genes[pid]
                            t = Transcript(pid + ".t", g.name, g.start, g.end, g.strand, "")
                            transcripts[pid] = t
                            g.transcripts.append(t)
                        else:
                            pending_children.append((ftl, pid, start, end))
                            continue
                    if ftl == "exon":
                        t.exons.append((start, end))
                    elif ftl == "cds":
                        t.cds.append((start, end))
    # children that arrived before their parents
    for ftl, pid, start, end in pending_children:
        t = transcripts.get(pid)
        if t is None:
            continue
        (t.exons if ftl == "exon" else t.cds if ftl == "cds" else []).append((start, end))
    # transcripts with CDS but no exons: use CDS as exons (or the pieces of a multi-row transcript)
    for tid, t in transcripts.items():
        if not t.exons and t.cds:
            t.exons = list(t.cds)
        elif not t.exons and len(tx_parts.get(tid, [])) > 1:
            t.exons = sorted(tx_parts[tid])
    # genes without transcripts: single exon = gene span (or the pieces of a multi-row gene)
    for g in genes.values():
        if not g.transcripts:
            parts = gene_parts.get(g.id) or [(g.start, g.end)]
            g.transcripts.append(Transcript(g.id + ".t", g.name, g.start, g.end, g.strand, g.biotype, sorted(parts)))
    ann.finalize()
    return ann


def load_bed(path: str) -> Annotation:
    ann = Annotation(); ann.source = path
    with _open(path) as fh:
        for k, ln in enumerate(fh):
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            seqid, s, e = p[0], int(p[1]), int(p[2])
            name = p[3] if len(p) > 3 else f"feat{k}"
            strand = -1 if len(p) > 5 and p[5] == "-" else 1
            t = Transcript(name, name, s, e, strand)
            if len(p) >= 12:
                ts, te = int(p[6]), int(p[7])
                sizes = [int(x) for x in p[10].rstrip(",").split(",")]
                starts = [int(x) for x in p[11].rstrip(",").split(",")]
                t.exons = [(s + st, s + st + sz) for st, sz in zip(starts, sizes)]
                if te > ts:
                    t.cds = [(max(a, ts), min(b, te)) for a, b in t.exons if b > ts and a < te]
            else:
                t.exons = [(s, e)]
            g = Gene(name, name, seqid, s, e, strand, "bed", [t])
            ann.add_gene(g)
    ann.finalize()
    return ann


def load_annotation(path: str, only_seqid: str | None = None) -> Annotation:
    low = path.lower().rstrip(".gz") if path.lower().endswith(".gz") else path.lower()
    if low.endswith(".bed"):
        return load_bed(path)
    if low.endswith((".gb", ".gbk", ".genbank")):
        return load_genbank(path)
    return load_gff(path, only_seqid)


def load_genbank(path: str) -> Annotation:
    """GenBank features -> genes. A bare `gene` feature is merged with the typed
    feature (CDS/tRNA/rRNA/…) that covers the same locus, so each gene appears once."""
    from Bio import SeqIO
    ann = Annotation(); ann.source = path
    with _open(path) as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            L = len(rec.seq)
            ann.lengths[rec.id] = L
            circ = str(rec.annotations.get("topology", "")).lower() == "circular"
            if circ:
                ann.circular[rec.id] = True
            typed: list[Gene] = []
            plain: list[Gene] = []
            for k, f in enumerate(rec.features):
                if f.type not in ("gene", "CDS", "rRNA", "tRNA", "ncRNA", "mRNA", "misc_feature", "D-loop"):
                    continue
                name = f.qualifiers.get("gene", f.qualifiers.get("product", f.qualifiers.get("locus_tag", [f.type])))[0]
                st = f.location.strand or 1
                attrs = {q: v[0] for q, v in f.qualifiers.items()
                         if q in ("product", "note", "transl_table", "db_xref")}
                # Biopython lists parts in transcription order (minus strand: descending);
                # put them in genomic order and unwrap a jump back past the origin, so a
                # feature such as the D-loop join(16024..16569,1..576) stays one feature
                # whose end lies beyond the sequence length instead of covering the genome.
                parts = [(int(p.start), int(p.end)) for p in f.location.parts]
                if st < 0:
                    parts = parts[::-1]
                span = max(e2 for _s2, e2 in parts) - min(s2 for s2, _e2 in parts)
                wrapped = False
                if len(parts) > 1 and (circ or span > 0.5 * L):
                    parts, wrapped = _unwrap_parts(parts, L)
                s2, e2 = min(a for a, _b in parts), max(b for _a, b in parts)
                t = Transcript(f"{name}.{k}", name, s2, e2, st, f.type, list(parts))
                if f.type == "CDS":
                    t.cds = list(t.exons)
                a2 = dict(attrs)
                if wrapped:
                    a2["wraps_origin"] = "true"
                g = Gene(f"{name}.{k}", name, rec.id, s2, e2, st, f.type, [t], attrs=a2)
                (plain if f.type == "gene" else typed).append(g)
            # keep typed features; add a plain gene only when no typed feature covers it
            for g in typed:
                ann.add_gene(g)
            for g in plain:
                if not any(t.name == g.name and t.start < g.end and t.end > g.start for t in typed):
                    ann.add_gene(g)
    ann.finalize()
    return ann


# ---------------------------------------------------------------------- PAF
def load_paf(path: str, min_len: int = 1000, min_mapq: int = 0, query: str | None = None) -> list[SyntenyBlock]:
    blocks = []
    with _open(path) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            if query and p[0] != query:
                continue
            alen = int(p[10])
            if alen < min_len or int(p[11]) < min_mapq:
                continue
            blocks.append(SyntenyBlock(p[0], int(p[2]), int(p[3]), -1 if p[4] == "-" else 1,
                                       p[5], int(p[7]), int(p[8]), int(p[9]), alen, int(p[11])))
    blocks.sort(key=lambda b: (b.qname, b.qstart))
    return blocks


def pack_lanes(items: Iterable[tuple[int, int, object]], gap: int = 0) -> list[list]:
    """Greedy interval packing: returns list of lanes, each a list of (start,end,obj)."""
    lanes: list[list] = []
    lane_end: list[int] = []
    for s, e, obj in sorted(items, key=lambda x: x[0]):
        for i, le in enumerate(lane_end):
            if s >= le + gap:
                lanes[i].append((s, e, obj)); lane_end[i] = e
                break
        else:
            lanes.append([(s, e, obj)]); lane_end.append(e)
    return lanes
