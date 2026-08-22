"""ORF finder with alternate genetic codes, MitoFinder-style options.

Positions are 0-based half-open on the *ungapped* sequence; `map_to_gapped`
converts back to alignment coordinates.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from Bio.Data import CodonTable
from Bio.Seq import Seq

from ..model.alignment import GAP_CHARS, Feature

GAPSET = set(GAP_CHARS)


@dataclass
class ORF:
    start: int            # 0-based, ungapped, on forward strand coordinates
    end: int              # exclusive
    strand: int           # +1/-1
    frame: int            # 1..3 (relative to the strand's 5' end)
    table: int
    start_codon: str
    stop_codon: str       # "" if partial at 3'
    partial5: bool
    partial3: bool
    nt: str
    aa: str
    row: int = 0
    name: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def length_nt(self) -> int:
        return self.end - self.start

    @property
    def length_aa(self) -> int:
        return len(self.aa.rstrip("*"))

    @property
    def wraps(self) -> bool:
        """True when the ORF continues past the origin of a circular molecule (end > length)."""
        return bool(self.extra.get("length")) and self.end > self.extra["length"]

    def label(self) -> str:
        return self.name or f"ORF {self.length_aa}aa {'+' if self.strand > 0 else '-'}{self.frame}"

    def _data(self, **extra) -> dict:
        d = {"aa": self.aa, "nt": self.nt, "table": self.table, "partial5": self.partial5, "partial3": self.partial3}
        d.update(extra)
        return d

    def to_feature(self, gap_map: list[int] | None = None, color="#f59e0b") -> Feature:
        """Single Feature in alignment columns (an ORF across the origin is clamped; see to_features)."""
        L = self.extra.get("length")
        s, e = self.start, min(self.end, L) if L else self.end
        if gap_map:
            s = gap_map[s] if s < len(gap_map) else len(gap_map)
            e = gap_map[e - 1] + 1 if e - 1 < len(gap_map) else len(gap_map)
        return Feature(self.row, s, e, self.strand, "ORF", self.label(), color, data=self._data())

    def to_features(self, gap_map: list[int] | None = None, color="#f59e0b") -> list[Feature]:
        """Feature piece(s): one normally; two for an ORF across the origin, each carrying the
        `phase` (bases to skip at its 5' end) so the amino-acid line translates in frame."""
        L = self.extra.get("length")
        if not L or self.end <= L:
            return [self.to_feature(gap_map, color)]
        def cols(a, b):
            if gap_map:
                return (gap_map[a] if a < len(gap_map) else len(gap_map),
                        gap_map[b - 1] + 1 if b - 1 < len(gap_map) else len(gap_map))
            return a, b
        p1, p2 = (self.start, L), (0, self.end - L)
        n1, n2 = L - self.start, self.end - L
        if self.strand > 0:
            ph1, ph2 = 0, (3 - n1 % 3) % 3
        else:
            ph1, ph2 = (3 - n2 % 3) % 3, 0
        lab = self.label()
        out = []
        for (a, b), ph, part in ((p1, ph1, 1), (p2, ph2, 2)):
            ca, cb = cols(a, b)
            out.append(Feature(self.row, ca, cb, self.strand, "ORF", lab, color,
                               data=self._data(phase=ph, wrap_part=part, wrap_len=self.end - self.start)))
        return out


def gap_map(seq: str) -> list[int]:
    """Index i in ungapped -> index in gapped sequence."""
    return [i for i, c in enumerate(seq) if c not in GAPSET]


def find_orfs(seq: str, table: int = 1, min_aa: int = 30,
              start_mode: str = "table",       # "atg" | "table" | "any"
              both_strands: bool = True,
              allow_partial: bool = True,       # ORFs running off the sequence ends
              nested: bool = False,             # also report ORFs starting inside a longer one (same frame)
              circular: bool = False,
              row: int = 0) -> list[ORF]:
    """Scan a sequence for ORFs.

    start_mode:
      "atg"   - only ATG starts
      "table" - any start codon of the genetic code (e.g. ATA/ATT/GTG for mito)
      "any"   - open reading frame from first codon after a stop (stop-to-stop)
    allow_partial:
      report frames that begin at the sequence 5' end without a start (partial5)
      and/or reach the 3' end without a stop (partial3). Mito genes frequently
      end with incomplete stop codons (T/TA) that get polyadenylated, so 3'
      partial ORFs are important.
    """
    tbl = CodonTable.unambiguous_dna_by_id[table]
    stops = set(tbl.stop_codons)
    starts = {"ATG"} if start_mode == "atg" else set(tbl.start_codons) | {"ATG"}
    fwd = tbl.forward_table

    clean = "".join(c for c in seq if c not in GAPSET).upper().replace("U", "T")
    n = len(clean)
    if circular and n >= 3:
        # scan the doubled molecule so frames run straight through the origin; there are
        # no molecule ends, so edge-partial ORFs make no sense here
        out = find_orfs(clean + clean, table, min_aa, start_mode, both_strands, allow_partial=False,
                        nested=nested, circular=False, row=row)
        kept: list[ORF] = []
        seen: set[tuple[int, int, int]] = set()
        for o in out:
            if o.end - o.start > n:                   # longer than the molecule: a frame with no stop at all
                continue
            st = o.start % n
            key = (st, o.end - o.start, o.strand)
            if key in seen:
                continue
            seen.add(key)
            o.start, o.end = st, st + (o.end - o.start)
            o.extra["length"] = n
            kept.append(o)
        kept.sort(key=lambda o: (-o.length_aa, o.start))
        return kept
    out: list[ORF] = []
    strands = [(+1, clean)]
    if both_strands:
        strands.append((-1, str(Seq(clean).reverse_complement())))

    def translate(s):
        aa = []
        for i in range(0, len(s) - 2, 3):
            c = s[i:i + 3]
            aa.append("*" if c in stops else fwd.get(c, "X"))
        return "".join(aa)

    for strand, s in strands:
        L = len(s)
        for frame in range(3):
            i = frame
            region_start = None   # start of current candidate (index of start codon)
            in_orf_from_edge = allow_partial and start_mode != "any"
            # candidate starts within a stop-to-stop segment (for nested reporting)
            seg_starts: list[int] = []
            seg_begin = frame
            while i + 3 <= L:
                codon = s[i:i + 3]
                is_stop = codon in stops
                if start_mode == "any":
                    # stop-to-stop: segment from seg_begin to stop
                    if is_stop:
                        _emit(out, s, seg_begin, i + 3, strand, frame, table, stops, fwd, min_aa,
                              partial5=(seg_begin == frame), row=row)
                        seg_begin = i + 3
                else:
                    if is_stop:
                        if seg_starts:
                            cands = seg_starts if nested else seg_starts[:1]
                            for st in cands:
                                _emit(out, s, st, i + 3, strand, frame, table, stops, fwd, min_aa, row=row)
                        elif in_orf_from_edge and allow_partial:
                            # 5' partial: from frame start to this stop
                            _emit(out, s, seg_begin, i + 3, strand, frame, table, stops, fwd, min_aa,
                                  partial5=True, row=row)
                        seg_starts = []
                        seg_begin = i + 3
                        in_orf_from_edge = False
                    elif codon in starts:
                        seg_starts.append(i)
                i += 3
            # ran off the 3' end
            if allow_partial:
                if start_mode == "any":
                    _emit(out, s, seg_begin, L - ((L - frame) % 3), strand, frame, table, stops, fwd, min_aa,
                          partial5=(seg_begin == frame), partial3=True, row=row)
                elif seg_starts:
                    cands = seg_starts if nested else seg_starts[:1]
                    for st in cands:
                        _emit(out, s, st, L - ((L - st) % 3), strand, frame, table, stops, fwd, min_aa,
                              partial3=True, row=row)
                elif in_orf_from_edge:
                    _emit(out, s, seg_begin, L - ((L - frame) % 3), strand, frame, table, stops, fwd, min_aa,
                          partial5=True, partial3=True, row=row)

    # convert minus-strand coordinates back to forward coordinates
    for o in out:
        if o.strand < 0:
            o.start, o.end = n - o.end, n - o.start
    out.sort(key=lambda o: (-o.length_aa, o.start))
    return out


def _emit(out, s, st, en, strand, frame, table, stops, fwd, min_aa, partial5=False, partial3=False, row=0):
    if en - st < 3:
        return
    nt = s[st:en]
    aa = []
    for i in range(0, len(nt) - 2, 3):
        c = nt[i:i + 3]
        aa.append("*" if c in stops else fwd.get(c, "X"))
    aa = "".join(aa)
    if len(aa.rstrip("*")) < min_aa:
        return
    stop_codon = nt[-3:] if aa.endswith("*") else ""
    out.append(ORF(st, en, strand, frame + 1, table, nt[:3], stop_codon, partial5, partial3, nt, aa, row=row))


def orfs_to_gff(orfs: Iterable[ORF], seqid: str) -> str:
    lines = ["##gff-version 3"]
    for k, o in enumerate(orfs, 1):
        attrs = f"ID=orf{k};Name={o.name or 'ORF' + str(k)};length_aa={o.length_aa};table={o.table}"
        if o.partial5:
            attrs += ";partial5=true"
        if o.partial3:
            attrs += ";partial3=true"
        L = o.extra.get("length")
        if L and o.end > L:     # across the origin: two rows sharing the ID, as NCBI does
            pieces = [(o.start, L), (0, o.end - L)]
            if o.strand < 0:
                pieces.reverse()
            for a, b in pieces:
                lines.append("\t".join([seqid, "neoedit", "CDS", str(a + 1), str(b), ".", "+" if o.strand > 0 else "-", "0", attrs]))
            continue
        lines.append("\t".join([seqid, "neoedit", "CDS", str(o.start + 1), str(o.end),
                                ".", "+" if o.strand > 0 else "-", "0", attrs]))
    return "\n".join(lines) + "\n"


def orfs_to_genbank_table(orfs: Iterable[ORF], seqid: str) -> str:
    """NCBI 5-column feature table."""
    lines = [f">Feature {seqid}"]
    for k, o in enumerate(orfs, 1):
        L = o.extra.get("length")
        if L and o.end > L:     # across the origin: two location lines (5-column table convention)
            if o.strand > 0:
                lines.append(f"{o.start + 1}\t{L}\tCDS"); lines.append(f"1\t{o.end - L}")
            else:
                lines.append(f"{o.end - L}\t1\tCDS"); lines.append(f"{L}\t{o.start + 1}")
            lines.append(f"\t\t\tproduct\t{o.name or 'hypothetical protein'}")
            lines.append(f"\t\t\ttransl_table\t{o.table}")
            continue
        a, b = (o.start + 1, o.end) if o.strand > 0 else (o.end, o.start + 1)
        a_s = f"<{a}" if (o.partial5 and o.strand > 0) or (o.partial3 and o.strand < 0) else str(a)
        b_s = f">{b}" if (o.partial3 and o.strand > 0) or (o.partial5 and o.strand < 0) else str(b)
        lines.append(f"{a_s}\t{b_s}\tCDS")
        lines.append(f"\t\t\tproduct\t{o.name or 'hypothetical protein'}")
        lines.append(f"\t\t\ttransl_table\t{o.table}")
    return "\n".join(lines) + "\n"


def orfs_to_fasta(orfs: Iterable[ORF], seqid: str, protein: bool = True) -> str:
    out = []
    for k, o in enumerate(orfs, 1):
        hdr = f">{seqid}_ORF{k} {o.start + 1}-{o.end} strand={'+' if o.strand > 0 else '-'} frame={o.frame} len={o.length_aa}aa table={o.table}"
        out.append(hdr)
        out.append(o.aa.rstrip("*") if protein else o.nt)
    return "\n".join(out) + "\n"
