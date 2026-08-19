"""Restriction enzyme site search (REBASE via Biopython), alignment-aware.

Beyond a plain map, sites can be summarised across all sequences of an alignment:
which enzymes cut every sequence, only some (diagnostic), or none.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

from Bio.Restriction import AllEnzymes, CommOnly, RestrictionBatch
from Bio.Seq import Seq

from ..model.alignment import GAP_CHARS, Feature

GAPSET = set(GAP_CHARS)


@dataclass
class EnzymeHit:
    enzyme: str
    site: str
    elucidate: str
    size: int
    overhang: str            # "blunt" | "5'" | "3'"
    suppliers: list[str]
    positions: list[int]     # 1-based positions of the base AFTER the top-strand cut (Biopython convention)
    fst5: int = 1            # cut offset within the recognition site (Biopython)
    row: int = 0

    def site_start(self, cut: int) -> int:
        """0-based start of the recognition site for a given cut position."""
        return cut - self.fst5 - 1

    @property
    def n_cuts(self) -> int:
        return len(self.positions)


@dataclass
class AlignmentSummary:
    enzyme: str
    site: str
    overhang: str
    size: int
    cuts_per_row: dict[int, list[int]] = field(default_factory=dict)

    def rows_cut(self) -> list[int]:
        return [r for r, p in self.cuts_per_row.items() if p]

    def is_universal(self, rows: Sequence[int]) -> bool:
        return all(self.cuts_per_row.get(r) for r in rows)

    def is_diagnostic(self, rows: Sequence[int]) -> bool:
        cut = [bool(self.cuts_per_row.get(r)) for r in rows]
        return any(cut) and not all(cut)


def enzyme_pool(commercial_only: bool = True, suppliers: Sequence[str] | None = None,
                min_site: int = 4, max_site: int = 8, blunt: bool | None = None,
                names: Sequence[str] | None = None) -> RestrictionBatch:
    """Build a RestrictionBatch honouring the usual filters."""
    if names:
        wanted = {n.strip().lower() for n in names if n.strip()}
        pool = [e for e in AllEnzymes if str(e).lower() in wanted]
    else:
        pool = list(CommOnly if commercial_only else AllEnzymes)
        if suppliers:
            want = {s.lower() for s in suppliers}
            pool = [e for e in pool if any(s.lower() in want for s in e.supplier_list())]
        pool = [e for e in pool if min_site <= e.size <= max_site]
        if blunt is True:
            pool = [e for e in pool if e.is_blunt()]
        elif blunt is False:
            pool = [e for e in pool if not e.is_blunt()]
    return RestrictionBatch(pool)


def all_suppliers() -> list[str]:
    out = set()
    for e in CommOnly:
        out.update(e.supplier_list())
    return sorted(out)


def _overhang(e) -> str:
    if e.is_blunt():
        return "blunt"
    return "5'" if e.is_5overhang() else "3'"


def search_sequence(seq: str, batch: RestrictionBatch, linear: bool = True,
                    row: int = 0, min_cuts: int = 1, max_cuts: int | None = None) -> list[EnzymeHit]:
    """Sites of `batch` in `seq` (gaps ignored). min_cuts=0 also reports non-cutters."""
    clean = "".join(c for c in seq if c not in GAPSET).upper().replace("U", "T")
    if not clean:
        return []
    res = batch.search(Seq(clean), linear=linear)
    hits = []
    for enz, pos in res.items():
        n = len(pos)
        if n < min_cuts or (max_cuts is not None and n > max_cuts):
            continue
        hits.append(EnzymeHit(str(enz), enz.site, enz.elucidate(), enz.size, _overhang(enz),
                              enz.supplier_list(), list(pos), getattr(enz, "fst5", 1), row))
    hits.sort(key=lambda h: (h.n_cuts, h.enzyme))
    return hits


def search_alignment(rows: Sequence[str], batch: RestrictionBatch, linear: bool = True,
                     row_indices: Sequence[int] | None = None) -> list[AlignmentSummary]:
    idx = list(row_indices) if row_indices is not None else list(range(len(rows)))
    summaries: dict[str, AlignmentSummary] = {}
    for r in idx:
        for h in search_sequence(rows[r], batch, linear, r):
            s = summaries.get(h.enzyme)
            if s is None:
                s = summaries[h.enzyme] = AlignmentSummary(h.enzyme, h.site, h.overhang, h.size)
            s.cuts_per_row[r] = h.positions
    for r in idx:
        for s in summaries.values():
            s.cuts_per_row.setdefault(r, [])
    return sorted(summaries.values(), key=lambda s: s.enzyme)


def ungapped_to_columns(seq: str) -> list[int]:
    return [i for i, c in enumerate(seq) if c not in GAPSET]


def hits_to_features(hits: Iterable[EnzymeHit], seq_by_row: dict[int, str], color: str = "#0ea5e9") -> list[Feature]:
    """Sites as features in alignment coordinates (site span, not the cut point)."""
    feats = []
    maps: dict[int, list[int]] = {}
    for h in hits:
        gm = maps.get(h.row)
        if gm is None:
            gm = maps[h.row] = ungapped_to_columns(seq_by_row.get(h.row, ""))
        for cut in h.positions:
            # place the feature over the recognition site (cut position minus the cut offset)
            start_u = max(0, h.site_start(cut))
            end_u = min(len(gm), start_u + h.size)
            if end_u <= start_u or start_u >= len(gm):
                continue
            c0 = gm[start_u]
            c1 = gm[min(end_u, len(gm)) - 1] + 1
            feats.append(Feature(h.row, c0, c1, 1, "restriction_site", f"{h.enzyme} ({h.site})", color,
                                 data={"enzyme": h.enzyme, "cut": cut, "overhang": h.overhang}))
    return feats


def digest_fragments(seq: str, enzymes: Sequence[str], linear: bool = True) -> list[tuple[int, int, int]]:
    """(start, end, length) of fragments produced by digesting with the given enzymes."""
    clean = "".join(c for c in seq if c not in GAPSET).upper()
    batch = enzyme_pool(names=enzymes)
    cuts = sorted({c for pos in batch.search(Seq(clean), linear=linear).values() for c in pos})
    if not cuts:
        return [(0, len(clean), len(clean))]
    frags = []
    if linear:
        bounds = [0] + [c - 1 for c in cuts] + [len(clean)]
        for a, b in zip(bounds, bounds[1:]):
            frags.append((a, b, b - a))
    else:
        cs = [c - 1 for c in cuts]
        for a, b in zip(cs, cs[1:] + [cs[0] + len(clean)]):
            frags.append((a, b % len(clean) if b != cs[0] + len(clean) else cs[0], (b - a) % len(clean) or len(clean)))
    return frags
