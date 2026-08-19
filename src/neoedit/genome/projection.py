"""Gap-aware mapping between reference coordinates and alignment columns.

The reference is one (possibly gapped) row of the alignment. Memory-light:
stores only the sorted positions of gap columns, so an ungapped 33-Mb
chromosome costs nothing and a gapped 16-kb mitogenome costs a few hundred ints.
"""
from __future__ import annotations

import bisect

GAPSET = set("-.~")


class RefProjection:
    def __init__(self, ref_seq: str):
        self.ncols = len(ref_seq)
        self.gap_cols = [i for i, c in enumerate(ref_seq) if c in GAPSET]
        self.ref_len = self.ncols - len(self.gap_cols)

    @property
    def identity(self) -> bool:
        return not self.gap_cols

    def col_to_ref(self, col: int) -> int:
        """Alignment column -> reference position (position of the residue at or
        left of this column). Clamped to [0, ref_len]."""
        col = max(0, min(self.ncols, col))
        return col - bisect.bisect_right(self.gap_cols, col - 1)

    def ref_to_col(self, u: int) -> int:
        """Reference position -> alignment column of that residue (u == ref_len ->
        one past the last residue's column)."""
        if not self.gap_cols:
            return max(0, min(self.ncols, u))
        u = max(0, min(self.ref_len, u))
        col = u
        while True:
            k = bisect.bisect_right(self.gap_cols, col)
            ncol = u + k
            if ncol == col:
                return col
            col = ncol

    def span_to_cols(self, s: int, e: int) -> tuple[int, int]:
        """Reference span [s,e) -> alignment column span [c0,c1)."""
        if e <= s:
            c = self.ref_to_col(s)
            return c, c
        return self.ref_to_col(s), self.ref_to_col(e - 1) + 1
