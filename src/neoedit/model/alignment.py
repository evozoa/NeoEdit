"""Pure-Python alignment model (no Qt) with undo/redo.

All edits go through methods that record an undoable snapshot of the rows they
touch, so undo/redo is uniform and cheap enough for interactive use.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable, Iterable, Optional

from Bio.Seq import Seq

GAP_CHARS = "-.~"
GAP = "-"

_GAP_TABLE = str.maketrans("", "", GAP_CHARS)
_DNA_RE = re.compile(r"^[ACGTUNRYKMSWBDHV\-\.\~\?]*$", re.I)


@dataclass
class SequenceRow:
    name: str
    seq: str
    description: str = ""
    id: str = ""

    def __len__(self) -> int:
        return len(self.seq)

    def ungapped(self) -> str:
        return self.seq.translate(_GAP_TABLE)

    def copy(self) -> "SequenceRow":
        return SequenceRow(self.name, self.seq, self.description, self.id)


@dataclass
class Feature:
    """A generic annotation on a row (ORF, primer, restriction site ...).

    Positions are 0-based, half-open, in *alignment* (gapped) coordinates.
    """
    row: int
    start: int
    end: int
    strand: int = 1           # +1 / -1
    type: str = "misc"
    label: str = ""
    color: str = "#3b82f6"
    data: dict = field(default_factory=dict)

    def __len__(self):
        return self.end - self.start


class _Edit:
    """Snapshot-based undo record."""

    def __init__(self, desc: str, before: dict, after: dict,
                 order_before: list[int] | None = None,
                 order_after: list[int] | None = None):
        self.desc = desc
        self.before = before      # {row_index: SequenceRow or None (row didn't exist)}
        self.after = after
        self.order_before = order_before
        self.order_after = order_after


class AlignmentModel:
    def __init__(self, rows: Iterable[SequenceRow] | None = None, seq_type: str | None = None):
        self.rows: list[SequenceRow] = list(rows or [])
        self.features: list[Feature] = []
        self._seq_type = seq_type
        self._detected_type = None
        self._undo: list[_Edit] = []
        self._redo: list[_Edit] = []
        self._listeners: list[Callable[[str], None]] = []
        self.path: Optional[str] = None
        self.format: Optional[str] = None
        self.dirty = False
        self._batch: _Edit | None = None

    # ------------------------------------------------------------------ basics
    @property
    def nrows(self) -> int:
        return len(self.rows)

    @property
    def width(self) -> int:
        return max((len(r.seq) for r in self.rows), default=0)

    @property
    def seq_type(self) -> str:
        if self._seq_type:
            return self._seq_type
        if self._detected_type is None:
            self._detected_type = self.detect_type()
        return self._detected_type

    @seq_type.setter
    def seq_type(self, v: str | None):
        self._seq_type = v
        self._detected_type = None
        self._emit("type")

    def detect_type(self) -> str:
        sample = "".join(r.seq[:3000].translate(_GAP_TABLE)[:2000] for r in self.rows[:20])
        if not sample:
            return "dna"
        if _DNA_RE.match(sample):
            return "rna" if "U" in sample.upper() and "T" not in sample.upper() else "dna"
        return "protein"

    def is_nucleotide(self) -> bool:
        return self.seq_type in ("dna", "rna")

    def residue(self, row: int, col: int) -> str:
        s = self.rows[row].seq
        return s[col] if 0 <= col < len(s) else ""

    def add_listener(self, fn: Callable[[str], None]):
        self._listeners.append(fn)

    def remove_listener(self, fn):
        if fn in self._listeners:
            self._listeners.remove(fn)

    def _emit(self, what: str = "data"):
        if what == "data":
            self._detected_type = None   # re-detected lazily on next access
        for fn in list(self._listeners):
            fn(what)

    # -------------------------------------------------------------------- undo
    def _record(self, desc: str, rows: Iterable[int], fn: Callable[[], None],
                order_change: bool = False):
        rows = sorted(set(rows))
        before = {i: (self.rows[i].copy() if i < len(self.rows) else None) for i in rows}
        order_before = list(range(len(self.rows))) if order_change else None
        fn()
        after = {i: (self.rows[i].copy() if i < len(self.rows) else None) for i in rows}
        edit = _Edit(desc, before, after)
        if self._batch is not None:
            # merge into batch: keep earliest 'before', latest 'after'
            for i, b in before.items():
                self._batch.before.setdefault(i, b)
            self._batch.after.update(after)
        else:
            self._undo.append(edit)
            self._redo.clear()
        self.dirty = True
        self._emit("data")

    def begin_batch(self, desc: str):
        if self._batch is None:
            self._batch = _Edit(desc, {}, {})

    def end_batch(self):
        if self._batch is not None:
            if self._batch.before:
                self._undo.append(self._batch)
                self._redo.clear()
            self._batch = None

    def _apply(self, snap: dict):
        for i in sorted(snap):
            r = snap[i]
            if r is None:
                if i < len(self.rows):
                    del self.rows[i]
            elif i < len(self.rows):
                self.rows[i] = r.copy()
            else:
                self.rows.append(r.copy())
        self.dirty = True
        self._emit("data")

    def can_undo(self) -> bool:
        return bool(self._undo)

    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo_text(self) -> str:
        return self._undo[-1].desc if self._undo else ""

    def redo_text(self) -> str:
        return self._redo[-1].desc if self._redo else ""

    def undo(self):
        if not self._undo:
            return
        e = self._undo.pop()
        # rows that were appended must be removed from the end first
        self._apply(e.before)
        self._redo.append(e)

    def redo(self):
        if not self._redo:
            return
        e = self._redo.pop()
        self._apply(e.after)
        self._undo.append(e)

    def clear_undo(self):
        self._undo.clear()
        self._redo.clear()

    # ------------------------------------------------------------ row editing
    def _pad(self, s: str, col: int) -> str:
        return s if len(s) >= col else s + GAP * (col - len(s))

    def insert_gaps(self, row: int, col: int, n: int = 1):
        if n <= 0:
            return
        def fn():
            r = self.rows[row]
            s = self._pad(r.seq, col)
            r.seq = s[:col] + GAP * n + s[col:]
        self._record("Insert gap", [row], fn)

    def delete_gaps(self, row: int, col: int, n: int = 1) -> int:
        """Delete up to n gap chars starting at col; stops at the first non-gap.
        Returns number deleted."""
        r = self.rows[row]
        k = 0
        while k < n and col + k < len(r.seq) and r.seq[col + k] in GAP_CHARS:
            k += 1
        if k == 0:
            return 0
        def fn():
            r.seq = r.seq[:col] + r.seq[col + k:]
        self._record("Delete gap", [row], fn)
        return k

    def delete_range(self, row: int, col: int, n: int = 1):
        """Delete n characters of any kind (residues included)."""
        r = self.rows[row]
        if col >= len(r.seq) or n <= 0:
            return
        def fn():
            r.seq = r.seq[:col] + r.seq[col + n:]
        self._record("Delete", [row], fn)

    def overwrite(self, row: int, col: int, text: str):
        def fn():
            r = self.rows[row]
            s = self._pad(r.seq, col)
            r.seq = s[:col] + text + s[col + len(text):]
        self._record("Overwrite", [row], fn)

    def insert_text(self, row: int, col: int, text: str):
        def fn():
            r = self.rows[row]
            s = self._pad(r.seq, col)
            r.seq = s[:col] + text + s[col:]
        self._record("Insert", [row], fn)

    def set_sequence(self, row: int, seq: str, desc: str = "Edit sequence"):
        def fn():
            self.rows[row].seq = seq
        self._record(desc, [row], fn)

    def rename(self, row: int, name: str):
        def fn():
            self.rows[row].name = name
        self._record("Rename", [row], fn)

    # ------------------------------------------------------- column editing
    def insert_gap_columns(self, col: int, n: int = 1, rows: Iterable[int] | None = None):
        idx = list(rows) if rows is not None else list(range(self.nrows))
        def fn():
            for i in idx:
                r = self.rows[i]
                s = self._pad(r.seq, col)
                r.seq = s[:col] + GAP * n + s[col:]
        self._record("Insert gap column", idx, fn)

    def delete_gap_columns(self, col: int, n: int = 1, rows: Iterable[int] | None = None) -> int:
        """Delete up to n gap-only columns at col across the given rows."""
        idx = list(rows) if rows is not None else list(range(self.nrows))
        k = 0
        while k < n:
            c = col + k
            if all((c >= len(self.rows[i].seq)) or self.rows[i].seq[c] in GAP_CHARS for i in idx) \
                    and any(c < len(self.rows[i].seq) for i in idx):
                k += 1
            else:
                break
        if k == 0:
            return 0
        def fn():
            for i in idx:
                r = self.rows[i]
                r.seq = r.seq[:col] + r.seq[col + k:]
        self._record("Delete gap column", idx, fn)
        return k

    def remove_gap_only_columns(self):
        w = self.width
        keep = [c for c in range(w)
                if any(c < len(r.seq) and r.seq[c] not in GAP_CHARS for r in self.rows)]
        if len(keep) == w:
            return
        def fn():
            for r in self.rows:
                r.seq = "".join(r.seq[c] for c in keep if c < len(r.seq))
        self._record("Remove gap-only columns", range(self.nrows), fn)

    def pad_to_equal_length(self):
        w = self.width
        def fn():
            for r in self.rows:
                r.seq = self._pad(r.seq, w)
        self._record("Pad sequences", range(self.nrows), fn)

    def block_shift(self, rows: Iterable[int], start: int, end: int, delta: int) -> bool:
        """Slide the block [start,end) in given rows by delta columns, consuming
        gaps on the destination side and leaving gaps behind (BioEdit drag-move).
        Returns False if a non-gap would have to be overwritten."""
        idx = list(rows)
        if delta == 0:
            return True
        new = {}
        for i in idx:
            s = self._pad(self.rows[i].seq, end)
            block = s[start:end]
            if delta > 0:
                region = s[end:end + delta]
                if len(region) < delta:
                    region = region + GAP * (delta - len(region))
                    s = s + GAP * (delta - len(s[end:end + delta]))
                if any(c not in GAP_CHARS for c in region):
                    return False
                ns = s[:start] + GAP * delta + block + s[end + delta:]
            else:
                d = -delta
                if start - d < 0:
                    return False
                region = s[start - d:start]
                if any(c not in GAP_CHARS for c in region):
                    return False
                ns = s[:start - d] + block + GAP * d + s[end:]
            new[i] = ns
        def fn():
            for i, ns in new.items():
                self.rows[i].seq = ns
        self._record("Move block", idx, fn)
        return True

    # ----------------------------------------------------- sequence ops
    def _transform(self, rows: Iterable[int], fn_seq: Callable[[str], str], desc: str):
        idx = list(rows)
        def fn():
            for i in idx:
                self.rows[i].seq = fn_seq(self.rows[i].seq)
        self._record(desc, idx, fn)

    def reverse_complement(self, rows: Iterable[int]):
        self._transform(rows, lambda s: str(Seq(s).reverse_complement()), "Reverse complement")

    def complement(self, rows: Iterable[int]):
        self._transform(rows, lambda s: str(Seq(s).complement()), "Complement")

    def reverse(self, rows: Iterable[int]):
        self._transform(rows, lambda s: s[::-1], "Reverse")

    def to_upper(self, rows: Iterable[int]):
        self._transform(rows, str.upper, "Uppercase")

    def to_lower(self, rows: Iterable[int]):
        self._transform(rows, str.lower, "Lowercase")

    def remove_gaps(self, rows: Iterable[int]):
        self._transform(rows, lambda s: "".join(c for c in s if c not in GAP_CHARS), "Remove gaps")

    def dna_to_rna(self, rows: Iterable[int]):
        self._transform(rows, lambda s: s.replace("T", "U").replace("t", "u"), "DNA -> RNA")

    def rna_to_dna(self, rows: Iterable[int]):
        self._transform(rows, lambda s: s.replace("U", "T").replace("u", "t"), "RNA -> DNA")

    # ------------------------------------------------------------ rows
    def add_row(self, row: SequenceRow, at: int | None = None):
        pos = self.nrows if at is None else at
        affected = range(pos, self.nrows + 1)
        def fn():
            self.rows.insert(pos, row)
        self._record("Add sequence", affected, fn)

    def remove_rows(self, rows: Iterable[int]):
        idx = sorted(set(rows))
        if not idx:
            return
        affected = range(idx[0], self.nrows)
        def fn():
            for i in reversed(idx):
                del self.rows[i]
        self._record("Delete sequence", affected, fn)
        self.features = [f for f in self.features if f.row not in idx]

    def move_rows(self, rows: Iterable[int], delta: int):
        idx = sorted(set(rows))
        if not idx:
            return
        if delta < 0 and idx[0] + delta < 0:
            return
        if delta > 0 and idx[-1] + delta >= self.nrows:
            return
        lo = min(idx[0], idx[0] + delta)
        hi = max(idx[-1], idx[-1] + delta)
        def fn():
            moving = [self.rows[i] for i in idx]
            rest = [r for i, r in enumerate(self.rows) if i not in idx]
            insert_at = idx[0] + delta
            self.rows = rest[:insert_at] + moving + rest[insert_at:]
        self._record("Move sequence", range(lo, hi + 1), fn)

    def duplicate_rows(self, rows: Iterable[int]):
        idx = sorted(set(rows))
        for k, i in enumerate(idx):
            r = self.rows[i + k].copy()
            r.name = r.name + "_copy"
            self.add_row(r, at=i + k + 1)

    # ------------------------------------------------------------ misc
    def column_slice(self, start: int, end: int) -> "AlignmentModel":
        m = AlignmentModel([SequenceRow(r.name, self._pad(r.seq, end)[start:end], r.description, r.id)
                            for r in self.rows], self._seq_type)
        return m

    def find(self, pattern: str, start_row: int = 0, start_col: int = 0,
             ignore_gaps: bool = True, regex: bool = False, case: bool = False):
        """Search rows from (start_row,start_col) onward, wrapping. Returns
        (row, start, end) in alignment coords or None."""
        flags = 0 if case else re.I
        if regex:
            pat = pattern
        elif ignore_gaps:
            pat = "[-.~]*".join(re.escape(c) for c in pattern)
        else:
            pat = re.escape(pattern)
        pattern = pat
        rx = re.compile(pattern, flags)
        order = list(range(start_row, self.nrows)) + list(range(0, start_row + 1))
        for n, i in enumerate(order):
            s = self.rows[i].seq
            off = start_col if n == 0 else 0
            if n == len(order) - 1:
                # wrapped back to the starting row: search only before start_col
                m = rx.search(s[:start_col] if start_col else "")
                if m:
                    return i, m.start(), m.end()
                continue
            m = rx.search(s, off)
            if m:
                return i, m.start(), m.end()
        return None
