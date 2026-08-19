"""Random access to large FASTA files via a samtools-style .fai index (pure Python)."""
from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass
class FaiRecord:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    """Open a (possibly multi-gigabyte) FASTA and fetch sub-sequences by coordinate.

    Uses `<path>.fai` if present, otherwise builds one (single pass) and writes it
    next to the FASTA when the directory is writable.
    """

    def __init__(self, path: str):
        self.path = path
        self.records: dict[str, FaiRecord] = {}
        self.order: list[str] = []
        fai = path + ".fai"
        if os.path.exists(fai) and os.path.getmtime(fai) >= os.path.getmtime(path):
            self._read_fai(fai)
        else:
            self._build_index()
            try:
                self._write_fai(fai)
            except OSError:
                pass
        self._fh = open(path, "rb")

    # ------------------------------------------------------------------ index
    def _read_fai(self, fai: str):
        with open(fai) as fh:
            for ln in fh:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                rec = FaiRecord(parts[0], int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
                self.records[rec.name] = rec
                self.order.append(rec.name)

    def _build_index(self):
        with open(self.path, "rb") as fh:
            name = None; length = 0; offset = 0; line_bases = 0; line_width = 0
            pos = 0
            for raw in fh:
                if raw.startswith(b">"):
                    if name is not None:
                        rec = FaiRecord(name, length, offset, line_bases, line_width)
                        self.records[name] = rec; self.order.append(name)
                    name = raw[1:].split()[0].decode() if raw[1:].split() else f"seq{len(self.order) + 1}"
                    pos += len(raw)
                    offset = pos; length = 0; line_bases = 0; line_width = 0
                    continue
                if name is not None:
                    bases = len(raw.rstrip(b"\r\n"))
                    if line_bases == 0:
                        line_bases, line_width = bases, len(raw)
                    length += bases
                pos += len(raw)
            if name is not None:
                rec = FaiRecord(name, length, offset, line_bases, line_width)
                self.records[name] = rec; self.order.append(name)

    def _write_fai(self, fai: str):
        with open(fai, "w") as fh:
            for n in self.order:
                r = self.records[n]
                fh.write(f"{r.name}\t{r.length}\t{r.offset}\t{r.line_bases}\t{r.line_width}\n")

    # ------------------------------------------------------------------ access
    def names(self) -> list[str]:
        return list(self.order)

    def length(self, name: str) -> int:
        return self.records[name].length

    def contigs(self) -> list[tuple[str, int]]:
        return [(n, self.records[n].length) for n in self.order]

    def fetch(self, name: str, start: int = 0, end: int | None = None) -> str:
        """0-based, half-open. Clamps to contig bounds."""
        r = self.records[name]
        start = max(0, start)
        end = r.length if end is None else min(end, r.length)
        if end <= start:
            return ""
        if r.line_bases == 0:
            return ""
        lines_before = start // r.line_bases
        first = r.offset + lines_before * r.line_width + (start % r.line_bases)
        lines_after = (end - 1) // r.line_bases
        last = r.offset + lines_after * r.line_width + ((end - 1) % r.line_bases) + 1
        self._fh.seek(first)
        raw = self._fh.read(last - first)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode()

    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass
