"""Translation, consensus, identity and conservation utilities (no Qt)."""
from __future__ import annotations

import math
from collections import Counter
from typing import Iterable

import numpy as np
from Bio.Data import CodonTable
from Bio.Seq import Seq

from ..model.alignment import GAP_CHARS

GAPSET = set(GAP_CHARS)


def codon_tables() -> list[tuple[int, str]]:
    """(id, name) for all NCBI genetic codes, sorted by id."""
    out = []
    for tid, tbl in CodonTable.unambiguous_dna_by_id.items():
        out.append((tid, tbl.names[0]))
    return sorted(out)


def translate_gapped(seq: str, table: int = 1, frame: int = 0, to_stop: bool = False) -> str:
    """Translate a (possibly gapped) DNA/RNA string in frame 0..2.

    Gaps are removed first; codons containing ambiguity translate to X.
    """
    s = "".join(c for c in seq if c not in GAPSET).upper().replace("U", "T")
    s = s[frame:]
    s = s[: len(s) - (len(s) % 3)]
    if not s:
        return ""
    try:
        return str(Seq(s).translate(table=table, to_stop=to_stop))
    except Exception:
        tbl = CodonTable.unambiguous_dna_by_id[table].forward_table
        stops = set(CodonTable.unambiguous_dna_by_id[table].stop_codons)
        out = []
        for i in range(0, len(s), 3):
            c = s[i:i + 3]
            if c in stops:
                if to_stop:
                    break
                out.append("*")
            else:
                out.append(tbl.get(c, "X"))
        return "".join(out)


def translate_aligned(seq: str, table: int = 1, frame: int = 0) -> str:
    """Translate keeping alignment: output has one char per codon *column triple*
    so it can be displayed beneath the DNA (each aa followed by two spaces).
    Codon gaps (---) produce '-'. Partial-gap codons produce 'X'."""
    s = seq[frame:]
    out = []
    tbl = CodonTable.unambiguous_dna_by_id[table]
    fwd, stops = tbl.forward_table, set(tbl.stop_codons)
    for i in range(0, len(s) - 2, 3):
        c = s[i:i + 3].upper().replace("U", "T")
        if all(ch in GAPSET for ch in c):
            out.append("-")
        elif any(ch in GAPSET for ch in c):
            out.append("X")
        elif c in stops:
            out.append("*")
        else:
            out.append(fwd.get(c, "X"))
    return "".join(out)


def six_frame(seq: str, table: int = 1) -> dict[str, str]:
    s = "".join(c for c in seq if c not in GAPSET)
    rc = str(Seq(s.upper().replace("U", "T")).reverse_complement())
    res = {}
    for f in range(3):
        res[f"+{f + 1}"] = translate_gapped(s, table, f)
        res[f"-{f + 1}"] = translate_gapped(rc, table, f)
    return res


def consensus(rows: Iterable[str], threshold: float = 0.5, ignore_gaps: bool = True,
              plurality: bool = True) -> str:
    rows = list(rows)
    if not rows:
        return ""
    w = max(len(r) for r in rows)
    out = []
    for c in range(w):
        col = [r[c].upper() for r in rows if c < len(r)]
        if ignore_gaps:
            res = [x for x in col if x not in GAPSET]
        else:
            res = col
        if not res:
            out.append("-")
            continue
        ch, n = Counter(res).most_common(1)[0]
        frac = n / (len(res) if ignore_gaps else len(col))
        if frac >= threshold or (plurality and frac >= 0.5):
            out.append(ch)
        else:
            out.append("N" if _looks_nt(res) else "X")
    return "".join(out)


def _looks_nt(chars):
    return all(c in "ACGTUNRYKMSWBDHV" for c in chars)


def identity_matrix(rows: list[str], ignore_gap_pairs: bool = True) -> np.ndarray:
    n = len(rows)
    mat = np.zeros((n, n))
    for i in range(n):
        a = rows[i].upper()
        for j in range(i, n):
            b = rows[j].upper()
            L = min(len(a), len(b))
            same = tot = 0
            for k in range(L):
                x, y = a[k], b[k]
                if x in GAPSET and y in GAPSET:
                    continue
                if ignore_gap_pairs and (x in GAPSET or y in GAPSET):
                    continue
                tot += 1
                if x == y:
                    same += 1
            v = same / tot if tot else 0.0
            mat[i, j] = mat[j, i] = v
    return mat


def column_entropy(rows: list[str]) -> np.ndarray:
    """Shannon entropy (bits) per column, gaps excluded."""
    w = max((len(r) for r in rows), default=0)
    ent = np.zeros(w)
    for c in range(w):
        col = [r[c].upper() for r in rows if c < len(r) and r[c] not in GAPSET]
        if not col:
            continue
        cnt = Counter(col)
        n = len(col)
        ent[c] = -sum((v / n) * math.log2(v / n) for v in cnt.values())
    return ent


def column_identity(rows: list[str]) -> np.ndarray:
    """Fraction of (non-gap) residues matching the plurality residue per column."""
    w = max((len(r) for r in rows), default=0)
    out = np.zeros(w)
    for c in range(w):
        col = [r[c].upper() for r in rows if c < len(r) and r[c] not in GAPSET]
        if col:
            out[c] = Counter(col).most_common(1)[0][1] / len(col)
    return out


def gc_content(seq: str) -> float:
    s = [c for c in seq.upper() if c not in GAPSET]
    if not s:
        return 0.0
    return sum(1 for c in s if c in "GCS") / len(s)


def composition(seq: str) -> dict[str, int]:
    return dict(Counter(c for c in seq.upper() if c not in GAPSET))


def molecular_weight(seq: str, seq_type: str) -> float | None:
    from Bio.SeqUtils import molecular_weight as mw
    s = "".join(c for c in seq.upper() if c not in GAPSET)
    try:
        return mw(Seq(s), seq_type="DNA" if seq_type == "dna" else "RNA" if seq_type == "rna" else "protein")
    except Exception:
        return None
