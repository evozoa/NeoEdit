"""Residue colour schemes.

Each scheme maps an uppercase residue char to a hex colour; None means "no
colour". The default schemes are BioEdit's own colour tables (parsed from the
`color.tab` files shipped with BioEdit, see resources/bioedit_tables/).
"""
from __future__ import annotations

import os
import re

_TABLE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resources", "bioedit_tables")


def _bgr_to_hex(v: int) -> str:
    """BioEdit stores colours as a 3-byte integer in BGR order (Windows COLORREF)."""
    v = int(v)
    if v < 0:
        v &= 0xFFFFFF
    b, g, r = (v >> 16) & 0xFF, (v >> 8) & 0xFF, v & 0xFF
    return f"#{r:02x}{g:02x}{b:02x}"


def parse_bioedit_table(text: str) -> dict[str, dict[str, str]]:
    """Parse a BioEdit .tab colour file -> {"protein": {...}, "dna": {...}}."""
    out = {"protein": {}, "dna": {}}
    cur = None
    lines = [ln.rstrip("\r\n") for ln in text.splitlines()]
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln == "/amino acids/":
            cur = "protein"; i += 1; continue
        if ln == "/nucleotides/":
            cur = "dna"; i += 1; continue
        if ln == "/////" or not ln:
            if ln == "/////":
                cur = None
            i += 1; continue
        if cur and i + 1 < len(lines):
            try:
                val = int(ln)
            except ValueError:
                try:
                    val = int(ln, 16)
                except ValueError:
                    i += 1; continue
            chars = lines[i + 1].strip()
            col = _bgr_to_hex(val)
            for ch in chars:
                if ch.isalpha():
                    out[cur][ch.upper()] = col
            i += 2
            continue
        i += 1
    d = out["dna"]
    if "T" in d and "U" not in d:
        d["U"] = d["T"]
    if "U" in d and "T" not in d:
        d["T"] = d["U"]
    return out


def _load_table(fname: str):
    try:
        with open(os.path.join(_TABLE_DIR, fname), "r", errors="replace") as fh:
            return parse_bioedit_table(fh.read())
    except OSError:
        return {"protein": {}, "dna": {}}


_BIOEDIT = _load_table("color.tab")          # the active default table shipped with BioEdit
_BIOEDIT_LEGACY = _load_table("defcolor.tab")
_BIOEDIT_BLOSUM = _load_table("BLOSUMcoloring.tab")
_BIOEDIT_PAM250 = _load_table("PAM250Coloring.tab")
_BIOEDIT_KD = _load_table("KyteDoolittleHydrophobicityColoring.tab")
_BIOEDIT_RUIZ = _load_table("ManuelRuizColorTable.tab")
_BIOEDIT_DEGEN = _load_table("codonDegeneracyColoring.tab")

# Fallback in case the resource files are missing: BioEdit's classic A/C/G/T
BIOEDIT_NT = _BIOEDIT["dna"] or {"A": "#008000", "C": "#0000ff", "G": "#000000", "T": "#ff0000", "U": "#ff0000", "N": "#808080"}
BIOEDIT_AA = _BIOEDIT["protein"]

# Brighter, background-friendly nucleotide palette
NUCLEOTIDE = {
    "A": "#5bd65b", "C": "#5b9bff", "G": "#f4c542", "T": "#ff6b6b", "U": "#ff6b6b",
    "N": "#c0c0c0", "R": "#bdbd7a", "Y": "#bd7abd", "K": "#7abdbd", "M": "#bd9a7a",
    "S": "#7a9abd", "W": "#9abd7a", "B": "#b0b0b0", "D": "#b0b0b0", "H": "#b0b0b0", "V": "#b0b0b0",
}

CLUSTALX = {
    "A": "#80a0f0", "I": "#80a0f0", "L": "#80a0f0", "M": "#80a0f0", "F": "#80a0f0", "W": "#80a0f0", "V": "#80a0f0", "C": "#80a0f0",
    "K": "#f01505", "R": "#f01505",
    "E": "#c048c0", "D": "#c048c0",
    "N": "#15c015", "Q": "#15c015", "S": "#15c015", "T": "#15c015",
    "G": "#f09048", "P": "#c0c000", "H": "#15a4a4", "Y": "#15a4a4",
}

ZAPPO = {
    "I": "#ffafaf", "L": "#ffafaf", "V": "#ffafaf", "A": "#ffafaf", "M": "#ffafaf",
    "F": "#ffc800", "W": "#ffc800", "Y": "#ffc800",
    "K": "#6464ff", "R": "#6464ff", "H": "#6464ff",
    "D": "#ff0000", "E": "#ff0000",
    "S": "#00ff00", "T": "#00ff00", "N": "#00ff00", "Q": "#00ff00",
    "P": "#ff00ff", "G": "#ff00ff", "C": "#ffff00",
}

TAYLOR = {
    "A": "#ccff00", "V": "#99ff00", "I": "#66ff00", "L": "#33ff00", "M": "#00ff00", "F": "#00ff66",
    "Y": "#00ffcc", "W": "#00ccff", "H": "#0066ff", "R": "#0000ff", "K": "#6600ff", "N": "#cc00ff",
    "Q": "#ff00cc", "E": "#ff0066", "D": "#ff0000", "S": "#ff3300", "T": "#ff6600", "G": "#ff9900",
    "P": "#ffcc00", "C": "#ffff00",
}

_KD = {"I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7,
       "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5,
       "N": -3.5, "K": -3.9, "R": -4.5}


def _lerp(a, b, t):
    return int(a + (b - a) * t)


def _hydro_color(v):
    t = (v + 4.5) / 9.0
    r, g, b = _lerp(0x40, 0xff, t), _lerp(0x80, 0x40, t), _lerp(0xff, 0x40, t)
    return f"#{r:02x}{g:02x}{b:02x}"


HYDROPHOBICITY = {k: _hydro_color(v) for k, v in _KD.items()}

# Ordered: first entry is the default
SCHEMES_NT = {
    "BioEdit": BIOEDIT_NT,
    "BioEdit (legacy default)": _BIOEDIT_LEGACY["dna"],
    "Bright": NUCLEOTIDE,
    "None": {},
}
SCHEMES_AA = {
    "BioEdit": BIOEDIT_AA,
    "BioEdit (legacy default)": _BIOEDIT_LEGACY["protein"],
    "BioEdit BLOSUM": _BIOEDIT_BLOSUM["protein"],
    "BioEdit PAM250": _BIOEDIT_PAM250["protein"],
    "BioEdit Kyte-Doolittle": _BIOEDIT_KD["protein"],
    "BioEdit Manuel Ruiz": _BIOEDIT_RUIZ["protein"],
    "BioEdit codon degeneracy": _BIOEDIT_DEGEN["protein"],
    "Clustal X": CLUSTALX,
    "Zappo": ZAPPO,
    "Taylor": TAYLOR,
    "Hydrophobicity": HYDROPHOBICITY,
    "None": {},
}
SCHEMES_NT = {k: v for k, v in SCHEMES_NT.items() if v or k == "None"}
SCHEMES_AA = {k: v for k, v in SCHEMES_AA.items() if v or k == "None"}


def schemes_for(seq_type: str) -> dict:
    return SCHEMES_NT if seq_type in ("dna", "rna") else SCHEMES_AA


def text_for(bg: str | None) -> str:
    if not bg:
        return "#000000"
    r, g, b = int(bg[1:3], 16), int(bg[3:5], 16), int(bg[5:7], 16)
    lum = 0.299 * r + 0.587 * g + 0.114 * b
    return "#000000" if lum > 140 else "#ffffff"
