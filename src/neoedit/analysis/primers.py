"""Primer design via primer3-py (bundled Primer3, no external binary)."""
from __future__ import annotations

from dataclasses import dataclass, field

import primer3

from ..model.alignment import GAP_CHARS, Feature

GAPSET = set(GAP_CHARS)


@dataclass
class Primer:
    seq: str
    start: int        # 0-based on ungapped template, 5' end of primer (for right primer: rightmost base index)
    length: int
    tm: float
    gc: float
    self_any: float
    self_end: float
    hairpin: float
    end_stability: float = 0.0

    @property
    def span(self):
        return self.start, self.start + self.length

    def pieces(self, length: int | None = None) -> list[tuple[int, int]]:
        """Ungapped [start,end) span(s) on the molecule: a right primer past the origin of a circular
        template is stored in unwrapped coordinates (start may be >= length) and comes back modulo
        length; a site that itself crosses the origin gives two pieces."""
        s, e = self.span
        if length:
            while s >= length:
                s -= length; e -= length
            if e > length:
                return [(s, length), (0, e - length)]
        return [(s, e)]

    def pos_text(self, length: int | None = None) -> str:
        pcs = self.pieces(length)
        if len(pcs) > 1:
            return f"{pcs[0][0] + 1}-{pcs[1][1]} (across origin)"
        s, e = pcs[0]
        return f"{s + 1}-{e}"


@dataclass
class PrimerPair:
    left: Primer
    right: Primer
    product_size: int
    penalty: float
    compl_any: float
    compl_end: float
    extra: dict = field(default_factory=dict)

    @property
    def circular_len(self) -> int | None:
        return self.extra.get("circular_len")

    def product_span(self) -> tuple[int, int]:
        """Ungapped [start,end) of the product; end may exceed the template length on a circle."""
        s = self.left.start
        return s, s + self.product_size

    @property
    def crosses_origin(self) -> bool:
        L = self.circular_len
        return bool(L) and self.product_span()[1] > L


def design_primers(template: str, target: tuple[int, int] | None = None,
                   product_range=((100, 300),), opt_size=20, min_size=18, max_size=25,
                   opt_tm=60.0, min_tm=57.0, max_tm=63.0, min_gc=30.0, max_gc=70.0,
                   num_return=5, excluded=None, max_poly_x=4, gc_clamp=0,
                   salt_monovalent=50.0, salt_divalent=1.5, dntp=0.6, dna_conc=50.0,
                   included: tuple[int, int] | None = None, circular: bool = False) -> list[PrimerPair]:
    """Design primer pairs on an (ungapped) template.

    target: (start, length) region that the product must span (0-based).
    included: (start, length) region to restrict primer picking to.
    circular: the template is a ring — Primer3 runs on the doubled sequence so primers and
    products may cross the origin; positions are reported modulo the template length
    (a primer/product past the origin has end > length; see Primer.pieces / PrimerPair.product_span).
    """
    tmpl = "".join(c for c in template if c not in GAPSET).upper().replace("U", "T")
    L = len(tmpl)
    maxprod = max((int(b) for _a, b in product_range), default=300)
    if target and not circular and int(target[0]) + int(target[1]) > L:
        raise ValueError(f"The target ({int(target[0]) + 1}-{int(target[0]) + int(target[1])}) runs past the end of the "
                         f"{L:,}-bp template. If the molecule is circular, mark it so (Topology) and the product may span the origin.")
    if circular and L:
        # One turn plus the longest product: every cross-origin primer/product is visible once.
        # With no target, restrict picking to [0, L) for the left primer so pairs are not found twice.
        ext = min(L, maxprod + 2 * max_size)
        tmpl2 = tmpl + tmpl[:ext]
        seq_args = {"SEQUENCE_ID": "template", "SEQUENCE_TEMPLATE": tmpl2}
        if target:
            seq_args["SEQUENCE_TARGET"] = [int(target[0]), int(target[1])]
        if included:
            seq_args["SEQUENCE_INCLUDED_REGION"] = [int(included[0]), int(included[1])]
        if excluded:
            # excluded regions that themselves cross the origin are expressed twice on the doubled template
            ex = []
            for a, b in excluded:
                ex.append([int(a), int(b)])
                if int(a) < ext:
                    ex.append([int(a) + L, min(int(b), ext)])
            seq_args["SEQUENCE_EXCLUDED_REGION"] = ex
    else:
        seq_args = {"SEQUENCE_ID": "template", "SEQUENCE_TEMPLATE": tmpl}
        if target:
            seq_args["SEQUENCE_TARGET"] = [int(target[0]), int(target[1])]
        if included:
            seq_args["SEQUENCE_INCLUDED_REGION"] = [int(included[0]), int(included[1])]
        if excluded:
            seq_args["SEQUENCE_EXCLUDED_REGION"] = [[int(a), int(b)] for a, b in excluded]
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1, "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_OPT_SIZE": opt_size, "PRIMER_MIN_SIZE": min_size, "PRIMER_MAX_SIZE": max_size,
        "PRIMER_OPT_TM": opt_tm, "PRIMER_MIN_TM": min_tm, "PRIMER_MAX_TM": max_tm,
        "PRIMER_MIN_GC": min_gc, "PRIMER_MAX_GC": max_gc,
        "PRIMER_MAX_POLY_X": max_poly_x, "PRIMER_GC_CLAMP": gc_clamp,
        "PRIMER_SALT_MONOVALENT": salt_monovalent, "PRIMER_SALT_DIVALENT": salt_divalent,
        "PRIMER_DNTP_CONC": dntp, "PRIMER_DNA_CONC": dna_conc,
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(r) for r in product_range],
        "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_MAX_SELF_ANY_TH": 45.0, "PRIMER_MAX_SELF_END_TH": 35.0, "PRIMER_MAX_HAIRPIN_TH": 24.0,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 35.0,
    }
    if circular and L and not target and not included:
        # ask for more than needed: candidates with a left primer on the second turn are duplicates
        global_args["PRIMER_NUM_RETURN"] = num_return * 3
    res = primer3.bindings.design_primers(seq_args, global_args)
    pairs = []
    seen = set()
    n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
    for i in range(n):
        ls, ll = res[f"PRIMER_LEFT_{i}"]
        rs, rl = res[f"PRIMER_RIGHT_{i}"]
        if circular and L:
            if ls >= L:
                continue                                   # same pair, one turn later
            key = (ls % L, (rs - rl + 1) % L)
            if key in seen:
                continue
            seen.add(key)
            rs = (rs - rl + 1) % L + rl - 1                # keep right.start (5' end on the top strand) in [0, L) + len
            if rs - rl + 1 < ls:                           # right primer wrapped back before left: unwrap it
                rs += L
        if len(pairs) >= num_return:
            break
        left = Primer(res[f"PRIMER_LEFT_{i}_SEQUENCE"], ls, ll, res[f"PRIMER_LEFT_{i}_TM"],
                      res[f"PRIMER_LEFT_{i}_GC_PERCENT"], res.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH", 0),
                      res.get(f"PRIMER_LEFT_{i}_SELF_END_TH", 0), res.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH", 0),
                      res.get(f"PRIMER_LEFT_{i}_END_STABILITY", 0))
        # primer3 reports right primer start as the 5' end on the top-strand coordinate = rightmost base
        right = Primer(res[f"PRIMER_RIGHT_{i}_SEQUENCE"], rs - rl + 1, rl, res[f"PRIMER_RIGHT_{i}_TM"],
                       res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], res.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH", 0),
                       res.get(f"PRIMER_RIGHT_{i}_SELF_END_TH", 0), res.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH", 0),
                       res.get(f"PRIMER_RIGHT_{i}_END_STABILITY", 0))
        extra = {"explain": {k: v for k, v in res.items() if k.endswith("_EXPLAIN")}}
        if circular and L:
            extra["circular_len"] = L
        pairs.append(PrimerPair(left, right, res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                                res[f"PRIMER_PAIR_{i}_PENALTY"],
                                res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH", 0),
                                res.get(f"PRIMER_PAIR_{i}_COMPL_END_TH", 0),
                                extra=extra))
    if not pairs:
        expl = {k: v for k, v in res.items() if k.endswith("_EXPLAIN")}
        raise ValueError("No primer pairs found. " + "; ".join(f"{k}: {v}" for k, v in expl.items()))
    return pairs


def primer_stats(seq: str, **kw) -> dict:
    s = seq.upper().replace("U", "T")
    return {
        "tm": primer3.bindings.calc_tm(s, **kw),
        "gc": 100.0 * sum(c in "GC" for c in s) / max(1, len(s)),
        "hairpin": primer3.bindings.calc_hairpin(s).tm,
        "homodimer": primer3.bindings.calc_homodimer(s).tm,
    }


def heterodimer_tm(a: str, b: str) -> float:
    return primer3.bindings.calc_heterodimer(a.upper(), b.upper()).tm


def primer_mismatches(primer: str, templates: list[str], strand: int, ungapped_start: int) -> list[int]:
    """Count mismatches of a primer vs each (ungapped) template at its position.
    For strand -1, the primer is given 5'->3' and compared to rev-comp of the site."""
    from Bio.Seq import Seq
    out = []
    p = primer.upper()
    for t in templates:
        t = "".join(c for c in t if c not in GAPSET).upper()
        site = t[ungapped_start:ungapped_start + len(p)]
        if strand < 0:
            site = str(Seq(site).reverse_complement())
        mm = sum(1 for x, y in zip(p, site) if x != y and y != "N") + abs(len(p) - len(site))
        out.append(mm)
    return out


def pair_to_features(pair: PrimerPair, row: int, gmap: list[int] | None, idx: int) -> list[Feature]:
    feats = []
    L = pair.circular_len
    for prm, strand, name in ((pair.left, 1, f"F{idx}"), (pair.right, -1, f"R{idx}")):
        for s, e in prm.pieces(L):                      # two pieces when the site crosses the origin
            if gmap:
                s2 = gmap[s] if s < len(gmap) else len(gmap)
                e2 = gmap[e - 1] + 1 if e - 1 < len(gmap) else len(gmap)
            else:
                s2, e2 = s, e
            feats.append(Feature(row, s2, e2, strand, "primer", f"{name} Tm {prm.tm:.1f}",
                                 "#10b981" if strand > 0 else "#ef4444",
                                 data={"seq": prm.seq, "tm": prm.tm, "gc": prm.gc, "pair": idx}))
    return feats
