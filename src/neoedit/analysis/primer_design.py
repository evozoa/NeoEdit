"""Alignment-aware primer design: conserved (universal) and discriminating (eDNA) primers.

The alignment is the source of truth. Candidate primers come from Primer3 designed
on a template (usually the reference row or the include-group consensus), restricted
to conserved windows, then every candidate is scored against every sequence in the
alignment with 3'-weighted mismatch penalties, giving an in-silico PCR table.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable, Sequence

from Bio.Seq import Seq

from ..model.alignment import GAP_CHARS
from . import primers as P

GAPSET = set(GAP_CHARS)

# IUPAC degenerate codes
IUPAC = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CG"): "S", frozenset("AT"): "W", frozenset("CGT"): "B", frozenset("AGT"): "D",
    frozenset("ACT"): "H", frozenset("ACG"): "V", frozenset("ACGT"): "N",
}
IUPAC_SETS = {v: set(k) for k, v in IUPAC.items()}
DEGENERACY = {c: len(s) for c, s in IUPAC_SETS.items()}


def matches(primer_base: str, template_base: str) -> bool:
    """True if a primer base can pair with a template base (IUPAC-aware)."""
    pb, tb = primer_base.upper(), template_base.upper()
    if tb in GAPSET or tb == "" :
        return False
    if tb == "N":
        return True
    return bool(IUPAC_SETS.get(pb, {pb}) & IUPAC_SETS.get(tb, {tb}))


# --------------------------------------------------------------------- columns
def column_stats(rows: Sequence[str], cols: range | Iterable[int] | None = None,
                 ignore_gaps: bool = True) -> list[tuple[str, float, int]]:
    """Per column: (plurality base, fraction agreeing, n informative)."""
    if not rows:
        return []
    w = max(len(r) for r in rows)
    cols = range(w) if cols is None else list(cols)
    out = []
    for c in cols:
        col = [r[c].upper() for r in rows if c < len(r)]
        obs = [x for x in col if x not in GAPSET and x != "N"] if ignore_gaps else col
        if not obs:
            out.append(("-", 0.0, 0)); continue
        base, n = Counter(obs).most_common(1)[0]
        out.append((base, n / len(obs), len(obs)))
    return out


def conservation(rows: Sequence[str], min_coverage: float = 0.5) -> list[float]:
    """Per-column conservation in [0,1]; columns covered by < min_coverage of the
    sequences are reported as 0 (unusable for primer placement)."""
    stats = column_stats(rows)
    n = len(rows)
    return [frac if cov >= max(1, min_coverage * n) else 0.0 for _base, frac, cov in stats]


def degenerate_consensus(rows: Sequence[str], threshold: float = 0.05, cols=None) -> str:
    """IUPAC consensus: any base present in > threshold of covered sequences is included."""
    if not rows:
        return ""
    w = max(len(r) for r in rows)
    cols = range(w) if cols is None else list(cols)
    out = []
    for c in cols:
        obs = [r[c].upper() for r in rows if c < len(r) and r[c] not in GAPSET]
        obs = [b for x in obs for b in IUPAC_SETS.get(x, set())]
        if not obs:
            out.append("N"); continue
        cnt = Counter(obs); tot = len(obs)
        keep = {b for b, k in cnt.items() if k / tot > threshold} or {cnt.most_common(1)[0][0]}
        out.append(IUPAC.get(frozenset(keep), "N"))
    return "".join(out)


def masked_regions(cons: Sequence[float], min_cons: float, min_run: int) -> list[tuple[int, int]]:
    """Regions (start,len) that primers must AVOID: everything not in a conserved run."""
    good = [i for i, v in enumerate(cons) if v >= min_cons]
    runs = []
    if good:
        s = prev = good[0]
        for i in good[1:]:
            if i != prev + 1:
                runs.append((s, prev + 1)); s = i
            prev = i
        runs.append((s, prev + 1))
    runs = [(a, b) for a, b in runs if b - a >= min_run]
    # complement of the runs
    bad = []
    pos = 0
    for a, b in runs:
        if a > pos:
            bad.append((pos, a - pos))
        pos = b
    if pos < len(cons):
        bad.append((pos, len(cons) - pos))
    return bad


# --------------------------------------------------------------------- scoring
@dataclass
class PrimerHit:
    """How one primer performs against one aligned sequence."""
    row: int
    name: str
    covered: bool                 # sequence actually spans the primer site
    mismatches: int
    mm_3prime: int                # mismatches in the last `three_prime_window` bases
    last_base_mm: bool
    site: str                     # template site, primer orientation
    tm: float | None = None

    def amplifies(self, max_mm: int = 3, max_3p: int = 0, require_last: bool = True) -> bool:
        if not self.covered:
            return False
        if require_last and self.last_base_mm:
            return False
        return self.mismatches <= max_mm and self.mm_3prime <= max_3p


@dataclass
class PairEvaluation:
    pair: P.PrimerPair
    left_hits: list[PrimerHit] = field(default_factory=list)
    right_hits: list[PrimerHit] = field(default_factory=list)
    include_rows: list[int] = field(default_factory=list)
    exclude_rows: list[int] = field(default_factory=list)
    left_cols: tuple[int, int] = (0, 0)     # alignment columns
    right_cols: tuple[int, int] = (0, 0)
    left_seq_deg: str = ""                  # degenerate versions (optional)
    right_seq_deg: str = ""

    def amplified_rows(self, **kw) -> list[int]:
        lh = {h.row: h for h in self.left_hits}
        rh = {h.row: h for h in self.right_hits}
        return [r for r in lh if r in rh and lh[r].amplifies(**kw) and rh[r].amplifies(**kw)]

    def stats(self, **kw) -> dict:
        amp = set(self.amplified_rows(**kw))
        inc = set(self.include_rows); exc = set(self.exclude_rows)
        return {
            "n_amplified": len(amp),
            "include_hit": len(amp & inc), "include_total": len(inc),
            "exclude_hit": len(amp & exc), "exclude_total": len(exc),
            "include_frac": (len(amp & inc) / len(inc)) if inc else 0.0,
            "exclude_frac": (len(amp & exc) / len(exc)) if exc else 0.0,
        }

    def score(self, discriminating: bool, **kw) -> float:
        """Higher is better. Universal: cover the include set. Discriminating: cover
        include, avoid exclude (heavily penalised), prefer 3'-end mismatches in exclude."""
        st = self.stats(**kw)
        s = 100.0 * st["include_frac"]
        if discriminating and self.exclude_rows:
            s -= 150.0 * st["exclude_frac"]
            # bonus for 3'-terminal mismatches against excluded rows (blocks extension)
            lh = {h.row: h for h in self.left_hits}; rh = {h.row: h for h in self.right_hits}
            bonus = 0.0
            for r in self.exclude_rows:
                for h in (lh.get(r), rh.get(r)):
                    if h and h.covered:
                        bonus += 6.0 if h.last_base_mm else (3.0 if h.mm_3prime else 0.0)
            s += min(30.0, bonus)
        s -= 2.0 * self.pair.penalty
        # prefer fewer total mismatches across the include set
        inc = set(self.include_rows)
        mm = sum(h.mismatches for h in self.left_hits + self.right_hits if h.row in inc)
        s -= 0.5 * mm
        return s


def _site_from_alignment(row_seq: str, c0: int, c1: int, strand: int) -> str:
    """Template site in primer orientation (gaps removed)."""
    site = "".join(ch for ch in row_seq[c0:c1] if ch not in GAPSET).upper()
    if strand < 0:
        site = str(Seq(site).reverse_complement())
    return site


def score_primer(primer: str, rows: Sequence[str], names: Sequence[str], c0: int, c1: int,
                 strand: int, three_prime_window: int = 5) -> list[PrimerHit]:
    """Compare `primer` (5'->3') to the site at alignment columns [c0,c1) in every row.

    For a forward primer (strand +1) the site is the top strand of the span; for a
    reverse primer (strand -1) it is the reverse complement, so in both cases primer
    index k lines up with site index k (primer 5' end first, 3' end last). Rows whose
    site is missing/short (gaps, sequence not covering the region) are marked
    `covered=False` and never count as amplifying.
    """
    hits: list[PrimerHit] = []
    L = len(primer)
    p = primer.upper()
    for i, (seq, name) in enumerate(zip(rows, names)):
        site = _site_from_alignment(seq, c0, c1, strand)
        covered = len(site) >= L - 2 and len(site) > 0
        mm = mm3 = 0
        last = False
        for k in range(L):
            tb = site[k] if k < len(site) else ""
            if not (tb and matches(p[k], tb)):
                mm += 1
                if k >= L - three_prime_window:
                    mm3 += 1
                if k == L - 1:
                    last = True
        hits.append(PrimerHit(i, name, covered, mm, mm3, last, site))
    return hits


def evaluate_pair(pair: P.PrimerPair, rows: Sequence[str], names: Sequence[str],
                  left_cols: tuple[int, int], right_cols: tuple[int, int],
                  include_rows: Sequence[int] = (), exclude_rows: Sequence[int] = (),
                  three_prime_window: int = 5) -> PairEvaluation:
    ev = PairEvaluation(pair, include_rows=list(include_rows), exclude_rows=list(exclude_rows),
                        left_cols=left_cols, right_cols=right_cols)
    ev.left_hits = score_primer(pair.left.seq, rows, names, *left_cols, 1, three_prime_window)
    ev.right_hits = score_primer(pair.right.seq, rows, names, *right_cols, -1, three_prime_window)
    return ev


def degenerate_primer_for(rows: Sequence[str], sel_rows: Sequence[int], c0: int, c1: int,
                          strand: int, threshold: float = 0.05, max_degeneracy: int = 64) -> tuple[str, int]:
    """IUPAC-degenerate primer covering the selected rows at this site. Returns
    (primer, degeneracy); falls back to the plurality base where degeneracy would explode."""
    sites = [_site_from_alignment(rows[i], c0, c1, strand) for i in sel_rows]
    L = max((len(s) for s in sites), default=0)
    sites = [s for s in sites if len(s) == L]
    if not sites:
        return "", 0
    cols = [[s[k] for s in sites] for k in range(L)]
    out = []
    for col in cols:
        obs = [b for x in col for b in IUPAC_SETS.get(x.upper(), set())]
        cnt = Counter(obs); tot = len(obs) or 1
        keep = {b for b, k in cnt.items() if k / tot > threshold} or {cnt.most_common(1)[0][0]}
        out.append(IUPAC.get(frozenset(keep), "N"))
    deg = 1
    for ch in out:
        deg *= DEGENERACY.get(ch, 4)
    while deg > max_degeneracy:
        # collapse the most degenerate column back to its plurality base
        worst = max(range(L), key=lambda k: DEGENERACY.get(out[k], 4))
        if DEGENERACY.get(out[worst], 1) <= 1:
            break
        col = [x.upper() for x in cols[worst]]
        out[worst] = Counter(col).most_common(1)[0][0]
        deg = 1
        for ch in out:
            deg *= DEGENERACY.get(ch, 4)
    return "".join(out), deg


@dataclass
class DesignResult:
    evaluations: list["PairEvaluation"]
    mask_applied: bool
    n_conserved_runs: int
    longest_run: int
    notes: str = ""

    def __iter__(self):
        return iter(self.evaluations)

    def __len__(self):
        return len(self.evaluations)

    def __getitem__(self, i):
        return self.evaluations[i]


def conserved_runs(cons: Sequence[float], min_cons: float, min_run: int) -> list[tuple[int, int]]:
    runs = []
    s = None
    for i, v in enumerate(cons):
        if v >= min_cons:
            s = i if s is None else s
        else:
            if s is not None:
                runs.append((s, i)); s = None
    if s is not None:
        runs.append((s, len(cons)))
    return [r for r in runs if r[1] - r[0] >= min_run]


def design_on_alignment(rows: Sequence[str], names: Sequence[str],
                        template_row: int = 0,
                        include_rows: Sequence[int] = (), exclude_rows: Sequence[int] = (),
                        region: tuple[int, int] | None = None,
                        discriminating: bool = False,
                        min_conservation: float = 0.9, min_conserved_run: int = 18,
                        product_range=((100, 500),), num_return: int = 20,
                        three_prime_window: int = 5,
                        max_mm: int = 3, max_3p: int = 0, require_last: bool = True,
                        primer3_kwargs: dict | None = None,
                        degenerate: bool = False, degeneracy_threshold: float = 0.05,
                        max_degeneracy: int = 64) -> list[PairEvaluation]:
    """Design primers on the alignment and rank them.

    rows/names are the *aligned* sequences. Primer3 designs on the template row
    (gaps removed), restricted to windows that are conserved across `include_rows`;
    every candidate is then scored against all rows in alignment coordinates.
    """
    inc = list(include_rows) or [i for i in range(len(rows)) if i != template_row] or [template_row]
    exc = list(exclude_rows)
    tmpl_aln = rows[template_row]
    # map: template ungapped index -> alignment column
    t2c = [i for i, ch in enumerate(tmpl_aln) if ch not in GAPSET]
    tmpl = "".join(tmpl_aln[i] for i in t2c).upper()
    if not tmpl:
        raise ValueError("Template row has no sequence.")
    c2t = {c: i for i, c in enumerate(t2c)}

    cons_aln = conservation([rows[i] for i in inc])
    # conservation projected onto template positions
    cons_t = [cons_aln[c] if c < len(cons_aln) else 0.0 for c in t2c]
    excluded = masked_regions(cons_t, min_conservation, min_conserved_run)

    # region of interest (alignment columns) -> template included region
    included = None
    if region:
        a, b = region
        ta = min((i for i, c in enumerate(t2c) if c >= a), default=0)
        tb = max((i for i, c in enumerate(t2c) if c < b), default=len(tmpl) - 1) + 1
        if tb > ta:
            included = (ta, tb - ta)

    runs = conserved_runs(cons_t, min_conservation, min_conserved_run)
    longest = max((b - a for a, b in runs), default=0)
    kw = dict(product_range=product_range, num_return=num_return)
    kw.update(primer3_kwargs or {})
    mask_applied = True
    notes = ""
    try:
        pairs = P.design_primers(tmpl, excluded=excluded, included=included, **kw)
    except ValueError as e:
        # nothing satisfied the conservation mask: retry unmasked, but say so loudly
        mask_applied = False
        notes = (f"No primer pair fits entirely inside conserved windows "
                 f"(conservation \u2265 {min_conservation:g}, run \u2265 {min_conserved_run} bp: "
                 f"{len(runs)} window(s), longest {longest} bp). "
                 f"Showing unrestricted candidates ranked by how well they match the include set \u2014 "
                 f"lower the conservation threshold or the minimum run, or widen the product range, "
                 f"for primers guaranteed to sit in conserved sites.")
        try:
            pairs = P.design_primers(tmpl, included=included, **kw)
        except ValueError:
            raise ValueError(str(e))

    evals = []
    for pr in pairs:
        lc = (t2c[pr.left.start], t2c[min(len(t2c) - 1, pr.left.start + pr.left.length - 1)] + 1)
        rs = pr.right.start
        rc = (t2c[rs], t2c[min(len(t2c) - 1, rs + pr.right.length - 1)] + 1)
        ev = evaluate_pair(pr, rows, names, lc, rc, inc, exc, three_prime_window)
        if degenerate:
            ev.left_seq_deg = degenerate_primer_for(rows, inc, *lc, 1, degeneracy_threshold, max_degeneracy)[0]
            ev.right_seq_deg = degenerate_primer_for(rows, inc, *rc, -1, degeneracy_threshold, max_degeneracy)[0]
        evals.append(ev)
    kwargs = dict(max_mm=max_mm, max_3p=max_3p, require_last=require_last)
    evals.sort(key=lambda e: -e.score(discriminating, **kwargs))
    return DesignResult(evals, mask_applied, len(runs), longest, notes)


def insilico_table(ev: PairEvaluation, max_mm: int = 3, max_3p: int = 0, require_last: bool = True) -> list[dict]:
    """One row per sequence: mismatches, 3' mismatches, predicted amplification."""
    lh = {h.row: h for h in ev.left_hits}
    rh = {h.row: h for h in ev.right_hits}
    inc, exc = set(ev.include_rows), set(ev.exclude_rows)
    out = []
    for r in sorted(set(lh) | set(rh)):
        L, R = lh.get(r), rh.get(r)
        kw = dict(max_mm=max_mm, max_3p=max_3p, require_last=require_last)
        amp = bool(L and R and L.amplifies(**kw) and R.amplifies(**kw))
        out.append({
            "row": r, "name": (L or R).name,
            "group": "include" if r in inc else ("exclude" if r in exc else ""),
            "fwd_mm": L.mismatches if L else None, "fwd_3p": L.mm_3prime if L else None,
            "fwd_last": L.last_base_mm if L else None, "fwd_cov": L.covered if L else False,
            "rev_mm": R.mismatches if R else None, "rev_3p": R.mm_3prime if R else None,
            "rev_last": R.last_base_mm if R else None, "rev_cov": R.covered if R else False,
            "amplifies": amp,
        })
    return out
