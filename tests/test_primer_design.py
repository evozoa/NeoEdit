import random
import pytest
from neoedit.analysis import primer_design as PD
from neoedit.analysis.primers import PrimerPair, Primer


def test_matches_iupac():
    assert PD.matches("A", "A") and not PD.matches("A", "C")
    assert PD.matches("R", "A") and PD.matches("R", "G") and not PD.matches("R", "C")
    assert PD.matches("A", "N") and PD.matches("N", "T")
    assert not PD.matches("A", "-")


def test_conservation_and_mask():
    rows = ["ACGTACGTAC", "ACGTACGTAC", "ACGTTCGTAC"]
    cons = PD.conservation(rows)
    assert cons[0] == 1.0 and abs(cons[4] - 2 / 3) < 1e-9
    bad = PD.masked_regions(cons, 0.9, 3)
    # position 4 is variable -> excluded; runs 0-3 and 5-9 are conserved
    assert (4, 1) in bad


def test_degenerate_consensus():
    rows = ["ACGT", "ACAT", "ACGT"]
    assert PD.degenerate_consensus(rows, threshold=0.1) == "ACRT"
    assert PD.degenerate_consensus(rows, threshold=0.5) == "ACGT"


def test_score_primer_orientation_and_3prime():
    # row: 20 bases; forward primer = first 10 bases
    seq = "ACGTACGTAAGGCCTTAACC"
    names = ["s"]
    hits = PD.score_primer(seq[:10], [seq], names, 0, 10, 1)
    assert hits[0].mismatches == 0 and hits[0].covered
    # mismatch at the 3' end
    p = seq[:9] + ("T" if seq[9] != "T" else "G")
    hits = PD.score_primer(p, [seq], names, 0, 10, 1)
    assert hits[0].mismatches == 1 and hits[0].mm_3prime == 1 and hits[0].last_base_mm
    assert not hits[0].amplifies()
    # mismatch at the 5' end is tolerated
    p = ("T" if seq[0] != "T" else "G") + seq[1:10]
    hits = PD.score_primer(p, [seq], names, 0, 10, 1)
    assert hits[0].mismatches == 1 and hits[0].mm_3prime == 0 and not hits[0].last_base_mm
    assert hits[0].amplifies()
    # reverse primer: revcomp of the last 10 bases
    from Bio.Seq import Seq
    rp = str(Seq(seq[10:20]).reverse_complement())
    hits = PD.score_primer(rp, [seq], names, 10, 20, -1)
    assert hits[0].mismatches == 0
    # gapped row still lines up
    gapped = seq[:5] + "---" + seq[5:]
    hits = PD.score_primer(seq[:10], [gapped], names, 0, 13, 1)
    assert hits[0].mismatches == 0
    # uncovered row
    hits = PD.score_primer(seq[:10], ["-" * 20], names, 0, 10, 1)
    assert not hits[0].covered and not hits[0].amplifies()


def _make_alignment():
    random.seed(7)
    core = "".join(random.choice("ACGT") for _ in range(600))
    rows = []
    # include group: 4 sequences differing only in a variable middle window
    for k in range(4):
        s = list(core)
        for i in range(250, 350, 3):
            if random.random() < 0.5:
                s[i] = random.choice("ACGT")
        rows.append("".join(s))
    # exclude group: 2 sequences that differ at the 3' end of a conserved window too
    for k in range(2):
        s = list(core)
        for i in range(250, 350, 3):
            s[i] = random.choice("ACGT")
        for i in range(100, 140):
            if random.random() < 0.25:
                s[i] = random.choice("ACGT")
        rows.append("".join(s))
    names = [f"inc{i}" for i in range(4)] + [f"exc{i}" for i in range(2)]
    return rows, names


def test_design_reports_mask_fallback():
    rows, names = _make_alignment()
    # impossible constraint: 60-bp fully conserved runs do not exist here
    res = PD.design_on_alignment(rows, names, template_row=0, include_rows=[0, 1, 2, 3],
                                 min_conservation=1.0, min_conserved_run=300,   # longer than any conserved run here
                                 product_range=((150, 400),), num_return=5)
    assert not res.mask_applied and "conserved" in res.notes and len(res) > 0
    res2 = PD.design_on_alignment(rows, names, template_row=0, include_rows=[0, 1, 2, 3],
                                  min_conservation=1.0, min_conserved_run=20,
                                  product_range=((150, 500),), num_return=10)
    assert res2.mask_applied and res2.n_conserved_runs >= 1


def test_design_universal_and_table():
    rows, names = _make_alignment()
    evs = PD.design_on_alignment(rows, names, template_row=0, include_rows=[0, 1, 2, 3],
                                 min_conservation=1.0, min_conserved_run=20,
                                 product_range=((150, 500),), num_return=10)
    assert evs
    best = evs[0]
    st = best.stats()
    assert st["include_total"] == 4 and st["include_hit"] == 4       # universal across include set
    tbl = PD.insilico_table(best)
    assert len(tbl) == 6 and all(r["amplifies"] for r in tbl if r["group"] == "include")
    # primers must sit in conserved columns: zero mismatches in the include set
    assert all(h.mismatches == 0 for h in best.left_hits if h.row < 4)
    assert all(h.mismatches == 0 for h in best.right_hits if h.row < 4)


def test_design_discriminating_prefers_exclusion():
    rows, names = _make_alignment()
    inc, exc = [0, 1, 2, 3], [4, 5]
    evs = PD.design_on_alignment(rows, names, template_row=0, include_rows=inc, exclude_rows=exc,
                                 discriminating=True, min_conservation=1.0, min_conserved_run=18,
                                 product_range=((80, 400),), num_return=30, max_mm=2)
    assert evs
    scored = [(e.stats(max_mm=2)["include_frac"], e.stats(max_mm=2)["exclude_frac"]) for e in evs]
    # the top-ranked pair should not be worse at excluding than the median candidate
    top_inc, top_exc = scored[0]
    assert top_inc == 1.0
    assert top_exc <= sorted(x[1] for x in scored)[len(scored) // 2]


def test_degenerate_primer_for():
    rows = ["ACGTACGTAC", "ACATACGTAC", "ACGTACGTAC"]
    p, deg = PD.degenerate_primer_for(rows, [0, 1, 2], 0, 10, 1, threshold=0.1)
    assert p == "ACRTACGTAC" and deg == 2
    p2, deg2 = PD.degenerate_primer_for(rows, [0, 1, 2], 0, 10, 1, threshold=0.1, max_degeneracy=1)
    assert deg2 == 1 and p2 == "ACGTACGTAC"
