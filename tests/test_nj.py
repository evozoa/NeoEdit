"""Neighbor joining: distances checked against R ape 5.8.1 (dist.dna), NJ against known trees."""
import csv
import os
import sys

import numpy as np
import pytest

from neoedit.analysis import nj
from neoedit.model.alignment import SequenceRow

# 5 sequences x 80 sites with gaps, N and R (IUPAC); expected values from ape::dist.dna
SMALL = [
    "ACTCAAC--GAGCTAGAACCACCCCCAAC-A-ATCCCTGATCCGCTAC-AAATTCTATGNAAGAGCGAGGAACAAATGCA",
    "ACTCAACAAGAGCTGGAARAACCCCCAACCAA-TCCCTGATC-GCT-GTAATTTRTAT-CAAGAGCGAGGTACAAAGGCA",
    "ACTCCAAAAGAGCTTGATCCAGCCCCAACC-AATCCCTGATCCGCTACTACATTCCATCGAATCGCTAGGTACAAACGCA",
    "ACT-AACAAG-GCTTGAACCACCCCCAACCAAATCC-TGATCCGCTACTAAATTCTAT-GAAGCGCGAGGTACAAATGCA",
    "ACTCGACAAGAGCTTGAACCACCCACAACCAAA-CCCTGNTCCGCTAGTAAATTCT-TGGCAGAGCGAGGAACAAATGAA",
]
APE = {  # (model, pairwise deletion) -> upper triangle, row by row
    ("p", True): [0.088235294118, 0.178082191781, 0.042857142857, 0.084507042254, 0.205479452055,
                  0.098591549296, 0.140845070423, 0.120000000000, 0.210526315789, 0.095890410959],
    ("jc69", True): [0.093872357216, 0.203308438050, 0.044130375017, 0.089658862987, 0.240125645697,
                     0.105702255473, 0.156000428409, 0.130765040359, 0.247109400848, 0.102598726327],
    ("k2p", True): [0.094083094159, 0.204666667407, 0.044304351002, 0.089849900058, 0.244653324501,
                    0.106805549143, 0.157264648282, 0.130956077429, 0.250028247215, 0.102932479870],
    ("tn93", True): [0.094519320023, 0.207446268027, 0.044360333031, 0.090246545619, 0.246895285404,
                     0.107141815616, 0.158078756969, 0.132697057290, 0.253472818618, 0.103375386468],
    ("p", False): [0.098360655738, 0.196721311475, 0.049180327869, 0.098360655738, 0.229508196721,
                   0.098360655738, 0.147540983607, 0.147540983607, 0.245901639344, 0.114754098361],
    ("jc69", False): [0.105436462966, 0.228158530802, 0.050866947254, 0.105436462966, 0.273974299787,
                      0.105436462966, 0.164290174547, 0.164290174547, 0.297976348102, 0.124545776942],
    ("k2p", False): [0.105706966667, 0.229688580599, 0.051100596549, 0.105706966667, 0.279997684513,
                     0.106533730986, 0.165576959325, 0.164606994988, 0.302223565798, 0.125054366826],
    ("tn93", False): [0.106262261811, 0.233624993108, 0.051175338101, 0.106262261811, 0.283062281533,
                      0.106868239858, 0.166541721315, 0.167458800302, 0.307772399928, 0.125717244098],
}


def _rows(seqs=SMALL):
    return [SequenceRow(f"t{i}", s) for i, s in enumerate(seqs)]


def _read_matrix(path):
    rows = list(csv.reader(open(path, encoding="utf-8")))
    return rows[0][1:], np.array([r[1:] for r in rows[1:]], float)


@pytest.mark.parametrize("model,pairwise", sorted(APE))
def test_distances_match_ape(tmp_path, model, pairwise):
    res = nj.run_nj(_rows(), str(tmp_path), "s", "dna", model=model,
                    gaps="pairwise" if pairwise else "complete")
    names, D = _read_matrix(res.distance_path)
    assert names == [f"t{i}" for i in range(5)]
    assert np.allclose(D[np.triu_indices(5, 1)], APE[(model, pairwise)], atol=1e-6)
    assert res.nsites == (80 if pairwise else 61)


def test_textbook_example_is_reproduced():
    # Wikipedia's neighbor-joining example: an additive matrix, so the tree fits it exactly
    D = np.array([[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7], [9, 10, 8, 0, 3], [8, 9, 7, 3, 0]], float)
    t = nj.neighbor_joining(D)
    assert np.allclose(t.path_lengths(), D)
    assert t.newick(list("abcde")) == "(((a:2,b:3):3,c:4):2,d:2,e:1);"
    assert t.newick(list("abcde"), outgroup=2) == "(c:4,(a:2,b:3):3,(d:2,e:1):2);"


def test_random_additive_trees_are_recovered():
    rng = np.random.default_rng(0)
    for _ in range(10):
        n = 25
        t = nj.Tree(n)
        pool, nxt = list(range(n)), n
        while len(pool) > 3:
            a, b = sorted(rng.choice(len(pool), 2, replace=False))
            t.add_edge(nxt, pool[a], rng.uniform(0.01, 1)); t.add_edge(nxt, pool[b], rng.uniform(0.01, 1))
            pool = [p for k, p in enumerate(pool) if k not in (a, b)] + [nxt]
            nxt += 1
        for p in pool:
            t.add_edge(nxt, p, rng.uniform(0.01, 1))
        t.center = nxt
        D = t.path_lengths()
        got = nj.neighbor_joining(D)
        assert np.allclose(got.path_lengths(), D)
        assert got.splits() == t.splits()


def test_negative_branch_lengths_are_zeroed():
    D = np.array([[0, 1, 5, 5], [1, 0, 5, 5], [5, 5, 0, 20], [5, 5, 20, 0]], float)
    t = nj.neighbor_joining(D)
    assert all(w >= 0 for edges in t.adj.values() for _, w in edges)


def test_protein_poisson():
    rows = [SequenceRow("a", "MKVLAAGIVG"), SequenceRow("b", "MKVLSAGIVA"), SequenceRow("c", "MRVLSAGLVA")]
    X, c = nj.site_patterns(nj.encode([r.seq for r in rows], "protein"))
    D = nj.distances(nj.pair_counts(X, c[None], dna=False), "poisson")[0]
    assert D[0, 1] == pytest.approx(-np.log(1 - 0.2))
    assert D[0, 2] == pytest.approx(-np.log(1 - 0.4))
    with pytest.raises(ValueError, match="not available"):
        nj.run_nj(rows, "unused", "x", "protein", model="k2p")


def test_undefined_distances_are_reported(tmp_path):
    far = ["AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG", "ACGTACGTAC"]
    with pytest.raises(ValueError, match="too\\s+different"):
        nj.run_nj(_rows(far), str(tmp_path), "x", "dna", model="jc69")
    no_overlap = ["ACGT------", "----ACGTAC", "ACGTACGTAC", "ACGAACGTAC"]
    with pytest.raises(ValueError, match="share no comparable sites"):
        nj.run_nj(_rows(no_overlap), str(tmp_path), "x", "dna", model="p")
    with pytest.raises(ValueError, match="complete deletion"):
        nj.run_nj(_rows(no_overlap), str(tmp_path), "x", "dna", model="p", gaps="complete")
    with pytest.raises(ValueError, match="not aligned"):
        nj.run_nj(_rows(SMALL[:3] + ["ACGT"]), str(tmp_path), "x", "dna")


def test_bootstrap_support_and_outputs(tmp_path):
    # two well-separated pairs plus two in between: the (t0,t1) and (t4,t5) splits are solid
    base = "ACGTTGCAAGCTTAGCCGATACGGATCCTAGCATGCAATGCCGTAGCTAGCTAGGATCCGATCGATTGCA"
    def mutate(s, sites, rng):
        s = list(s)
        for i in sites:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
        return "".join(s)
    rng = np.random.default_rng(1)
    a = mutate(base, range(0, 60, 3), rng)
    b = mutate(base, range(1, 60, 3), rng)
    seqs = [mutate(a, [5], rng), mutate(a, [9], rng), base, mutate(base, [30, 40], rng),
            mutate(b, [7], rng), mutate(b, [11], rng)]
    rows = [SequenceRow(n, s) for n, s in zip(["A one", "A two", "mid", "mid 2", "B one", "B's two"], seqs)]
    r1 = nj.run_nj(rows, str(tmp_path / "r1"), "demo", "dna", model="k2p", replicates=200, seed=7, outgroup=2)
    r2 = nj.run_nj(rows, str(tmp_path / "r2"), "demo", "dna", model="k2p", replicates=200, seed=7, outgroup=2)
    assert r1.newick == r2.newick                   # same seed, same supports
    assert r1.newick.startswith("(mid:")            # outgroup first
    assert "(A_one:" in r1.newick or "(A_two:" in r1.newick
    assert "B_s_two:" in r1.newick and "'" not in r1.newick
    assert ")100:" in r1.newick
    from Bio import Phylo
    tree = Phylo.read(r1.tree_path, "newick")
    assert sorted(t.name for t in tree.get_terminals()) == ["A_one", "A_two", "B_one", "B_s_two", "mid", "mid_2"]
    assert open(str(tmp_path / "r1" / "demo.names.tsv")).read().splitlines()[-1] == "B_s_two\tB's two"
    report = open(r1.report_path, encoding="utf-8").read()
    assert "200 replicates, random seed 7" in report and "Outgroup:         mid" in report
    assert r1.seed == 7


def test_cancel_and_progress(tmp_path):
    seen = []
    with pytest.raises(nj.Cancelled):
        nj.run_nj(_rows(), str(tmp_path), "x", "dna", replicates=500, seed=1,
                  progress=lambda d, t: seen.append(d), cancelled=lambda: len(seen) >= 3)
    assert 3 <= len(seen) < 500


def test_nj_dialog(tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtCore import QSettings
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication
    from neoedit.model.alignment import AlignmentModel
    from neoedit.ui.dialogs.tree_dialogs import NJDialog
    QApplication.instance() or QApplication(sys.argv)
    settings = QSettings(str(tmp_path / "s.ini"), QSettings.IniFormat)
    model = AlignmentModel(_rows())
    d = NJDialog(None, settings, model, [], None)
    d.out_dir.setText(str(tmp_path))
    d.dist_box.setCurrentIndex(d.dist_box.findData("tn93"))
    d.reps.setValue(50)
    d.outgroup.setCurrentIndex(d.outgroup.findData(3))
    d.start()
    for _ in range(300):
        if not d.running():
            break
        QTest.qWait(50)
    QTest.qWait(50)
    assert d.tree is not None, d.summary.text()
    assert d.out_folder == str(tmp_path / "alignment_nj")
    assert d.tree.newick.startswith("(t3:")
    assert "Tamura–Nei" in d.log.toPlainText()
    assert d.b_copy.isEnabled() and d.b_run.isEnabled() and not d.b_stop.isEnabled()
    assert settings.value("nj/model_dna") == "tn93"
    d.close()
