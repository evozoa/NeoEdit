"""IQ-TREE wrapper: input writing, name mapping, command line; a real run when IQ-TREE is installed."""
import os
import sys

import pytest

from neoedit.analysis import phylo as P
from neoedit.model.alignment import SequenceRow

HERE = os.path.dirname(__file__)
EXAMPLE = os.path.join(HERE, "..", "examples", "cox1_demo.fasta")


def _rows():
    # two names that collide once IQ-TREE sanitizes them, plus Newick-special characters
    return [SequenceRow("taxon A (COI); x", "ATGCATGCATGCAAGG.TCCAAGGTTAACCGG"),
            SequenceRow("taxon A (COI): x", "ATGCATGCATGCAAGGTTCCAAGGTTAACCGA"),
            SequenceRow("Bob's fish", "ATGCATGCAAGCAAGGTTCCAAGG~~AACCGG"),
            SequenceRow("plain_name", "ATGGATGCAAGCAAGGTTCCTAGGTTAACCGG")]


def test_newick_label_quoting():
    assert P.newick_label("plain_name") == "plain_name"
    assert P.newick_label("a b") == "'a b'"
    assert P.newick_label("Bob's fish") == "'Bob''s fish'"
    assert P.newick_label("x:1") == "'x:1'"


def test_rename_newick_only_touches_leaves():
    names = ["A", "B c", "D"]
    tree = "(s0:0.1,(s1:0.2,s2:0.3)95/100:0.05,s10:0.1);"
    assert P.rename_newick(tree, names) == "(A:0.1,('B c':0.2,D:0.3)95/100:0.05,s10:0.1);"


def test_prepare_run_writes_safe_ids(tmp_path):
    run = P.prepare_run(_rows(), str(tmp_path / "out"), "demo", "dna")
    text = open(run.input_path).read()
    assert text.splitlines()[0::2] == [">s0", ">s1", ">s2", ">s3"]
    assert "." not in text and "~" not in text          # BioEdit gap characters become '-'
    names = open(run.prefix + ".names.tsv").read().splitlines()
    assert names[1] == "s0\ttaxon A (COI); x"
    sub = P.prepare_run(_rows(), str(tmp_path / "cols"), "demo", "dna", columns=(4, 10))
    assert sub.nsites == 6


def test_prepare_run_rejects_bad_input(tmp_path):
    with pytest.raises(ValueError, match="three"):
        P.prepare_run(_rows()[:2], str(tmp_path), "x", "dna")
    rows = _rows(); rows[0].seq += "A"
    with pytest.raises(ValueError, match="not aligned"):
        P.prepare_run(rows, str(tmp_path), "x", "dna")


def test_iqtree_args(tmp_path):
    run = P.prepare_run(_rows(), str(tmp_path), "demo", "protein")
    args = P.iqtree_args(run, model="LG+G4", bootstrap=1, replicates=200, alrt=1000,
                         threads=2, seed=7, outgroup=[2], extra_args=["--quiet"])
    assert args[:6] == ["-s", run.input_path, "--prefix", run.prefix, "-st", "AA"]
    for pair in (["-m", "LG+G4"], ["-T", "2"], ["-B", "1000"], ["--alrt", "1000"],
                 ["--seed", "7"], ["-o", "s2"]):
        i = args.index(pair[0]); assert args[i:i + 2] == pair
    assert args[-1] == "--quiet"
    plain = P.iqtree_args(run)
    assert "-B" not in plain and "-b" not in plain and plain[plain.index("-T") + 1] == "AUTO"
    cap = int(plain[plain.index("--threads-max") + 1])
    assert cap == max(1, (os.cpu_count() or 1) - 1)
    assert "--threads-max" not in args                 # an explicit thread count is not capped


def test_unique_dir(tmp_path):
    p = str(tmp_path / "a_iqtree")
    assert P.unique_dir(p) == p
    os.makedirs(p)
    assert P.unique_dir(p) == p + "_2"


def test_version_problem():
    assert P.version_problem("3.1.4") is None
    assert P.version_problem("2.4.0") is None
    assert "too old" in P.version_problem("1.6.12")
    assert P.version_problem("unknown")


needs_iqtree = pytest.mark.skipif(not P.find_iqtree(), reason="IQ-TREE not installed")


@needs_iqtree
def test_real_run_restores_names(tmp_path):
    run = P.prepare_run(_rows(), str(tmp_path / "run"), "demo", "dna")
    res = P.run_iqtree(run, model="JC", threads=1, seed=1, timeout=300)
    assert res.best_model == "JC"
    assert float(res.log_likelihood) < 0
    for r in _rows():
        assert P.newick_label(r.name) in res.newick
    assert open(run.tree_path).read().strip() == res.newick
    from Bio import Phylo
    tree = Phylo.read(run.tree_path, "newick")
    assert sorted(t.name for t in tree.get_terminals()) == sorted(r.name for r in _rows())


@needs_iqtree
def test_dialog_runs_iqtree(tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtCore import QSettings
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication
    from neoedit.model import io as mio
    from neoedit.ui.dialogs.tree_dialogs import IQTreeDialog
    QApplication.instance() or QApplication(sys.argv)
    settings = QSettings(str(tmp_path / "s.ini"), QSettings.IniFormat)
    model = mio.load(EXAMPLE)
    d = IQTreeDialog(None, settings, model, [0, 1, 2], (0, 100))
    d.out_dir.setText(str(tmp_path))
    d.model_box.setCurrentIndex(d.model_box.findData("HKY+F+G4"))
    d.boot.setCurrentIndex(0)
    d.alrt.setChecked(False)
    d.threads.setValue(1)
    d.cols_scope.setCurrentIndex(1)
    d.start()
    for _ in range(600):
        if not d.running():
            break
        QTest.qWait(100)
    assert d.tree is not None, d.summary.text() + d.log.toPlainText()[-500:]
    assert d.run.nsites == 100
    assert d.run.out_dir == str(tmp_path / "cox1_demo_iqtree")
    assert d.tree.best_model == "HKY+F+G4"
    assert d.b_copy.isEnabled() and not d.b_stop.isEnabled()
    assert settings.value("iqtree/model_dna") == "HKY+F+G4"
    d.close()
