"""Offscreen GUI smoke test: drives the main window through key operations."""
import os
import sys
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt, QPoint
from PySide6.QtTest import QTest

HERE = os.path.dirname(__file__)
EXAMPLE = os.path.join(HERE, "..", "examples", "cox1_demo.fasta")


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def test_main_window_flow(app, tmp_path):
    from neoedit.ui.main_window import MainWindow
    from neoedit.analysis.orf_finder import find_orfs
    w = MainWindow()
    w.show()
    w.open_path(EXAMPLE)
    assert w.model.nrows == 5 and w.model.seq_type == "dna"
    v = w.view
    v.resize(900, 500)
    v.viewport().repaint()
    # keyboard: move, insert gap, undo
    v.set_cursor(0, 5)
    QTest.keyClick(v, Qt.Key_Space)
    assert w.model.rows[0].seq[5] == "-"
    w.undo()
    assert w.model.rows[0].seq[5] != "-"
    # insert mode typing
    v.set_mode("insert")
    v.set_cursor(0, 5)
    QTest.keyClick(v, Qt.Key_N)
    assert w.model.rows[0].seq[5] == "N"
    w.undo()
    v.set_mode("select")
    # selection + column gap
    v.select_region(0, 4, 10, 12)
    QTest.keyClick(v, Qt.Key_Space)
    assert all(r.seq[10] == "-" for r in w.model.rows)
    QTest.keyClick(v, Qt.Key_Backspace)  # removes the gap column again
    assert all(r.seq[10] != "-" for r in w.model.rows)
    # toggle translation overlay & colour modes; repaint
    w.a_translation.trigger(); v.viewport().repaint()
    w._set_color_mode("identity"); v.viewport().repaint()
    w.a_dots.trigger(); v.viewport().repaint()
    w._set_color_mode("scheme")
    # ORF finder: table 2
    orfs = find_orfs(w.model.rows[0].seq, table=2, min_aa=100, both_strands=True)
    assert orfs and orfs[0].length_aa >= 100
    from neoedit.ui.dialogs.orf_dialog import ORFFinderDialog
    d = ORFFinderDialog(w.model, [0, 1], w, 2)
    d.featuresReady.connect(w._add_features)
    d.run()
    assert d.table_w.rowCount() > 0
    d.table_w.selectRow(0)
    d.add_features()
    assert w.model.features
    v.viewport().repaint()
    # primer dialog
    from neoedit.ui.dialogs.primer_dialog import PrimerDialog
    p = PrimerDialog(w.model, [0], (200, 260), w)
    p.featuresReady.connect(w._add_features)
    p.run()
    assert p.pairs, p.status.text()
    p.table_w.selectRow(0)
    p.add_features()
    p.check_alignment()
    assert "mm" in p.detail.toPlainText()
    # revcomp + save roundtrip + export features
    w.model.reverse_complement([1]); w.undo()
    out = tmp_path / "out.aln"
    from neoedit.model import io as mio
    mio.save(w.model, str(out), "clustal")
    assert out.exists()
    # stats/identity/plot dialogs construct
    from neoedit.ui.dialogs.misc_dialogs import StatsDialog, IdentityDialog, PlotDialog
    StatsDialog(w.model, range(5), w); IdentityDialog(w.model, range(5), w); PlotDialog(w.model, range(5), w)
    # translate to new window
    from neoedit.analysis import translate as T
    aa = T.translate_aligned(w.model.rows[0].seq, 2, 0)
    assert aa.startswith("MF")
    w.model.dirty = False
    w.close()


def test_right_click_actions(app):
    from neoedit.ui.main_window import MainWindow
    from PySide6.QtCore import QPoint
    w = MainWindow(); w.show(); w.open_path(EXAMPLE)
    v = w.view; v.resize(800, 400)
    rect = v.cell_rect(1, 4)
    pos = rect.center()
    before = [r.seq for r in w.model.rows]
    # insert gap in clicked sequence (row 1)
    w._set_right_click("ins_sel")
    QTest.mouseClick(v.viewport(), Qt.RightButton, Qt.NoModifier, pos)
    assert w.model.rows[1].seq[4] == "-" and w.model.rows[0].seq == before[0]
    # delete it again
    w._set_right_click("del_sel")
    QTest.mouseClick(v.viewport(), Qt.RightButton, Qt.NoModifier, pos)
    assert w.model.rows[1].seq == before[1]
    # insert gap in all unselected (everything except row 1)
    w._set_right_click("ins_other")
    QTest.mouseClick(v.viewport(), Qt.RightButton, Qt.NoModifier, pos)
    assert w.model.rows[1].seq == before[1]
    assert all(w.model.rows[i].seq[4] == "-" for i in (0, 2, 3, 4))
    w._set_right_click("del_other")
    QTest.mouseClick(v.viewport(), Qt.RightButton, Qt.NoModifier, pos)
    assert [r.seq for r in w.model.rows] == before
    w.model.dirty = False
    w.close()
