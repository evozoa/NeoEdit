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
    v.set_mode("edit"); v.set_edit_submode("insert")
    v.set_cursor(0, 5)
    QTest.keyClick(v, Qt.Key_N)
    assert w.model.rows[0].seq[5] == "N"
    w.undo()
    v.set_mode("slide")
    # selection + column gap
    v.select_region(0, 4, 10, 12)
    QTest.keyClick(v, Qt.Key_Space)
    assert all(r.seq[10] == "-" for r in w.model.rows)
    QTest.keyClick(v, Qt.Key_Backspace)  # removes the gap column again
    assert all(r.seq[10] != "-" for r in w.model.rows)
    # toggle translation overlay & color modes; repaint
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


def test_modes_slide_grab_downstream(app):
    from neoedit.ui.main_window import MainWindow
    from neoedit.model import AlignmentModel, SequenceRow
    w = MainWindow(); w.show()
    w._set_model(AlignmentModel([SequenceRow("a", "AC--GTAC"), SequenceRow("b", "ACGTGTAC")]))
    v = w.view; v.resize(800, 300)
    # --- Select/Slide: select AC in row 0 and drag right by 1 -> crunch gap ahead, open behind
    v.select_region(0, 0, 0, 1)
    start = v.cell_rect(0, 1).center()
    QTest.mousePress(v.viewport(), Qt.LeftButton, Qt.NoModifier, start)
    QTest.mouseMove(v.viewport(), v.cell_rect(0, 2).center())
    QTest.mouseRelease(v.viewport(), Qt.LeftButton, Qt.NoModifier, v.cell_rect(0, 2).center())
    assert w.model.rows[0].seq == "-AC-GTAC"
    w.undo(); assert w.model.rows[0].seq == "AC--GTAC"
    # --- Shift+drag = move whole downstream sequence (inserts gap at selection start)
    v.select_region(0, 0, 0, 1)
    QTest.mousePress(v.viewport(), Qt.LeftButton, Qt.ShiftModifier, v.cell_rect(0, 1).center())
    QTest.mouseMove(v.viewport(), v.cell_rect(0, 2).center())
    QTest.mouseRelease(v.viewport(), Qt.LeftButton, Qt.ShiftModifier, v.cell_rect(0, 2).center())
    assert w.model.rows[0].seq == "-AC--GTAC"
    w.undo()
    # --- Grab & Drag: grab the G (col 4) in row 0 and drag left by 2 over the gaps
    w._set_mode("grab")
    QTest.mousePress(v.viewport(), Qt.LeftButton, Qt.NoModifier, v.cell_rect(0, 4).center())
    QTest.mouseMove(v.viewport(), v.cell_rect(0, 2).center())
    QTest.mouseRelease(v.viewport(), Qt.LeftButton, Qt.NoModifier, v.cell_rect(0, 2).center())
    assert w.model.rows[0].seq == "ACG--TAC"
    w.undo()
    # grabbing a gap does nothing but start a selection
    QTest.mouseClick(v.viewport(), Qt.LeftButton, Qt.NoModifier, v.cell_rect(0, 2).center())
    assert w.model.rows[0].seq == "AC--GTAC"
    # --- Edit mode hides/shows insert/overwrite and Insert key toggles
    w._set_mode("edit")
    assert w.a_mode_insert.isVisible()
    QTest.keyClick(v, Qt.Key_Insert)
    assert v.edit_submode == "overwrite"
    w._set_mode("slide")
    assert not w.a_mode_insert.isVisible()
    # model-level downstream left move fails if non-gaps precede
    assert not w.model.move_downstream([1], 2, -1)
    w.model.dirty = False; w.close()


def test_bioedit_shortcuts(app):
    """Ctrl+Shift+R reverse-complements the selected rows, as in BioEdit (Ctrl+R kept as alias);
    no two actions may share a key binding."""
    from PySide6.QtGui import QAction
    from neoedit.ui.main_window import MainWindow
    w = MainWindow()
    w.show()
    w.open_path(EXAMPLE)
    assert [k.toString() for k in w.a_revcomp.shortcuts()] == ["Ctrl+Shift+R", "Ctrl+R"]
    seen = {}
    for a in w.findChildren(QAction):
        for k in a.shortcuts():
            if k.toString():
                seen.setdefault(k.toString(), []).append(a.text())
    assert not {k: v for k, v in seen.items() if len(v) > 1}
    w.activateWindow(); w.view.setFocus(); app.processEvents()
    rows = w.view.target_rows()
    seqs = [w.model.rows[r].seq for r in rows]
    from Bio.Seq import Seq
    rc = [str(Seq(s).reverse_complement()) for s in seqs]
    QTest.keyClick(w.view, Qt.Key_R, Qt.ControlModifier | Qt.ShiftModifier); app.processEvents()
    assert [w.model.rows[r].seq for r in rows] == rc
    QTest.keyClick(w.view, Qt.Key_R, Qt.ControlModifier); app.processEvents()   # alias undoes it again
    assert [w.model.rows[r].seq for r in rows] == seqs


def test_pinned_reference_strip_and_consensus_tool(app):
    """The header strip shows the pinned reference (not a consensus); the consensus is a tool."""
    from neoedit.ui.main_window import MainWindow
    from neoedit.ui.dialogs.misc_dialogs import ConsensusDialog
    from neoedit.model import AlignmentModel, SequenceRow
    w = MainWindow(); w.show()
    v = w.view
    assert v.pinned_row() == -1 and v.header_h == v.ruler_h          # empty: ruler only
    w.open_path(EXAMPLE)
    # all five rows on screen: the reference is visible, so no sticky copy (it would look like a 6th row)
    v.resize(900, 400); v.viewport().repaint()
    assert v.pinned_row() == -1 and v.header_h == v.ruler_h
    # shrink the viewport so only ~2 rows fit and scroll the reference off-screen -> the strip appears
    v.resize(900, v.ruler_h + 2 * v.row_h + 30); v.viewport().repaint()
    w.model.ref_row = 2
    v.verticalScrollBar().setValue(3)
    assert v.pinned_row() == 2 and v.header_h == v.ruler_h + v.cell_h
    w.a_pinned_ref.setChecked(False); w._toggle_pinned_ref(False)
    assert v.pinned_row() == -1
    w._toggle_pinned_ref(True)
    # clicking the strip's name jumps to the reference row (which scrolls back into view, hiding the strip)
    v.viewport().repaint()
    QTest.mouseClick(v.viewport(), Qt.LeftButton, Qt.NoModifier, QPoint(20, v.ruler_h + 4))
    assert v.sel_rows == {2} and v.verticalScrollBar().value() <= 2 and v.pinned_row() == -1
    # a lone sequence needs no pinned copy of itself
    w._set_model(AlignmentModel([SequenceRow("a", "ACGT")]))
    assert v.pinned_row() == -1
    # consensus tool: majority vs IUPAC, scope, add as row
    m = AlignmentModel([SequenceRow("a", "ACGT-A"), SequenceRow("b", "ACGA-A"), SequenceRow("c", "ACGTTA"), SequenceRow("d", "AC-TTG")])
    w._set_model(m)
    d = ConsensusDialog(m, [0, 1], 0.5, w)
    assert d.scope.currentData() == "sel" and d._cons == "ACGT-A"          # two rows: T/A tie at 50 % -> plurality keeps T
    d.plurality.setChecked(False); d.thr.setValue(0.6)
    assert d._cons == "ACGN-A"                                            # strict: the tie becomes N
    d.scope.setCurrentIndex(0); d.thr.setValue(0.5); d.plurality.setChecked(True)
    assert d._cons == "ACGTTA"                                            # all four: plurality
    d.thr.setValue(0.9); d.plurality.setChecked(False)
    assert d._cons == "ACGNTN"                                            # 90 %: only unanimous columns keep a residue
    d.method.setCurrentIndex(1); d.thr.setValue(0.1)
    assert d._cons == "ACGWTR"                                            # IUPAC: T/A -> W, A/G -> R
    d.name.setText("cons"); d.add_row()
    assert m.nrows == 5 and m.rows[-1].name == "cons" and m.rows[-1].seq == "ACGWTR"
    w.model.dirty = False; w.close()


def test_name_column_divider(app):
    """The name column is resizable: drag the divider, double-click to fit, menu to reset."""
    from neoedit.ui.main_window import MainWindow
    w = MainWindow(); w.show(); w.open_path(EXAMPLE)
    v = w.view; v.resize(900, 400); v.viewport().repaint()
    v.fit_name_width()                                   # start from the automatic width whatever was saved
    w0 = v.name_w
    assert v.on_divider(QPoint(w0, 50)) and not v.on_divider(QPoint(w0 + 40, 50))
    # drag the divider 60 px to the right
    QTest.mousePress(v.viewport(), Qt.LeftButton, Qt.NoModifier, QPoint(w0, 60))
    QTest.mouseMove(v.viewport(), QPoint(w0 + 30, 60)); QTest.mouseMove(v.viewport(), QPoint(w0 + 60, 60))
    QTest.mouseRelease(v.viewport(), Qt.LeftButton, Qt.NoModifier, QPoint(w0 + 60, 60))
    assert v.name_w == w0 + 60 and v.name_w_user == w0 + 60
    assert not v.sel_rows                                # the drag did not select a row
    # dragging cannot hide the grid
    v.set_name_width(5000); assert v.name_w <= v.viewport().width() - 120
    # reset / fit via the View actions
    w.a_fit_names_auto.trigger(); assert v.name_w_user is None and v.name_w == w0 <= 340
    w.model.rename(0, "A" * 80); w.a_fit_names.trigger()
    assert v.name_w > 340 and v.name_w >= v._fm.horizontalAdvance("A" * 80)
    w.model.dirty = False; w.close()


def test_gene_model_bars_off_by_default(app):
    """The reference's gene models are not drawn as bars under its residues unless asked."""
    import os
    from neoedit.ui.main_window import MainWindow
    gb = os.path.join(HERE, "..", "examples", "mito", "NC_012920_MDP.gb")
    if not os.path.exists(gb):
        pytest.skip("rCRS example missing")
    w = MainWindow(); w.show(); w.open_path(gb)
    v = w.view
    assert w.annotation and v.feature_provider is not None
    assert not v.show_gene_models and not w.a_gene_models.isChecked()
    assert v._features_at(0, 6000) == []                        # inside COX1, but no bar/tooltip
    w.a_gene_models.setChecked(True); w._toggle_gene_models(True)
    assert any(f.type == "CDS" for f in v._features_at(0, 6000))
    w._toggle_gene_models(False)
    w.model.dirty = False; w.close()
