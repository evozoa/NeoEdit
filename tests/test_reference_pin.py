"""Pinned reference row, insertion carets, and .bio round-trip of the pin."""
import os
import pytest

from neoedit.model import AlignmentModel, SequenceRow
from neoedit.model import io as mio
from neoedit.model import bioedit_format as B

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _model():
    #                       ref has gaps at cols 4-6 -> an insertion in the samples
    return AlignmentModel([SequenceRow("refseq", "ACGT---ACGTACGT"),
                           SequenceRow("sample1", "ACGTTTTACGTACGT"),
                           SequenceRow("sample2", "ACGTTT-ACGTACGT"),
                           SequenceRow("maskrow", "****---********")])


def test_set_reference_and_mask_detection():
    m = _model()
    assert m.ref_row == 0 and m.reference().name == "refseq"
    m.set_reference(1)
    assert m.ref_row == 1 and m.reference().name == "sample1"
    assert m.is_mask_row(3) and not m.is_mask_row(0)
    # deleting rows keeps the pin on the same sequence
    m.set_reference(2)
    m.remove_rows([0])
    assert m.reference().name == "sample2"


def test_bio_roundtrip_of_pin(tmp_path):
    m = _model()
    m.set_reference(2)
    m.mask_row = 3
    p = tmp_path / "x.bio"
    mio.save(m, str(p), "bio")
    _rows, info = B.read_bio(str(p))
    assert info["numbering_index"] == 2 and info["mask_index"] == 3
    m2 = mio.load(str(p))
    assert m2.ref_row == 2 and m2.mask_row == 3
    # a numbering mask pointing at a synthetic mask row must not become the reference
    m3 = _model(); m3.set_reference(3)
    p2 = tmp_path / "y.bio"
    mio.save(m3, str(p2), "bio")
    m4 = mio.load(str(p2))
    assert m4.ref_row == 0


def test_insertion_carets_and_projection():
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from neoedit.ui.main_window import MainWindow
    w = MainWindow(); w.show()
    m = _model(); m.remove_rows([3]); m.dirty = False
    w._set_model(m)
    ins = w.insertion_carets()
    assert len(ins) == 1
    ref_pos, ncols, nseq, col0 = ins[0]
    assert (ref_pos, ncols, nseq, col0) == (4, 3, 2, 4)     # after ref base 4, 3 columns, both samples
    # pinning a sample removes the insertion (it now has the bases)
    w.model.set_reference(1)
    w._proj = None; w._proj_row = -1
    assert w.insertion_carets() == []
    # coordinates follow the pin
    w.model.set_reference(0); w._proj = None; w._proj_row = -1
    assert w.proj().ref_len == 12 and w.ref_ungapped() == "ACGTACGTACGT"
    w.model.set_reference(1); w._proj = None; w._proj_row = -1; w._ref_ungapped = None
    assert w.proj().ref_len == 15
    # clicking a caret scrolls the grid to those columns
    w.model.set_reference(0); w._proj = None; w._proj_row = -1; w._ref_ungapped = None
    w._goto_columns(4, 7)
    sel = w.view.selection()
    assert sel and sel[2] == 4 and sel[3] == 6
    w.model.dirty = False; w.close()


def test_circular_from_genbank():
    HERE = os.path.dirname(__file__)
    gb = os.path.join(HERE, "..", "examples", "mito", "NC_002333.gb")
    m = mio.load(gb, "genbank")
    assert m.circular is True
