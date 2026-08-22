"""Importing an annotated GenBank record must show its chromosome map.

Regression guard for the importer losing the genome view when it stopped routing an
import into an empty alignment through open_path (commit 6ec37b7). Also pins the rule
that annotations belong to the *reference* sequence: the map and the gene models follow
whichever row is pinned, and switching the reference switches the features on screen.
"""
import os
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZEBRAFISH = os.path.join(ROOT, "examples", "mito", "NC_002333.gb")
HUMAN = os.path.join(ROOT, "examples", "mito", "NC_012920.gb")
PLAIN = os.path.join(ROOT, "examples", "cox1_demo.fasta")

pytestmark = pytest.mark.skipif(not os.path.exists(ZEBRAFISH), reason="example GenBank records not present")


@pytest.fixture
def win():
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from neoedit.ui.main_window import MainWindow
    w = MainWindow(); w.show()
    yield w
    w.model.dirty = False        # closeEvent would otherwise raise a modal save prompt
    w.close()


def _ref_name(w):
    return w.model.rows[w.model.ref_row].accession      # names carry the full definition line; the accession is the key


def _shown(w):
    """Genes currently addressed by the genome view (i.e. the reference's)."""
    if not w.annotation or not w.genome_contig:
        return 0
    return len(w.annotation.genes_by_seq.get(w.genome_contig, []))


def test_import_into_empty_alignment_shows_the_map(win):
    win.import_path(ZEBRAFISH)
    assert win.genome_panel.isVisible()
    assert win.genome_contig == "NC_002333.2"
    assert _ref_name(win) == "NC_002333.2"
    assert _shown(win) == 38


def test_import_into_populated_alignment_pins_the_record(win):
    win.open_path(PLAIN)
    assert not win.genome_panel.isVisible()
    win.import_path(ZEBRAFISH)
    assert win.genome_panel.isVisible()
    assert _ref_name(win) == "NC_002333.2"
    assert _shown(win) == 38


def test_features_follow_the_reference(win):
    win.import_path(ZEBRAFISH)
    win.import_path(HUMAN)
    assert _ref_name(win) == "NC_012920.1"       # newest annotated record takes the pin
    human = _shown(win)
    assert human > 0

    zf = next(i for i, r in enumerate(win.model.rows) if r.accession == "NC_002333.2")
    win.view.cur_row = zf
    win.pin_reference()
    assert win.genome_contig == "NC_002333.2"
    assert _shown(win) == 38

    # the grid provider hands features to the reference row and to no other row
    prov = win.view.feature_provider
    assert prov is not None
    for i in range(win.model.nrows):
        got = prov(i, 0, 2000)
        assert (len(got) > 0) == (i == win.model.ref_row), f"row {i} disagrees with the reference"


def test_annotated_import_does_not_duplicate_features_in_the_grid(win):
    win.import_path(ZEBRAFISH)
    # the annotation drives the display, so the per-row copies are dropped
    assert [f for f in win.model.features if f.row == win.model.ref_row] == []


def test_unannotated_import_leaves_the_reference_alone(win):
    win.import_path(ZEBRAFISH)
    win.import_path(PLAIN)
    assert _ref_name(win) == "NC_002333.2"
    assert win.genome_panel.isVisible()


def test_opening_a_plain_file_clears_a_stale_map(win):
    win.open_path(ZEBRAFISH)
    assert win.genome_panel.isVisible()
    win.open_path(PLAIN)
    assert not win.genome_panel.isVisible()
    assert win.annotation is None
    assert win.genome_contig is None
