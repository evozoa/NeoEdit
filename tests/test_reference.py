import pytest
from neoedit.genome.projection import RefProjection
from neoedit.model import SequenceRow
from neoedit.analysis.external import mafft_add, find_mafft


def test_projection_identity():
    p = RefProjection("ACGTACGT")
    assert p.identity and p.ref_to_col(5) == 5 and p.col_to_ref(5) == 5
    assert p.span_to_cols(2, 4) == (2, 4)


def test_projection_gapped():
    #      col: 0123456789
    #      seq: A-CG--TAC-
    #      ref: 0 12  345
    p = RefProjection("A-CG--TAC-")
    assert p.ref_len == 6 and not p.identity
    assert [p.ref_to_col(u) for u in range(6)] == [0, 2, 3, 6, 7, 8]
    assert [p.col_to_ref(c) for c in range(10)] == [0, 1, 1, 2, 3, 3, 3, 4, 5, 6]
    assert p.span_to_cols(1, 3) == (2, 4)      # residues C,G
    assert p.span_to_cols(3, 6) == (6, 9)      # T,A,C
    assert p.span_to_cols(0, 6) == (0, 9)
    # round trip for every residue
    for u in range(6):
        assert p.col_to_ref(p.ref_to_col(u)) == u


@pytest.mark.skipif(not find_mafft(), reason="mafft not installed")
def test_mafft_add_keeplength():
    ref = [SequenceRow("ref", "ATGGCCATTGTACTGAGCCATCCGTATGCAAGCTTGGACTACGGCTAA")]
    # barcode = middle slice of ref with one substitution and an insertion (to be trimmed)
    bc = ref[0].seq[10:40].replace("CCG", "CCA", 1)
    bc = bc[:15] + "GGGG" + bc[15:]
    rows = mafft_add(ref, [SequenceRow("bc1", bc)])
    assert rows[0].seq == ref[0].seq                      # reference untouched
    assert len(rows) == 2 and len(rows[1].seq) == len(ref[0].seq)   # keeplength
    assert rows[1].seq.startswith("-")                    # leading gap where barcode absent


@pytest.mark.skipif(not find_mafft(), reason="mafft not installed")
def test_reference_workflow_gui(tmp_path):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from neoedit.ui.main_window import MainWindow
    from neoedit.model import io as mio
    from neoedit.genome import annotations as GA
    HERE = os.path.dirname(__file__)
    gb = os.path.join(HERE, "..", "examples", "mito", "NC_002333.gb")
    coi = os.path.join(HERE, "..", "examples", "mito", "coi_barcodes.fasta")
    w = MainWindow(); w.resize(1200, 700); w.show()
    # open reference
    model = mio.load(gb, "genbank"); model.features = []; model.dirty = False
    ann = GA.load_genbank(gb)
    w._set_model(model); w._enter_reference_mode(ann)
    assert w.genome_panel.isVisible() and w.proj().identity
    L = w.proj().ref_len
    assert L > 16000 and ann.count() > 30
    coi_genes = [g for g in ann.find("COX1") + ann.find("COI") + ann.find("co1")]
    if not coi_genes:
        coi_genes = [g for genes in ann.genes_by_seq.values() for g in genes if "COX1" in g.name.upper() or "COI" in g.name.upper()]
    assert coi_genes, [g.name for genes in ann.genes_by_seq.values() for g in genes][:20]
    cox1 = coi_genes[0]
    # anchor barcodes with MAFFT --add --keeplength
    from neoedit.analysis.external import mafft_add
    new = mio.load(coi)
    rows = mafft_add(w.model.rows, new.rows)
    w.model.begin_batch("Anchor")
    for r in rows[1:] if len(rows) == len(new.rows) + 1 else rows[w.model.nrows:]:
        w.model.add_row(r)
    w.model.end_batch()
    assert w.model.nrows == 1 + len(new.rows)
    assert all(len(r.seq) == L for r in w.model.rows)     # keeplength: no new columns
    # barcodes should land inside COX1
    bc = w.model.rows[1].seq
    first = len(bc) - len(bc.lstrip("-")); last = len(bc.rstrip("-"))
    assert cox1.start - 50 <= first and last <= cox1.end + 50, (first, last, cox1.start, cox1.end)
    # variation is nonzero inside the barcode region, zero outside coverage
    var = w._variation(first, first + 200)
    assert var and any(v > 0 for v in var)
    out = w._variation(0, 100)
    assert out is not None and all(v == 0 for v in out[:50])
    # projection stays valid after inserting a gap into the reference row
    w.model.insert_gaps(0, first + 10, 3)
    assert not w.proj().identity
    feats = w._genome_features(0, first, first + 300)
    assert feats
    # gene view coordinate of cox1 start maps to a column with a real residue
    c = w.proj().ref_to_col(cox1.start)
    assert w.model.rows[0].seq[c] not in "-.~"
    w.model.undo()
    w.model.dirty = False; w.close()
