"""True circular support: features and ORFs across the origin, topology choice, exports."""
import io
import os
import random

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from neoedit.analysis import orf_finder as OF, genbank_export as GX
from neoedit.analysis.translate import translate_region
from neoedit.genome import annotations as GA
from neoedit.model import io as mio

HERE = os.path.dirname(__file__)
RCRS = os.path.join(HERE, "..", "examples", "mito", "NC_012920_MDP.gb")


def _molecule():
    """756 bp with a 150-aa + ORF and a 60-aa − ORF at known places."""
    rnd = random.Random(5)
    rand = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    core = "ATG" + "".join(rnd.choice(["GCT", "GGC", "AAA", "CCT", "GAT", "TTC", "TGG", "CAT", "GTA", "AGC"]) for _ in range(149)) + "TAA"
    minus = str(Seq("ATG" + "".join(rnd.choice(["GCT", "GGC", "AAA", "CCT", "GAT"]) for _ in range(59)) + "TAG").reverse_complement())
    mol = rand(50) + core + rand(40) + minus + rand(30)
    return mol, core, minus


def test_orf_finder_across_origin():
    mol, core, minus = _molecule()
    L = len(mol)
    lin = OF.find_orfs(mol, table=1, min_aa=50, start_mode="atg", allow_partial=False)
    plus = [o for o in lin if o.strand > 0][0]
    # rotate so the + ORF straddles the origin
    rot = mol[200:] + mol[:200]
    circ = OF.find_orfs(rot, table=1, min_aa=50, start_mode="atg", circular=True)
    o = [o for o in circ if o.strand > 0 and o.length_aa == 150][0]
    assert o.wraps and o.start == (50 - 200) % L and o.end - o.start == len(core) and o.aa == plus.aa
    assert o.extra["length"] == L
    # linear scan of the rotated molecule cannot see it; circular reports it exactly once
    assert not any(x.length_aa == 150 for x in OF.find_orfs(rot, table=1, min_aa=50, start_mode="atg", allow_partial=False))
    assert sum(1 for x in circ if x.start == o.start and x.strand > 0) == 1
    # pieces + phases translate in frame on both sides of the origin
    f1, f2 = o.to_features(None)
    assert (f1.start, f1.end, f1.data["wrap_part"]) == (o.start, L, 1) and (f2.start, f2.end, f2.data["wrap_part"]) == (0, o.end - L, 2)
    aa = translate_region(rot, f1.start, f1.end, 1, 1, f1.data["phase"]) + translate_region(rot, f2.start, f2.end, 1, 1, f2.data["phase"])
    assert aa.rstrip("*") == o.aa.rstrip("*")
    # minus strand across the origin: 5' piece is the one after the origin
    mstart = 50 + len(core) + 40; mid = mstart + len(minus) // 2
    rot2 = mol[mid:] + mol[:mid]
    o2 = [x for x in OF.find_orfs(rot2, table=1, min_aa=50, start_mode="atg", circular=True) if x.strand < 0 and x.length_aa == 60][0]
    assert o2.wraps
    g1, g2 = o2.to_features(None)
    n2 = g2.end - g2.start
    a_head = translate_region(rot2, g2.start, g2.end, 1, -1, g2.data["phase"])
    a_tail = translate_region(rot2, g1.start, g1.end, 1, -1, g1.data["phase"])
    k = n2 // 3
    assert a_head == o2.aa[:k] and a_tail == o2.aa[k + (1 if n2 % 3 else 0):]
    # exports
    gff = OF.orfs_to_gff([o], "m")
    assert gff.count("\tCDS\t") == 2 and f"\t{o.start + 1}\t{L}\t" in gff and "\t1\t" + str(o.end - L) + "\t" in gff
    tbl = OF.orfs_to_genbank_table([o2], "m")
    assert f"{o2.end - L}\t1\tCDS" in tbl and f"{L}\t{o2.start + 1}" in tbl
    # GenBank: join across the origin, round-trips and translates
    rec = GX.build_record(rot2, [GX.ORFAnnotation(o2, table=1)], record_id="ROT", topology="circular", default_table=1)
    txt = GX.dumps_genbank(rec)
    assert f"complement(join({o2.start + 1}..{L},1..{o2.end - L}))" in txt
    back = next(SeqIO.parse(io.StringIO(txt), "genbank"))
    cds = [f for f in back.features if f.type == "CDS"][0]
    assert str(cds.extract(back.seq).translate(table=1)) == o2.aa
    # a linear molecule is unaffected
    assert [(x.start, x.end) for x in OF.find_orfs(mol, table=1, min_aa=50, start_mode="atg", allow_partial=False, circular=False)] == \
           [(x.start, x.end) for x in lin]


def test_genbank_features_across_origin_stay_whole():
    if not os.path.exists(RCRS):
        pytest.skip("rCRS example missing")
    ann = GA.load_genbank(RCRS)
    d = ann.find("D-loop")[0]
    assert (d.start, d.end) == (16023, 16569 + 576) and d.attrs["wraps_origin"] == "true"
    assert d in ann.overlapping("NC_012920.1", 16500, 16569) and d in ann.overlapping("NC_012920.1", 0, 10)
    assert d not in ann.overlapping("NC_012920.1", 600, 16000)
    # the grid's per-row feature copies: two pieces, never a whole-genome span
    m = mio.load(RCRS)
    pieces = sorted((f.start, f.end) for f in m.features if f.label.startswith("D-loop"))
    assert pieces == [(0, 576), (16023, 16569)] and m.circular and m.topology_known
    assert max(f.end - f.start for f in m.features) < 5000


def test_gff_two_row_gene_and_is_circular(tmp_path):
    gff = ("##gff-version 3\n##sequence-region chrM 1 1000\n"
           "chrM\tx\tregion\t1\t1000\t.\t+\t.\tID=r;Is_circular=true\n"
           "chrM\tx\tgene\t901\t1000\t.\t+\t.\tID=dl;Name=Dloop\n"
           "chrM\tx\tgene\t1\t50\t.\t+\t.\tID=dl;Name=Dloop\n"
           "chrM\tx\tgene\t200\t400\t.\t-\t.\tID=g2;Name=other\n")
    p = tmp_path / "c.gff3"; p.write_text(gff)
    ann = GA.load_gff(str(p))
    assert ann.is_circular("chrM") and ann.length_of("chrM") == 1000
    d = ann.find("Dloop")[0]
    assert (d.start, d.end) == (900, 1050) and d.transcripts[0].exons == [(900, 1000), (1000, 1050)]
    assert [g.name for g in ann.overlapping("chrM", 0, 10)] == ["Dloop"]
    assert GA.split_span(900, 1050, 1000) == [(900, 1000), (0, 50)]
    assert GA.fmt_span(900, 1050, 1000) == "901-50 (across origin)"
    # the same file without the topology stays as two pieces until the user declares it circular
    p2 = tmp_path / "l.gff3"; p2.write_text(gff.replace("chrM\tx\tregion\t1\t1000\t.\t+\t.\tID=r;Is_circular=true\n", ""))
    ann2 = GA.load_gff(str(p2))
    assert ann2.find("Dloop")[0].end == 1000
    ann2.set_topology("chrM", 1000, True)
    assert ann2.find("Dloop")[0].end == 1050


def test_topology_choice_gui(tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")
    import sys
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    from neoedit.ui.main_window import MainWindow
    from neoedit.model import AlignmentModel, SequenceRow
    if not os.path.exists(RCRS):
        pytest.skip("rCRS example missing")
    w = MainWindow(); w.resize(1100, 700); w.show()
    # a GenBank record declares its topology: adopted, shown, no question asked
    w.open_path(RCRS); app.processEvents()
    assert w.model.circular and w.model.topology_known and not w.topo_bar.isVisible()
    assert w.topo_combo.currentData() is True and "circular" in w.lbl_topo.text()
    assert w.annotation.is_circular("NC_012920.1")
    # region view draws the D-loop at both ends; the grid overlay gets both pieces
    w.genome_panel.set_window(0, 1000); app.processEvents(); w.genome_panel.region.grab()
    assert any(g.name.startswith("D-loop") for _r, g, _t in w.genome_panel.region._gene_hits)
    assert [(f.start, f.end) for f in w._genome_features(0, 0, 600) if "D-loop" in f.label] == [(0, 576)]
    assert [(f.start, f.end) for f in w._genome_features(0, 16000, 16569) if "D-loop" in f.label] == [(16023, 16569)]
    # a plain sequence: topology not set; the Topology box sets it
    seq = w.model.rows[0].seq
    w.model.dirty = False
    m = AlignmentModel([SequenceRow("rot", seq[6500:] + seq[:6500])], "dna"); w._set_model(m)
    assert not m.topology_known and w.topo_combo.currentData() is None
    w.topo_combo.setCurrentIndex(2); app.processEvents()
    assert m.circular and m.topology_known
    w.orf_finder(); d = w._children[-1]
    assert d.circ_cb.isChecked()
    d.run(); assert any(o.wraps for o in d.orfs); d.close()
    # annotation of unknown topology onto a FASTA asks once; answering closes the bar and flags the annotation
    gff = tmp_path / "u.gff3"
    gff.write_text("##gff-version 3\nNC_012920.1\tx\tgene\t16024\t16569\t.\t+\t.\tID=dl;Name=Dloop\nNC_012920.1\tx\tgene\t1\t576\t.\t+\t.\tID=dl;Name=Dloop\n")
    w.model.dirty = False
    m3 = AlignmentModel([SequenceRow("NC_012920.1", seq)], "dna"); w._set_model(m3)
    w.annotation = None; w.genome_contig = None
    w.load_annotation_path(str(gff)); w._enter_reference_mode(w.annotation); app.processEvents()
    assert w.topo_bar.isVisible() and not m3.topology_known
    assert w.annotation.find("Dloop")[0].end == 16569
    w.set_topology(True); app.processEvents()
    assert not w.topo_bar.isVisible() and m3.circular and m3.topology_known
    assert w.annotation.find("Dloop")[0].end == 16569 + 576
    # translation regions carry a phase for the piece after the origin
    assert MainWindow._wrap_pieces(16000, 16569 + 100, 1, 16569) == [(16000, 16569, 0), (0, 100, (3 - 569 % 3) % 3)]
    assert MainWindow._wrap_pieces(16000, 16569 + 100, -1, 16569) == [(16000, 16569, (3 - 100 % 3) % 3), (0, 100, 0)]
    w.model.dirty = False; w.close()
