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


def test_rotate_model_and_annotation():
    from neoedit.model import AlignmentModel, SequenceRow, Feature
    seq = "AAAACCCCGGGGTTTT"                  # 16 bp
    m = AlignmentModel([SequenceRow("a", seq), SequenceRow("b", seq.lower())], "dna")
    m.features = [Feature(0, 2, 6, 1, "ORF", "x"), Feature(0, 12, 16, 1, "ORF", "y")]
    m.rotate(4)                               # old column 4 -> 0
    assert m.rows[0].seq == "CCCCGGGGTTTTAAAA" and m.rows[1].seq == "ccccggggttttaaaa"
    fx = [f for f in m.features if f.label == "x"]; fy = [f for f in m.features if f.label == "y"]
    assert [(f.start, f.end) for f in fy] == [(8, 12)]
    assert sorted((f.start, f.end) for f in fx) == [(0, 2), (14, 16)]      # x now crosses the origin: two pieces
    m.rotate(14)                              # pieces meet again and re-join
    assert [(f.start, f.end) for f in m.features if f.label == "x"] == [(0, 4)]
    assert m.rows[0].seq == "AACCCCGGGGTTTTAA"
    m.rotate(0, flip=True)
    assert m.rows[0].seq == str(Seq("AACCCCGGGGTTTTAA").reverse_complement())
    assert [(f.start, f.end, f.strand) for f in m.features if f.label == "x"] == [(12, 16, -1)]
    # annotation: a gene across the origin is whole again once the origin sits at its start
    ann = GA.Annotation(); ann.lengths["c"] = 1000; ann.circular["c"] = True
    t = GA.Transcript("t", "g", 900, 1050, 1, "CDS", [(900, 1000), (1000, 1050)], [(900, 1000), (1000, 1050)])
    ann.add_gene(GA.Gene("g", "g", "c", 900, 1050, 1, "CDS", [t], attrs={"wraps_origin": "true"}))
    ann.add_gene(GA.Gene("h", "h", "c", 100, 200, -1, "CDS", [GA.Transcript("th", "h", 100, 200, -1, "CDS", [(100, 200)], [(100, 200)])]))
    ann.finalize()
    GA.rotate_annotation(ann, "c", 900, 1000)
    g = ann.find("g")[0]; h = ann.find("h")[0]
    assert (g.start, g.end, "wraps_origin" in g.attrs) == (0, 150, False) and g.transcripts[0].cds == [(0, 100), (100, 150)]
    assert (h.start, h.end) == (200, 300)
    # flip: minus-strand h reads forward when the origin is put at its 5' end (= its old end) and the strand is flipped
    GA.rotate_annotation(ann, "c", 300, 1000, flip=True)
    h = ann.find("h")[0]; g = ann.find("g")[0]
    assert (h.start, h.end, h.strand) == (0, 100, 1)
    assert (g.start, g.end, g.strand) == (150, 300, -1)          # (0,150)+ -> rotated to 700 -> mirrored
    assert GA.largest_gap_origin([(900, 1050), (100, 200)], 1000) == 900       # widest gap 200..900 -> next feature starts at 900
    assert GA.largest_gap_origin([(0, 1000)], 1000) is None


def test_set_origin_gui():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")
    import sys
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    from neoedit.ui.main_window import MainWindow
    from neoedit.ui.dialogs.misc_dialogs import OriginDialog
    if not os.path.exists(RCRS):
        pytest.skip("rCRS example missing")
    w = MainWindow(); w.resize(1100, 700); w.show()
    w.open_path(RCRS); app.processEvents()
    seq0 = w.model.rows[0].seq
    L = len(seq0)
    co = w.annotation.find("COX1")[0]
    cox_start = co.start
    # rotate so that COX1 is split across the origin, then use Set origin at COX1 -> whole again
    w.rotate_origin(cox_start + 300)
    assert w.model.rows[0].seq == seq0[cox_start + 300:] + seq0[:cox_start + 300]
    co = w.annotation.find("COX1")[0]
    assert co.end > L and co.attrs.get("wraps_origin") == "true"
    assert [f for f in w._translation_regions(0, 0, 10) if f[4] == "COX1"]            # visible right after the origin...
    d = OriginDialog(L, [(g.name, g.start, g.end, g.strand) for g in w.annotation.genes_by_seq["NC_012920.1"]], 0, w)
    idx = next(i for i in range(d.feat.count()) if d.feat.itemText(i).startswith("COX1 "))
    d.r_feat.setChecked(True); d.feat.setCurrentIndex(idx)
    pos, flip = d.values()
    assert pos == co.start and not flip
    w.rotate_origin(pos, flip)
    co = w.annotation.find("COX1")[0]
    assert (co.start, co.end) == (0, 1542) and "wraps_origin" not in co.attrs
    assert w.model.rows[0].seq[:12] == "ATGTTCGCCGAC"                                  # COX1 begins at position 1
    regs = [r for r in w._translation_regions(0, 0, 2000) if r[4] == "COX1"]
    assert len(regs) == 1 and regs[0][:2] == (0, 1542) and regs[0][5] == 0
    # default suggestion: tRNA-Phe (vertebrate mitogenome convention) is preselected
    d2 = OriginDialog(L, [(g.name, g.start, g.end, g.strand) for g in w.annotation.genes_by_seq["NC_012920.1"]], 0, w)
    assert d2.feat.currentText().upper().startswith("TRNF")
    # minus-strand ND6 at the origin with a flip reads forward from position 1
    nd6 = w.annotation.find("ND6")[0]
    d3 = OriginDialog(L, [(g.name, g.start, g.end, g.strand) for g in w.annotation.genes_by_seq["NC_012920.1"]], 0, w)
    j = next(i for i in range(d3.feat.count()) if d3.feat.itemText(i).startswith("ND6 "))
    d3.r_feat.setChecked(True); d3.feat.setCurrentIndex(j)
    assert d3.flip.isChecked()
    pos, flip = d3.values(); assert flip and pos == nd6.end % L
    w.rotate_origin(pos, flip)
    nd6 = w.annotation.find("ND6")[0]
    assert (nd6.start, nd6.end, nd6.strand) == (0, 525, 1)
    assert str(Seq(w.model.rows[0].seq[:525]).translate(table=2)).startswith("MMYALF")
    w.model.dirty = False; w.close()


def test_primer_design_circular():
    from neoedit.analysis import primers as P, primer_design as PD
    if not os.path.exists(RCRS):
        pytest.skip("rCRS example missing")
    seq = mio.load(RCRS).rows[0].seq; L = len(seq)
    rot = seq[6500:] + seq[:6500]; dbl = rot + rot
    # target across the origin: every product must cross it; primers verified on the ring
    pairs = P.design_primers(rot, target=(L - 100, 200), product_range=((200, 400),), num_return=5, circular=True)
    assert pairs and all(p.crosses_origin and p.circular_len == L for p in pairs)
    p = pairs[0]
    ls, le = p.left.span; rs, re_ = p.right.span
    assert dbl[ls:le] == p.left.seq and str(Seq(dbl[rs:re_]).reverse_complement()) == p.right.seq
    assert p.product_size == re_ - ls and p.left.start < L <= p.right.start
    assert p.right.pieces(L) == [(rs - L, re_ - L)] and "across" not in p.right.pos_text(L)
    feats = P.pair_to_features(p, 0, None, 1)
    assert [(f.start, f.end) for f in feats] == [(ls, le), (rs - L, re_ - L)]
    # the same target on a linear template is refused with a hint, not a Primer3 crash
    with pytest.raises(ValueError, match="circular"):
        P.design_primers(rot, target=(L - 100, 200), product_range=((200, 400),), circular=False)
    # no target: candidates are unique and reported once, on the first turn
    p2 = P.design_primers(rot, product_range=((100, 300),), num_return=6, circular=True)
    assert len(p2) == 6 and all(x.left.start < L for x in p2) and len({(x.left.start, x.right.start) for x in p2}) == 6
    # a primer site straddling the origin: two pieces, positions read 'a-b (across origin)'
    pr = P.Primer("ACGT" * 5, L - 8, 20, 60.0, 50.0, 0, 0, 0)
    assert pr.pieces(L) == [(L - 8, L), (0, 12)] and pr.pos_text(L) == f"{L - 7}-12 (across origin)"
    # scoring a site given as pieces reads the ring in order (forward and reverse)
    site = rot[L - 8:] + rot[:12]
    hit = PD.score_primer(site, [rot], ["a"], L - 8, L, 1, pieces=[(L - 8, L), (0, 12)])[0]
    assert hit.mismatches == 0 and hit.covered
    rhit = PD.score_primer(str(Seq(site).reverse_complement()), [rot], ["a"], L - 8, L, -1, pieces=[(L - 8, L), (0, 12)])[0]
    assert rhit.mismatches == 0
    # alignment design across the origin: pieces, in-silico table, features on every row
    rnd = random.Random(1)
    def mut(s, n):
        s = list(s)
        for _ in range(n):
            i = rnd.randrange(len(s)); s[i] = rnd.choice("ACGT")
        return "".join(s)
    rows = [rot, mut(rot, 30), mut(rot, 30)]
    res = PD.design_on_alignment(rows, ["a", "b", "c"], template_row=0, include_rows=[0, 1, 2], product_range=((150, 400),),
                                 num_return=5, circular=True, min_conservation=0.8, primer3_kwargs=dict(target=(L - 100, 200)))
    evs = list(res)
    assert evs and all(e.pair.crosses_origin for e in evs)
    e = evs[0]
    assert e.left_pieces[0][0] < L and e.right_pieces[0][1] <= L and e.cols_text("right") == f"{e.right_pieces[0][0] + 1:,}-{e.right_pieces[0][1]:,}"
    tbl = PD.insilico_table(e)
    assert tbl[0]["fwd_mm"] == 0 and tbl[0]["rev_mm"] == 0 and tbl[0]["amplifies"]


def test_primer_dialogs_follow_topology():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")
    import sys
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    from neoedit.ui.main_window import MainWindow
    if not os.path.exists(RCRS):
        pytest.skip("rCRS example missing")
    w = MainWindow(); w.show(); w.open_path(RCRS); app.processEvents()
    assert w.model.circular
    w.primer_design(); d = w._children[-1]; assert d.circ_cb.isChecked(); d.close()
    w.model.duplicate_rows([0])                       # across-alignment design needs two sequences
    w.design_primers(); d2 = w._children[-1]; assert d2.circ_cb.isChecked(); d2.close()
    w.set_topology(False)
    w.primer_design(); d3 = w._children[-1]; assert not d3.circ_cb.isChecked(); d3.close()
    w.model.dirty = False; w.close()
