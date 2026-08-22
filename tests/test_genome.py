import os
from neoedit.genome.fasta_index import IndexedFasta
from neoedit.genome.annotations import load_gff, load_bed, load_paf, pack_lanes


def test_indexed_fasta(tmp_path):
    p = tmp_path / "g.fa"
    seq1 = "ACGT" * 25   # 100 bp
    seq2 = "GGGGCCCCAAAATTTT" * 3
    p.write_text(">chr1 desc\n" + "\n".join(seq1[i:i + 30] for i in range(0, 100, 30)) + "\n>chr2\n" + seq2 + "\n")
    fa = IndexedFasta(str(p))
    assert os.path.exists(str(p) + ".fai")
    assert fa.contigs() == [("chr1", 100), ("chr2", 48)]
    assert fa.fetch("chr1", 0, 10) == seq1[:10]
    assert fa.fetch("chr1", 28, 35) == seq1[28:35]     # crosses a line break
    assert fa.fetch("chr1", 95, 200) == seq1[95:]
    assert fa.fetch("chr2") == seq2
    # re-open uses the index
    fa2 = IndexedFasta(str(p)); assert fa2.fetch("chr1", 97, 100) == seq1[97:]


GFF = """##gff-version 3
chr1\tx\tgene\t101\t1000\t.\t+\t.\tID=gene1;Name=abc1;gene_biotype=protein_coding
chr1\tx\tmRNA\t101\t1000\t.\t+\t.\tID=tx1;Parent=gene1;Name=abc1-201
chr1\tx\texon\t101\t300\t.\t+\t.\tID=e1;Parent=tx1
chr1\tx\texon\t501\t1000\t.\t+\t.\tID=e2;Parent=tx1
chr1\tx\tCDS\t201\t300\t.\t+\t0\tID=c1;Parent=tx1
chr1\tx\tCDS\t501\t800\t.\t+\t2\tID=c2;Parent=tx1
chr1\tx\tgene\t2001\t2500\t.\t-\t.\tID=gene2;Name=xyz2
chr1\tx\ttRNA\t2001\t2500\t.\t-\t.\tID=t2;Parent=gene2
chr1\tx\texon\t2001\t2500\t.\t-\t.\tParent=t2
chr2\tx\tgene\t1\t50\t.\t+\t.\tID=gene3
"""


def test_gff3(tmp_path):
    p = tmp_path / "a.gff3"; p.write_text(GFF)
    ann = load_gff(str(p))
    assert ann.count() == 3
    g = ann.find("abc1")[0]
    assert (g.start, g.end, g.strand) == (100, 1000, 1)
    t = g.transcripts[0]
    assert t.exons == [(100, 300), (500, 1000)] and t.cds == [(200, 300), (500, 800)]
    assert t.utrs() == [(100, 200), (800, 1000)]
    assert [x.name for x in ann.overlapping("chr1", 0, 150)] == ["abc1"]
    assert [x.name for x in ann.overlapping("chr1", 1500, 2100)] == ["xyz2"]
    assert ann.overlapping("chr1", 1200, 1900) == []
    assert ann.find("gene3")[0].transcripts[0].exons == [(0, 50)]
    ann2 = load_gff(str(p), only_seqid="chr2"); assert ann2.count() == 1


def test_gtf_and_bed(tmp_path):
    gtf = ('chr1\tx\tgene\t101\t1000\t.\t+\t.\tgene_id "g1"; gene_name "foo";\n'
           'chr1\tx\ttranscript\t101\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
           'chr1\tx\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
           'chr1\tx\tCDS\t201\t300\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n')
    p = tmp_path / "a.gtf"; p.write_text(gtf)
    ann = load_gff(str(p)); g = ann.find("foo")[0]
    assert g.transcripts[0].exons == [(100, 300)] and g.transcripts[0].cds == [(200, 300)]
    bed = "chr1\t10\t100\tfeatA\t0\t-\t20\t90\t0\t2\t30,40\t0,50\n"
    p = tmp_path / "a.bed"; p.write_text(bed)
    ann = load_bed(str(p)); g = ann.find("featA")[0]
    assert g.strand == -1 and g.transcripts[0].exons == [(10, 40), (60, 100)] and g.transcripts[0].cds == [(20, 40), (60, 90)]


def test_paf_and_lanes(tmp_path):
    paf = "q1\t1000\t0\t500\t+\tt1\t5000\t100\t600\t450\t500\t60\nq1\t1000\t600\t900\t-\tt2\t5000\t0\t300\t200\t300\t10\n"
    p = tmp_path / "a.paf"; p.write_text(paf)
    b = load_paf(str(p), min_len=100)
    assert len(b) == 2 and b[1].strand == -1 and abs(b[0].identity - 0.9) < 1e-9
    assert len(load_paf(str(p), min_len=100, min_mapq=30)) == 1
    lanes = pack_lanes([(0, 10, "a"), (5, 15, "b"), (12, 20, "c")])
    assert [[o for _, _, o in l] for l in lanes] == [["a", "c"], ["b"]]


def test_genome_panel_gui(tmp_path):
    import os, random
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from neoedit.ui.main_window import MainWindow
    random.seed(1)
    L = 60000
    seq = "".join(random.choice("ACGT") for _ in range(L))
    fa = tmp_path / "g.fa"
    fa.write_text(">chrA\n" + "\n".join(seq[i:i + 80] for i in range(0, L, 80)) + "\n>chrB\n" + "ACGT" * 500 + "\n")
    gff = tmp_path / "g.gff3"
    lines = ["##gff-version 3"]
    for k in range(12):
        s = 1000 + k * 4500; e = s + 3000; st = "+" if k % 2 == 0 else "-"
        lines += [f"chrA\tx\tgene\t{s}\t{e}\t.\t{st}\t.\tID=g{k};Name=gene{k};gene_biotype=protein_coding",
                  f"chrA\tx\tmRNA\t{s}\t{e}\t.\t{st}\t.\tID=t{k};Parent=g{k};Name=gene{k}-201",
                  f"chrA\tx\texon\t{s}\t{s + 500}\t.\t{st}\t.\tParent=t{k}", f"chrA\tx\texon\t{s + 1500}\t{e}\t.\t{st}\t.\tParent=t{k}",
                  f"chrA\tx\tCDS\t{s + 200}\t{s + 500}\t.\t{st}\t0\tParent=t{k}", f"chrA\tx\tCDS\t{s + 1500}\t{e - 300}\t.\t{st}\t0\tParent=t{k}"]
    gff.write_text("\n".join(lines) + "\n")
    paf = tmp_path / "s.paf"
    paf.write_text("chrA\t60000\t0\t30000\t+\tNC_1\t100000\t5000\t35000\t28000\t30000\t60\nchrA\t60000\t31000\t59000\t-\tNC_2\t90000\t0\t28000\t25000\t28000\t60\n")
    w = MainWindow(); w.resize(1100, 700); w.show()
    w.open_genome_path(str(fa), "chrA")
    assert w.model.nrows == 1 and len(w.model.rows[0].seq) == L and w.genome_panel.isVisible()
    w.load_annotation_path(str(gff))
    assert w.annotation.count() == 12
    from neoedit.genome.annotations import load_paf
    w.synteny_blocks = load_paf(str(paf), min_len=100); w.genome_panel.set_synteny(w.synteny_blocks)
    w.genome_panel.set_window(0, 30000); app.processEvents()
    # grid overlay features via provider
    feats = w._genome_features(0, 0, 10000)
    assert any(f.type == "CDS" for f in feats) and any(f.type == "UTR" for f in feats)
    # go to a gene by name -> focus moves the grid
    w.genome_panel.region_edit.setText("gene5"); w.genome_panel._goto_text(); app.processEvents()
    assert w.view.horizontalScrollBar().value() > 0
    # switch contig
    w._genome_load_contig("chrB"); assert w.model.rows[0].name == "chrB"
    w._genome_load_contig("chrA")
    # open region in new window with features
    w._genome_open_region(1000, 5000)
    child = w._children[-1]
    assert child.model.rows[0].seq == seq[1000:5000] and child.model.features
    w.genome_panel.grab(); w.view.viewport().repaint()
    child.close(); w.model.dirty = False; w.close()


def test_circular_view(tmp_path):
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from neoedit.ui.circular_view import CircularView, gc_series, gene_color
    from neoedit.genome.annotations import load_genbank
    from neoedit.model import io as mio
    HERE = os.path.dirname(__file__)
    gb = os.path.join(HERE, "..", "examples", "mito", "NC_002333.gb")
    gc, skew = gc_series("GGGGCCCCAAAATTTT", 2)
    assert gc == [1.0, 0.0] and skew[0] == 0.0
    gc, skew = gc_series("GGGGGGGGCCCCCCCC", 2)
    assert skew == [1.0, -1.0]
    ann = load_genbank(gb)
    seq = mio.load(gb, "genbank").rows[0].seq
    v = CircularView(); v.resize(600, 600)
    v.set_data("NC_002333.2", len(seq), ann, lambda s, e: seq[s:e], "zebrafish mtDNA")
    v.set_focus(6400, 8000)
    v.grab()                      # paints without error and populates hit paths
    assert v._hits and any(g.name.upper() == "COX1" for _p, g in v._hits)
    # colors distinguish feature classes
    genes = {g.name.upper(): g for _p, g in v._hits}
    assert gene_color(genes["COX1"]) != gene_color(genes["TRNA"]) if "TRNA" in genes else True
    png = tmp_path / "c.png"; svg = tmp_path / "c.svg"
    v.export_image(str(png), 500); v.export_image(str(svg), 500)
    assert png.stat().st_size > 5000 and svg.stat().st_size > 5000
    assert "<svg" in svg.read_text()[:400]


def test_origin_spanning_and_mdp_features():
    """rCRS: the D-loop wraps the origin and must not swallow the genome; MDP CDSs
    carry transl_table=1 and are flagged cytoplasmic."""
    import os
    from neoedit.genome.annotations import load_genbank
    HERE = os.path.dirname(__file__)
    gb = os.path.join(HERE, "..", "examples", "mito", "NC_012920_MDP.gb")
    if not os.path.exists(gb):
        import pytest; pytest.skip("rCRS+MDP example not built")
    ann = load_genbank(gb)
    dl = [g for gs in ann.genes_by_seq.values() for g in gs if g.name.startswith("D-loop")]
    # one feature, kept whole in unwrapped coordinates: join(16024..16569,1..576) -> [16023, 17145)
    assert len(dl) == 1 and dl[0].attrs.get("wraps_origin") == "true"
    d = dl[0]
    assert (d.start, d.end, len(d)) == (16023, 16569 + 576, 1122) and ann.lengths["NC_012920.1"] == 16569
    assert ann.is_circular("NC_012920.1")
    # found from both sides of the origin, not by a window in the middle of the genome
    assert d in ann.overlapping("NC_012920.1", 16500, 16569)
    assert d in ann.overlapping("NC_012920.1", 0, 100)
    assert d not in ann.overlapping("NC_012920.1", 5000, 6000)
    from neoedit.genome.annotations import split_span, fmt_span
    assert split_span(d.start, d.end, 16569) == [(16023, 16569), (0, 576)]
    assert fmt_span(d.start, d.end, 16569) == "16,024-576 (across origin)"
    mdp = sorted(g.name for gs in ann.genes_by_seq.values() for g in gs if g.cytoplasmic)
    assert mdp == ["MOTS-c", "SHLP1", "SHLP2", "SHLP3", "SHLP4", "SHLP5", "SHLP6", "humanin"]
    # humanin sits inside MT-RNR2 and uses the standard code
    hn = [g for gs in ann.genes_by_seq.values() for g in gs if g.name == "humanin"][0]
    assert (hn.start, hn.end, hn.strand) == (2632, 2704, 1)
    assert hn.attrs["transl_table"] == "1"


def test_minor_features_hidden_by_default(tmp_path):
    """A bare misc_feature (e.g. the rCRS placeholder at 3107) is a real annotation
    but not a gene, so the region view hides it unless asked."""
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from neoedit.ui.genome_panel import RegionView, MINOR_TYPES
    from neoedit.genome.annotations import load_genbank
    HERE = os.path.dirname(__file__)
    gb = os.path.join(HERE, "..", "examples", "mito", "NC_012920_MDP.gb")
    if not os.path.exists(gb):
        import pytest; pytest.skip("rCRS+MDP example not built")
    ann = load_genbank(gb)
    misc = [g for gs in ann.genes_by_seq.values() for g in gs if g.biotype == "misc_feature"]
    assert misc and len(misc[0]) == 1 and misc[0].start == 3106     # rCRS 3107 placeholder
    assert "misc_feature" in MINOR_TYPES
    v = RegionView(); v.resize(900, 300)
    v.length = 16569; v.seqid = misc[0].seqid; v.ann = ann
    v.set_window(3000, 3200)
    v.grab()
    drawn = {g.name for _r, g, _t in v._gene_hits}
    assert "misc_feature" not in drawn
    v.show_minor = True
    v.grab()
    assert "misc_feature" in {g.name for _r, g, _t in v._gene_hits}


def test_region_view_hideable():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from neoedit.ui.main_window import MainWindow
    HERE = os.path.dirname(__file__)
    gb = os.path.join(HERE, "..", "examples", "mito", "NC_012920_MDP.gb")
    if not os.path.exists(gb):
        import pytest; pytest.skip("rCRS example missing")
    w = MainWindow(); w.resize(1100, 700); w.show(); w.open_path(gb); app.processEvents()
    gp = w.genome_panel
    assert gp.isVisible() and gp.region_visible() and w.a_g_region.isChecked()
    w.a_g_region.trigger(); app.processEvents()                 # menu -> hidden, overview stays
    assert not gp.region_visible() and gp.overview.isVisible() and not gp.region_btn.isChecked()
    assert gp.preferred_height() < 120
    gp.region_btn.setChecked(True); app.processEvents()          # panel button -> shown, menu follows
    assert gp.region_visible() and w.a_g_region.isChecked()
    w.model.dirty = False; w.close()
