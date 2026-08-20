"""NCBI / Ensembl importers. Offline by default (fake fetcher); set NEOEDIT_NET_TESTS=1 for live checks."""
import json
import os
import urllib.parse

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from neoedit.remote import RemoteError
from neoedit.remote import ncbi as N, ensembl as E
from neoedit.remote.http import write_download, safe_filename

HERE = os.path.dirname(__file__)
LIVE = os.environ.get("NEOEDIT_NET_TESTS") == "1"


# ---------------------------------------------------------------- helpers
class FakeFetch:
    """Routes URLs to canned responses; records every URL it was asked for."""

    def __init__(self, routes):
        self.routes = routes          # list of (substring, body or callable(url) -> body)
        self.urls = []

    def __call__(self, url, **kw):
        self.urls.append(url)
        for key, body in self.routes:
            if key in url:
                if callable(body):
                    body = body(url)
                if isinstance(body, Exception):
                    raise body
                return body.encode() if isinstance(body, str) else body
        raise RemoteError(f"unrouted {url}")


def test_parse_helpers():
    assert N.parse_ids("NC_012920.1, MN996528.1\n>KX353935 NC_012920.1") == ["NC_012920.1", "MN996528.1", "KX353935"]
    assert E.parse_region("17:43,044,295-43,170,245") == ("17", 43044295, 43170245, 1)
    assert E.parse_region("MT:1..16569:-1") == ("MT", 1, 16569, -1)
    assert E.parse_region("chrX:100-50:-") == ("chrX", 50, 100, -1)
    with pytest.raises(RemoteError):
        E.parse_region("chr1 100 200")
    assert E.looks_like_id("ENSG00000012048.23") and E.looks_like_id("ENSDARG00000000001") and E.looks_like_id("ENST00000357654")
    assert not E.looks_like_id("BRCA1") and not E.looks_like_id("MT-CO1")
    assert E.normalize_species("Homo sapiens") == "homo_sapiens"
    assert safe_filename("17:43044295-43170245(+)") == "17_43044295-43170245_+"
    assert N._unspace("Error: F a i l e d  t o  u n d e r s t a n d  i d :  N O P E") == "Error: Failed to understand id: NOPE"


def test_write_download_dedupes(tmp_path):
    p1 = write_download("A", str(tmp_path), "x", ".gb")
    p2 = write_download("A", str(tmp_path), "x", ".gb")       # identical → same file
    p3 = write_download("B", str(tmp_path), "x", ".gb")       # different → numbered
    assert p1 == p2 and p3.endswith("x_2.gb") and open(p3).read() == "B"


# ---------------------------------------------------------------- Ensembl record building
def _ov(ft, id_, s, e, strand, **kw):
    d = {"feature_type": ft, "id": id_, "start": s, "end": e, "strand": strand, "seq_region_name": "1"}
    d.update(kw)
    return d


def _synthetic():
    """A 300 bp + region with: a + gene (2 CDS exons), a − gene clipped at the region's right edge."""
    import random
    rnd = random.Random(7)
    seq = "".join(rnd.choice("ACGT") for _ in range(300))
    feats = [
        _ov("gene", "G1", 1011, 1120, 1, external_name="ALPHA", biotype="protein_coding", description="alpha gene"),
        _ov("transcript", "T1", 1011, 1120, 1, Parent="G1", biotype="protein_coding", is_canonical=1, version=3, external_name="ALPHA-201"),
        _ov("exon", "X1", 1011, 1040, 1, Parent="T1"), _ov("exon", "X2", 1081, 1120, 1, Parent="T1"),
        _ov("cds", "T1", 1021, 1040, 1, Parent="T1", phase=0, protein_id="P1"),
        _ov("cds", "T1", 1081, 1110, 1, Parent="T1", phase="1", protein_id="P1"),
        _ov("gene", "G2", 1250, 1400, -1, external_name="BETA", biotype="lncRNA"),
        _ov("transcript", "T2", 1250, 1400, -1, Parent="G2", biotype="lncRNA", version=1),
        _ov("exon", "X3", 1250, 1280, -1, Parent="T2"), _ov("exon", "X4", 1350, 1400, -1, Parent="T2"),
    ]
    return seq, feats


def test_build_record_plus_and_roundtrip(tmp_path):
    seq, feats = _synthetic()
    rec = E.build_record(seq, "1", 1001, 1300, feats, species="homo_sapiens", assembly="GRCh38", release=116)
    types = [(f.type, int(f.location.start), int(f.location.end), f.location.strand) for f in rec.features]
    assert types[0] == ("source", 0, 300, 1)
    assert ("gene", 10, 120, 1) in types and ("mRNA", 10, 120, 1) in types and ("CDS", 20, 110, 1) in types
    cds = [f for f in rec.features if f.type == "CDS"][0]
    assert [(int(p.start), int(p.end)) for p in cds.location.parts] == [(20, 40), (80, 110)]
    assert cds.qualifiers["codon_start"] == [1] and cds.qualifiers["protein_id"] == ["P1"]
    assert cds.qualifiers["transcript_id"] == ["T1.3"] and "transl_table" not in cds.qualifiers
    # the − lncRNA runs past the region end (1400 > 1300): clipped and marked partial at its 5' (high) end
    nc = [f for f in rec.features if f.type == "ncRNA"][0]
    assert nc.qualifiers["ncRNA_class"] == ["lncRNA"] and nc.location.strand == -1
    assert str(nc.location) == "[249:>280](-)"       # its 5' exon lies outside the region: partial
    # survives a GenBank round-trip and is what NeoEdit's loader will see
    p = E.save_record(rec, str(tmp_path))
    back = next(SeqIO.parse(p, "genbank"))
    assert str(back.seq) == seq and len(back.features) == len(rec.features)
    assert "Ensembl 116" in back.description and back.annotations["organism"] == "Homo sapiens"


def test_build_record_minus_orientation_and_phase():
    """Oriented to −: sequence reverse-complemented, features mirrored, biological order kept."""
    seq, feats = _synthetic()
    plus = E.build_record(seq, "1", 1001, 1300, feats)
    minus = E.build_record(seq, "1", 1001, 1300, feats, strand=-1)
    assert str(minus.seq) == str(Seq(seq).reverse_complement())
    cds_p = [f for f in plus.features if f.type == "CDS"][0]
    cds_m = [f for f in minus.features if f.type == "CDS"][0]
    assert cds_m.location.strand == -1
    # extracting the CDS from either orientation yields the same coding sequence
    assert str(cds_p.extract(plus.seq)) == str(cds_m.extract(minus.seq))
    # clip the + CDS at its 5' end inside exon 2 (phase 1): 1084 removes 3 bases of exon 2 -> codon_start = ((1-3)%3)+1 = 2
    rec = E.build_record(seq, "1", 1084, 1300, feats)
    cds = [f for f in rec.features if f.type == "CDS"][0]
    assert cds.qualifiers["codon_start"] == [2] and str(cds.location).startswith("[<0:")
    # region start 1085 removes 4 bases -> ((1-4)%3)+1 = 1
    rec = E.build_record(seq, "1", 1085, 1300, feats)
    assert [f for f in rec.features if f.type == "CDS"][0].qualifiers["codon_start"] == [1]
    # MT chromosome: vertebrate mito code on CDS, organelle on source
    rec = E.build_record(seq, "MT", 1001, 1300, feats)
    assert [f for f in rec.features if f.type == "CDS"][0].qualifiers["transl_table"] == ["2"]
    assert rec.features[0].qualifiers["organelle"] == ["mitochondrion"]


# ---------------------------------------------------------------- Ensembl client (offline)
def test_ensembl_resolve_and_fetch_offline(tmp_path):
    fasta = ">chromosome:GRCh38:MT:5901:5960:1\nCTGATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATAC\n"
    routes = [
        ("e116.rest.ensembl.org/info/software", RemoteError("archive down")),
        ("rest.ensembl.org/info/software", json.dumps({"release": 116})),
        ("/info/species", json.dumps({"species": [
            {"name": "homo_sapiens", "display_name": "Human", "common_name": "human", "assembly": "GRCh38"},
            {"name": "danio_rerio", "display_name": "Zebrafish", "common_name": "zebrafish", "assembly": "GRCz11"}]})),
        ("/lookup/symbol/homo_sapiens/MT-CO1", json.dumps({"id": "ENSG00000198804", "object_type": "Gene", "species": "homo_sapiens",
                                                            "seq_region_name": "MT", "start": 5904, "end": 7445, "strand": 1,
                                                            "display_name": "MT-CO1", "biotype": "protein_coding", "assembly_name": "GRCh38"})),
        ("/lookup/id/ENSP00000354499", json.dumps({"id": "ENSP00000354499", "object_type": "Translation", "Parent": "ENST00000361624", "start": 5904, "end": 7445})),
        ("/lookup/id/ENST00000361624", json.dumps({"id": "ENST00000361624", "object_type": "Transcript", "species": "homo_sapiens",
                                                   "seq_region_name": "MT", "start": 5904, "end": 7445, "strand": 1, "Parent": "ENSG00000198804"})),
        ("/sequence/region/homo_sapiens/MT:5901..5960:1", fasta),
        ("/overlap/region/homo_sapiens/MT:5901..5960", json.dumps([
            _ov("gene", "ENSG00000198804", 5904, 7445, 1, external_name="MT-CO1", biotype="protein_coding"),
            _ov("transcript", "ENST00000361624", 5904, 7445, 1, Parent="ENSG00000198804", biotype="protein_coding", version=2),
            _ov("exon", "ENSE1", 5904, 7445, 1, Parent="ENST00000361624"),
            _ov("cds", "ENST00000361624", 5904, 7445, 1, Parent="ENST00000361624", phase=0, protein_id="ENSP00000354499")])),
        ("/sequence/id/ENSG00000198804", ">ENST00000361624.2\nATGTTCGCCGACCGTTGA\n"),
    ]
    ff = FakeFetch(routes)
    c = E.EnsemblClient(fetch=ff)
    srv, rel = c.resolve()
    assert srv == E.LIVE_SERVER and rel == 116 and "release 116" in c.release_note()
    assert ff.urls[0].startswith("https://e116.rest.ensembl.org")           # archive tried first
    sp = c.species()
    assert [s.name for s in sp] == ["homo_sapiens", "danio_rerio"] and sp[0].label() == "Human — homo_sapiens [GRCh38]"
    lk = c.lookup("MT-CO1", "Homo sapiens")
    assert (lk.id, lk.seq_region_name, lk.start, lk.strand) == ("ENSG00000198804", "MT", 5904, 1)
    lk2 = c.lookup("ENSP00000354499", "homo_sapiens")            # translation -> parent transcript
    assert lk2.id == "ENST00000361624" and lk2.start == 5904
    rec = c.fetch_genomic("homo_sapiens", "MT", 5901, 5960, gene=lk)
    assert rec.id == "MT-CO1" and len(rec.seq) == 60
    cds = [f for f in rec.features if f.type == "CDS"][0]
    assert str(cds.location) == "[3:>60](+)" and cds.qualifiers["transl_table"] == ["2"]
    assert str(cds.extract(rec.seq).translate(table=2)) == "MFADRWLFSTNHKDIGTLY"
    recs = c.fetch_sequences("ENSG00000198804", "cds", "MT-CO1", "homo_sapiens")
    assert recs[0].id == "ENST00000361624.2" and "MT-CO1 cds" in recs[0].description
    p = E.save_records(recs, str(tmp_path), "MT-CO1_cds")
    assert p.endswith("MT-CO1_cds.fasta")
    # a pinned server that serves another release is reported, not silently accepted
    c2 = E.EnsemblClient(server="https://rest.ensembl.org", fetch=FakeFetch([("/info/software", json.dumps({"release": 117}))]))
    assert "not available" in c2.release_note()
    # region too big for gene models
    with pytest.raises(RemoteError):
        c.overlap_region("homo_sapiens", "1", 1, 6_000_000)


# ---------------------------------------------------------------- NCBI client (offline)
def test_ncbi_offline(tmp_path):
    gb = "LOCUS       NC_012920              16569 bp    DNA     circular PRI 01-JAN-2020\n//\n"
    routes = [
        ("esearch.fcgi", json.dumps({"esearchresult": {"count": "19", "idlist": ["1", "2"]}})),
        ("esummary.fcgi", json.dumps({"result": {"uids": ["1", "2"],
                                                 "1": {"caption": "NC_012920", "accessionversion": "NC_012920.1", "title": "Homo sapiens mitochondrion",
                                                       "slen": 16569, "organism": "Homo sapiens", "topology": "circular"},
                                                 "2": {"error": "cannot get document summary"}}})),
        ("efetch.fcgi", lambda url: ("Error: F a i l e d  t o  u n d e r s t a n d  i d :  N O P E" if "NOPE" in url else gb)),
    ]
    ff = FakeFetch(routes)
    c = N.NCBIClient(email="me@example.org", fetch=ff)
    c._throttle = lambda: None
    assert c.search("nuccore", "x[Organism]", 5) == (19, ["1", "2"])
    assert "email=me%40example.org" in ff.urls[-1] and "tool=neoedit" in ff.urls[-1]
    sums = c.summaries("nuccore", ["1", "2"])
    assert len(sums) == 1 and sums[0].accession == "NC_012920.1" and sums[0].length == 16569 and sums[0].topology == "circular"
    path, text = c.download("nuccore", ["NC_012920.1"], "gb", str(tmp_path), seq_start=577, seq_stop=647)
    assert path.endswith("NC_012920.1_577-647.gb") and N.count_records(text, "gb") == 1
    q = urllib.parse.parse_qs(urllib.parse.urlsplit(ff.urls[-1]).query)
    assert q["rettype"] == ["gbwithparts"] and q["seq_start"] == ["577"] and q["seq_stop"] == ["647"]
    path, _ = c.download("protein", ["YP_003024028.1"], "gb", str(tmp_path))
    assert "rettype=gp" in ff.urls[-1]
    with pytest.raises(RemoteError, match="Failed to understand id: NOPE"):
        c.fetch("nuccore", ["NOPE_000000"], "fasta")


# ---------------------------------------------------------------- GUI wiring
def test_import_dialog_wiring(tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")
    import sys
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    from neoedit.ui.main_window import MainWindow
    from neoedit.ui.dialogs.import_dialog import ImportDialog
    w = MainWindow(); w.show()
    assert w.a_import_remote.shortcut().toString() == "Ctrl+Shift+I"
    w.import_remote()
    dlg = w._import_dialog
    assert isinstance(dlg, ImportDialog) and dlg.isVisible()
    dlg.autoload = False
    dlg.dir_edit.setText(str(tmp_path))
    # NCBI job is built from the text box + ticked search rows + sub-range
    dlg.n_ids.setPlainText("NC_012920.1 MN996528.1")
    job, busy = dlg._ncbi_job()
    assert "2 record(s)" in busy
    dlg.n_ids.setPlainText("NC_012920.1"); dlg.n_from.setValue(577); dlg.n_to.setValue(647)
    job, busy = dlg._ncbi_job(); assert "1 record(s)" in busy
    dlg.n_ids.setPlainText("")
    with pytest.raises(RemoteError):
        dlg._ncbi_job()
    # Ensembl jobs: gene / region / type validation
    dlg.tabs.setCurrentIndex(1)
    dlg.e_gene.setText("MT-CO1"); job, busy = dlg._ensembl_job(); assert "MT-CO1" in busy
    dlg.e_r_region.setChecked(True); dlg.e_region.setText("MT:1-16569:-1"); job, busy = dlg._ensembl_job(); assert "MT:1-16569" in busy
    dlg.e_type.setCurrentIndex(2)
    with pytest.raises(RemoteError):
        dlg._ensembl_job()
    # a finished fetch hands the file to the main window: "add" appends rows, "open" replaces the alignment
    example = os.path.join(HERE, "..", "examples", "cox1_demo.fasta")
    dlg.r_add.setChecked(True)
    dlg._fetched((example, "done"))
    assert w.model.nrows == 5
    dlg._fetched((example, "done"))
    assert w.model.nrows == 10
    w.model.dirty = False
    dlg.r_open.setChecked(True)
    gbp = os.path.join(HERE, "..", "examples", "mito", "NC_012920_MDP.gb")
    if os.path.exists(gbp):
        dlg._fetched((gbp, "done"))
        assert w.model.nrows == 1 and w.model.circular
    dlg.close()
    w.close()


# ---------------------------------------------------------------- live (opt-in)
@pytest.mark.skipif(not LIVE, reason="set NEOEDIT_NET_TESTS=1 to hit NCBI/Ensembl")
def test_live_ensembl_and_ncbi(tmp_path):
    c = E.EnsemblClient()
    srv, rel = c.resolve()
    assert rel == 116, c.release_note()
    lk = c.lookup("MT-ND6", "homo_sapiens")
    rec = c.fetch_genomic(lk.species, lk.seq_region_name, lk.start, lk.end, strand=lk.strand, gene=lk)
    cds = [f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"] == ["MT-ND6"]][0]
    aa = str(cds.extract(rec.seq).translate(table=2))
    assert aa.startswith("MMYALF") and aa.endswith("*") and "*" not in aa[:-1]
    n = N.NCBIClient()
    path, text = n.download("nuccore", ["NC_012920.1"], "fasta", str(tmp_path), seq_start=577, seq_stop=647)
    assert text.startswith(">NC_012920.1:577-647") and N.count_records(text, "fasta") == 1
