import os
from Bio import SeqIO
from neoedit.analysis import genbank_export as GX
from neoedit.analysis.orf_finder import find_orfs
from neoedit.model import io as mio

HERE = os.path.dirname(__file__)
GB = os.path.join(HERE, "..", "examples", "mito", "NC_002333.gb")


def test_translate_and_feature():
    orf = find_orfs("ATGAAACCCGGGTAA", table=1, min_aa=3, both_strands=False, allow_partial=False)[0]
    a = GX.ORFAnnotation(orf, name="testpep", product="test peptide", table=1)
    f = GX.orf_to_feature(a)
    assert f.type == "CDS"
    assert f.qualifiers["transl_table"] == ["1"]
    assert f.qualifiers["translation"] == ["MKPG"]        # stop trimmed
    assert any("cytoplasmic ribosomes" in n for n in f.qualifiers["note"])
    # same ORF read with the vertebrate mito code
    a2 = GX.ORFAnnotation(orf, table=2, cytoplasmic=False)
    f2 = GX.orf_to_feature(a2)
    assert f2.qualifiers["transl_table"] == ["2"]
    assert not any("cytoplasmic" in n for n in f2.qualifiers["note"])


def test_merge_into_real_record(tmp_path):
    src = GX.load_source_record(GB)
    seq = str(src.seq)
    n_before = len(src.features)
    # find candidate MDP-like ORFs with the STANDARD code inside the 16S rRNA region
    rrna = [f for f in src.features if f.type == "rRNA"]
    assert rrna
    big = max(rrna, key=lambda f: len(f))
    s, e = int(big.location.start), int(big.location.end)
    orfs = find_orfs(seq[s:e], table=1, min_aa=20, both_strands=False, allow_partial=False, start_mode="atg")
    assert orfs, "expected at least one standard-code ORF inside the rRNA gene"
    host = big.qualifiers.get("product", ["rRNA"])[0]
    anns = []
    for o in orfs[:3]:
        o.start += s; o.end += s          # back to record coordinates
        name, note = GX.suggest_mdp_name(o, host)
        anns.append(GX.ORFAnnotation(o, name=name, product=name, table=1, cytoplasmic=True, note=note))
    rec = GX.build_record(seq, anns, source=src, organism="Danio rerio")
    assert len(rec.features) == n_before + len(anns)
    path = tmp_path / "out.gb"
    GX.write_genbank(rec, str(path))
    # re-read: both translation tables coexist, canonical genes untouched
    back = SeqIO.read(str(path), "genbank")
    tables = {q for f in back.features for q in f.qualifiers.get("transl_table", [])}
    assert "1" in tables and "2" in tables
    cds1 = [f for f in back.features if f.qualifiers.get("transl_table") == ["1"]]
    assert len(cds1) == len(anns)
    # translations are correct under their own table
    for f in cds1:
        nt = str(f.extract(back.seq))
        assert GX.translate_orf(nt, 1) == f.qualifiers["translation"][0]
    # a canonical mito CDS still translates under table 2
    cox1 = [f for f in back.features if f.type == "CDS" and "COX1" in str(f.qualifiers.get("gene", ""))]
    if cox1:
        assert cox1[0].qualifiers["transl_table"] == ["2"]
    assert "NeoEdit" in back.annotations.get("comment", "")


def test_build_without_source():
    seq = "ATG" + "AAACCCGGG" * 5 + "TAA"
    orf = find_orfs(seq, table=1, min_aa=5, both_strands=False, allow_partial=False)[0]
    rec = GX.build_record(seq, [GX.ORFAnnotation(orf, name="pep1", table=1)],
                          record_id="TEST1", organism="Testus fishus", topology="circular")
    text = GX.dumps_genbank(rec)
    assert "TEST1" in text and "transl_table=1" in text and "Testus fishus" in text
    assert "source" in text and "CDS" in text


def test_orf_dialog_genbank_export(tmp_path, monkeypatch):
    """End-to-end: ORF finder -> named MDP ORFs -> GenBank with both codes."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from neoedit.ui.main_window import MainWindow
    from neoedit.ui.dialogs.orf_dialog import ORFFinderDialog
    from neoedit.genome import annotations as GA
    w = MainWindow(); w.show()
    model = mio.load(GB, "genbank"); model.features = []; model.dirty = False
    ann = GA.load_genbank(GB)
    w._set_model(model); w._enter_reference_mode(ann)
    d = ORFFinderDialog(w.model, [0], w, 1)          # standard code = MDP-style ORFs
    d.table.setCurrentIndex(d.table.findData(1))
    d.min_aa.setValue(20); d.both.setChecked(False)
    d.start_mode.setCurrentIndex(d.start_mode.findData("atg"))
    d.run()
    assert d.table_w.rowCount() > 0
    # ORFs inside rRNA genes get a host-gene column and an auto-suggested name
    hosts = [d.table_w.item(i, 2).text() for i in range(d.table_w.rowCount())]
    assert any(h not in ("—", "") for h in hosts)
    named = [o for o in d.orfs if o.name]
    assert named
    # export without opening the modal options dialog
    out = tmp_path / "mdp.gb"
    monkeypatch.setattr("neoedit.ui.dialogs.orf_dialog.GenBankExportOptions.exec", lambda self: True)
    monkeypatch.setattr("neoedit.ui.dialogs.orf_dialog.GenBankExportOptions.values", lambda self: dict(
        table=1, default_table=2, cytoplasmic=True, feature_type="CDS", use_source=True,
        source_path=GB, record_id="NC_002333_MDP", organism="Danio rerio", topology="circular"))
    monkeypatch.setattr("neoedit.ui.dialogs.orf_dialog.QFileDialog.getSaveFileName",
                        staticmethod(lambda *a, **k: (str(out), "")))
    d.export_genbank()
    assert out.exists()
    back = SeqIO.read(str(out), "genbank")
    tabs = {q for f in back.features for q in f.qualifiers.get("transl_table", [])}
    assert {"1", "2"} <= tabs, tabs
    mdp = [f for f in back.features if f.qualifiers.get("transl_table") == ["1"]]
    assert mdp and all("cytoplasmic ribosomes" in " ".join(f.qualifiers.get("note", [])) for f in mdp)
    # canonical genes survived untouched
    assert any(f.type == "rRNA" for f in back.features)
    assert len([f for f in back.features if f.type == "CDS"]) > 13
    w.model.dirty = False; w.close()
