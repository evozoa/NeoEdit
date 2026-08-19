from __future__ import annotations

from PySide6.QtGui import QFont
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QDialog, QFormLayout, QComboBox, QSpinBox, QCheckBox, QPushButton, QHBoxLayout,
                               QVBoxLayout, QLabel, QGroupBox, QTableWidgetItem, QMessageBox, QFileDialog,
                               QSplitter, QPlainTextEdit, QWidget, QLineEdit, QDialogButtonBox, QGridLayout)

from ...analysis import orf_finder as OF
from ...analysis import genbank_export as GX
from ...model.alignment import Feature
from .common import codon_table_combo, make_table, NumItem, save_text


class ORFFinderDialog(QDialog):
    """MitoFinder-style ORF finder with alternate genetic codes."""
    orfSelected = Signal(int, int, int)        # row, start(gapped), end(gapped)
    featuresReady = Signal(list)               # list[Feature]
    orfsFound = Signal(list, int)              # list[ORF], genetic code table

    def __init__(self, model, rows: list[int], parent=None, default_table=1):
        super().__init__(parent)
        self.setWindowTitle("ORF Finder")
        self.resize(900, 600)
        self.model = model
        self.rows = rows
        self.orfs: list[OF.ORF] = []
        self.annotation = getattr(parent, "annotation", None)
        self.source_path = None
        self._loading = False
        self.setModal(False)

        lay = QVBoxLayout(self)
        top = QHBoxLayout()
        params = QGroupBox("Parameters")
        form = QFormLayout(params)
        self.seq_combo = QComboBox()
        for r in rows:
            self.seq_combo.addItem(model.rows[r].name, r)
        self.seq_combo.addItem("All selected sequences", -1)
        self.table = codon_table_combo(default_table)
        self.min_aa = QSpinBox(); self.min_aa.setRange(1, 100000); self.min_aa.setValue(30)
        self.start_mode = QComboBox()
        self.start_mode.addItem("Start codons of genetic code (ATG + alternatives)", "table")
        self.start_mode.addItem("ATG only", "atg")
        self.start_mode.addItem("Any codon (stop-to-stop open frames)", "any")
        self.both = QCheckBox("Both strands"); self.both.setChecked(True)
        self.partial = QCheckBox("Include partial ORFs at sequence ends (incomplete stops)"); self.partial.setChecked(True)
        self.nested = QCheckBox("Report nested ORFs (alternative starts in same frame)")
        form.addRow("Sequence", self.seq_combo)
        form.addRow("Genetic code", self.table)
        form.addRow("Min length (aa)", self.min_aa)
        form.addRow("Start codons", self.start_mode)
        form.addRow("", self.both)
        form.addRow("", self.partial)
        form.addRow("", self.nested)
        top.addWidget(params, 1)
        lay.addLayout(top)

        btns = QHBoxLayout()
        self.run_btn = QPushButton("Find ORFs")
        self.run_btn.setDefault(True)
        self.run_btn.clicked.connect(self.run)
        self.annot_btn = QPushButton("Add as features")
        self.annot_btn.clicked.connect(self.add_features)
        self.exp_gff = QPushButton("Export GFF3")
        self.exp_gff.clicked.connect(lambda: self.export("gff"))
        self.exp_tbl = QPushButton("Export GenBank table")
        self.exp_tbl.clicked.connect(lambda: self.export("tbl"))
        self.exp_aa = QPushButton("Export proteins (FASTA)")
        self.exp_aa.clicked.connect(lambda: self.export("aa"))
        self.exp_nt = QPushButton("Export ORFs (nt FASTA)")
        self.exp_nt.clicked.connect(lambda: self.export("nt"))
        self.name_btn = QPushButton("Name from reference peptides…")
        self.name_btn.setToolTip("Load a FASTA of known peptides (e.g. humanin, MOTS-c, SHLP1-6) and name ORFs "
                                 "that are similar to them — homology by sequence, not by length")
        self.name_btn.clicked.connect(self.name_from_references)
        self.exp_gb = QPushButton("Export annotated GenBank…")
        self.exp_gb.setToolTip("Write a GenBank file containing the existing annotation plus these ORFs, "
                               "each with its own /transl_table (e.g. cytoplasmically-translated MDPs)")
        self.exp_gb.clicked.connect(self.export_genbank)
        for b in (self.run_btn, self.annot_btn, self.name_btn, self.exp_gb, self.exp_gff, self.exp_tbl, self.exp_aa, self.exp_nt):
            btns.addWidget(b)
        btns.addStretch()
        lay.addLayout(btns)

        split = QSplitter(Qt.Vertical)
        self.table_w = make_table(["Seq", "Name", "Within", "Start", "End", "Strand", "Frame", "Length (nt)",
                                   "Length (aa)", "Start codon", "Stop codon", "Partial", "Protein"])
        from PySide6.QtWidgets import QAbstractItemView
        self.table_w.setEditTriggers(QAbstractItemView.DoubleClicked | QAbstractItemView.SelectedClicked)
        self.table_w.itemChanged.connect(self._name_edited)
        self.table_w.itemSelectionChanged.connect(self._on_select)
        split.addWidget(self.table_w)
        self.detail = QPlainTextEdit(); self.detail.setReadOnly(True)
        f = self.detail.font(); f.setFamily("DejaVu Sans Mono"); f.setStyleHint(QFont.Monospace); self.detail.setFont(f)
        split.addWidget(self.detail)
        split.setSizes([400, 150])
        lay.addWidget(split, 1)
        self.status = QLabel("")
        lay.addWidget(self.status)

    def _host_gene(self, o: OF.ORF) -> str:
        """Name of the annotated gene the ORF sits inside (if an annotation is loaded)."""
        ann = self.annotation
        if ann is None:
            return ""
        seqid = getattr(self.parent(), "genome_contig", None) or self.model.rows[o.row].name
        try:
            genes = ann.overlapping(seqid, o.start, o.end)
        except Exception:
            return ""
        if not genes:
            return ""
        best = max(genes, key=lambda g: min(g.end, o.end) - max(g.start, o.start))
        return best.name

    def _name_edited(self, item):
        if self._loading or item.column() != 1:
            return
        i = item.data(Qt.UserRole + 1)
        if i is not None and 0 <= i < len(self.orfs):
            self.orfs[i].name = item.text().strip()

    def _targets(self):
        r = self.seq_combo.currentData()
        return self.rows if r == -1 else [r]

    def run(self):
        self.orfs = []
        self.table_w.setSortingEnabled(False)
        self.table_w.setRowCount(0)
        for r in self._targets():
            seq = self.model.rows[r].seq
            found = OF.find_orfs(seq, table=self.table.currentData(), min_aa=self.min_aa.value(),
                                 start_mode=self.start_mode.currentData(), both_strands=self.both.isChecked(),
                                 allow_partial=self.partial.isChecked(), nested=self.nested.isChecked(), row=r)
            self.orfs.extend(found)
        self._loading = True
        for i, o in enumerate(self.orfs):
            self.table_w.insertRow(i)
            partial = ("5'" if o.partial5 else "") + ("3'" if o.partial3 else "")
            host = self._host_gene(o)
            if not o.name:
                o.name, o.extra["note"] = GX.suggest_mdp_name(o, host)
            vals = [self.model.rows[o.row].name, o.name, host or "—", NumItem(o.start + 1), NumItem(o.end),
                    "+" if o.strand > 0 else "-", NumItem(o.frame), NumItem(o.length_nt), NumItem(o.length_aa),
                    o.start_codon, o.stop_codon or "(none)", partial or "-",
                    o.aa[:60] + ("…" if len(o.aa) > 60 else "")]
            for j, v in enumerate(vals):
                item = v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v))
                item.setData(Qt.UserRole + 1, i)
                if j != 1:
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.table_w.setItem(i, j, item)
        self._loading = False
        self.table_w.setSortingEnabled(True)
        self.table_w.sortItems(8, Qt.DescendingOrder)
        self.status.setText(f"{len(self.orfs)} ORF(s) found")
        self.orfsFound.emit(self.orfs, self.table.currentData())

    def _current_orf(self):
        items = self.table_w.selectedItems()
        if not items:
            return None
        return self.orfs[items[0].data(Qt.UserRole + 1)]

    def _on_select(self):
        o = self._current_orf()
        if not o:
            return
        gm = OF.gap_map(self.model.rows[o.row].seq)
        f = o.to_feature(gm)
        self.orfSelected.emit(o.row, f.start, f.end)
        self.detail.setPlainText(f"{self.model.rows[o.row].name}  {o.start + 1}-{o.end} ({'+' if o.strand > 0 else '-'}{o.frame})  "
                                 f"table {o.table}\n\nProtein ({o.length_aa} aa):\n{o.aa}\n\nNucleotide:\n{o.nt}")

    def add_features(self):
        feats = []
        for k, o in enumerate(self.orfs, 1):
            gm = OF.gap_map(self.model.rows[o.row].seq)
            if not o.name:
                o.name = f"ORF{k}"
            feats.append(o.to_feature(gm))
        if feats:
            self.featuresReady.emit(feats)
            self.status.setText(f"Added {len(feats)} features")

    def name_from_references(self):
        if not self.orfs:
            QMessageBox.information(self, "ORF Finder", "Run the search first."); return
        path, _ = QFileDialog.getOpenFileName(self, "Reference peptides (FASTA)", "",
                                              "FASTA (*.fasta *.fa *.faa);;All files (*)")
        if not path:
            return
        try:
            refs = GX.load_reference_peptides(path)
        except Exception as e:
            QMessageBox.critical(self, "Reference peptides", str(e)); return
        if not refs:
            QMessageBox.information(self, "Reference peptides", "No sequences found in that file."); return
        matches = GX.name_by_similarity(self.orfs, refs, min_identity=0.5)
        for i, (name, ident, note) in matches.items():
            self.orfs[i].name = name
            self.orfs[i].extra["note"] = note
        self._refill_names()
        self.orfsFound.emit(self.orfs, self.table.currentData())
        self.status.setText(f"Named {len(matches)} of {len(self.orfs)} ORFs from {len(refs)} reference peptide(s) "
                            f"(\u2265 50% identity); the rest keep descriptive names")

    def _refill_names(self):
        self._loading = True
        for i in range(self.table_w.rowCount()):
            it = self.table_w.item(i, 1)
            k = it.data(Qt.UserRole + 1)
            if k is not None and 0 <= k < len(self.orfs):
                it.setText(self.orfs[k].name)
        self._loading = False

    def export_genbank(self):
        if not self.orfs:
            QMessageBox.information(self, "ORF Finder", "Run the search first."); return
        rows = sorted({o.row for o in self.orfs})
        if len(rows) > 1:
            QMessageBox.information(self, "GenBank export",
                                    "Select a single sequence (Sequence: …) before exporting GenBank."); return
        row = rows[0]
        dlg = GenBankExportOptions(self, table=self.table.currentData(),
                                   seq_name=self.model.rows[row].name,
                                   source_path=getattr(self.parent(), "model", None) and getattr(self.parent().model, "path", None))
        if not dlg.exec():
            return
        opts = dlg.values()
        seq = self.model.rows[row].seq
        source = None
        if opts["use_source"] and opts["source_path"]:
            try:
                source = GX.load_source_record(opts["source_path"])
            except Exception as e:
                QMessageBox.critical(self, "GenBank", f"Could not read source record:\n{e}"); return
            if len("".join(c for c in seq if c not in "-.~")) != len(source.seq):
                r = QMessageBox.question(self, "GenBank",
                                         "The sequence length differs from the source record. Continue and keep the "
                                         "source record's sequence and features?")
                if r != QMessageBox.Yes:
                    return
        anns = []
        for o in self.orfs:
            if o.row != row:
                continue
            anns.append(GX.ORFAnnotation(o, name=o.name, product=o.name, table=opts["table"],
                                         cytoplasmic=opts["cytoplasmic"], feature_type=opts["feature_type"],
                                         note=o.extra.get("note", "")))
        rec = GX.build_record(seq, anns, source=source, record_id=opts["record_id"],
                              organism=opts["organism"], topology=opts["topology"],
                              default_table=opts["default_table"], keep_source_features=True)
        path, _ = QFileDialog.getSaveFileName(self, "Export annotated GenBank",
                                              f"{opts['record_id'] or self.model.rows[row].name}_annotated.gb",
                                              "GenBank (*.gb *.gbk)")
        if not path:
            return
        try:
            GX.write_genbank(rec, path)
        except Exception as e:
            QMessageBox.critical(self, "GenBank", str(e)); return
        n1 = sum(1 for f in rec.features if f.qualifiers.get("transl_table") == [str(opts["table"])])
        self.status.setText(f"Wrote {path}: {len(rec.features)} features "
                            f"({len(anns)} added with transl_table={opts['table']}, total {n1} using that code)")

    def export(self, kind):
        if not self.orfs:
            QMessageBox.information(self, "ORF Finder", "Run the search first.")
            return
        by_row = {}
        for o in self.orfs:
            by_row.setdefault(o.row, []).append(o)
        chunks = []
        for r, lst in by_row.items():
            sid = self.model.rows[r].name
            if kind == "gff":
                chunks.append(OF.orfs_to_gff(lst, sid))
            elif kind == "tbl":
                chunks.append(OF.orfs_to_genbank_table(lst, sid))
            elif kind == "aa":
                chunks.append(OF.orfs_to_fasta(lst, sid, protein=True))
            else:
                chunks.append(OF.orfs_to_fasta(lst, sid, protein=False))
        filt = {"gff": "GFF3 (*.gff3 *.gff)", "tbl": "Feature table (*.tbl)", "aa": "FASTA (*.faa *.fasta)", "nt": "FASTA (*.fna *.fasta)"}[kind]
        save_text(self, "".join(chunks), "Export ORFs", filt)


class GenBankExportOptions(QDialog):
    """How to write the ORFs into a GenBank record."""

    def __init__(self, parent=None, table=1, seq_name="", source_path=None):
        super().__init__(parent)
        self.setWindowTitle("Export annotated GenBank")
        self.resize(620, 340)
        lay = QVBoxLayout(self)
        info = QLabel("Mitochondrial-derived peptides (humanin, MOTS-c, SHLPs …) are encoded in mtDNA but "
                      "translated by <b>cytoplasmic</b> ribosomes, so they use the standard genetic code, "
                      "while the surrounding mitochondrial genes use the mitochondrial code. Each feature is "
                      "written with its own <tt>/transl_table</tt>, so both classes display correctly in any "
                      "GenBank viewer.")
        info.setWordWrap(True)
        lay.addWidget(info)
        form = QFormLayout()
        self.table = codon_table_combo(table)
        self.default_table = codon_table_combo(2)
        self.cytoplasmic = QCheckBox("Add note: translated by cytoplasmic ribosomes")
        self.cytoplasmic.setChecked(table == 1)
        self.ftype = QComboBox(); self.ftype.addItems(["CDS", "misc_feature", "mat_peptide", "gene"])
        self.use_source = QCheckBox("Merge into an existing GenBank record (keeps its genes)")
        self.use_source.setChecked(bool(source_path and str(source_path).lower().endswith((".gb", ".gbk", ".genbank"))))
        self.source_path = QLineEdit(source_path or "")
        browse = QPushButton("…"); browse.clicked.connect(self._browse)
        srow = QHBoxLayout(); srow.addWidget(self.source_path); srow.addWidget(browse)
        sw = QWidget(); sw.setLayout(srow)
        self.record_id = QLineEdit(seq_name)
        self.organism = QLineEdit("")
        self.topology = QComboBox(); self.topology.addItems(["circular", "linear"])
        form.addRow("Genetic code for these ORFs", self.table)
        form.addRow("", self.cytoplasmic)
        form.addRow("Feature type", self.ftype)
        form.addRow("Genome's own code (new records)", self.default_table)
        form.addRow("", self.use_source)
        form.addRow("Source record", sw)
        form.addRow("Record ID / LOCUS", self.record_id)
        form.addRow("Organism", self.organism)
        form.addRow("Topology", self.topology)
        lay.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Export…")
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def _browse(self):
        p, _ = QFileDialog.getOpenFileName(self, "Source GenBank record", "", "GenBank (*.gb *.gbk *.genbank)")
        if p:
            self.source_path.setText(p); self.use_source.setChecked(True)

    def values(self):
        return dict(table=self.table.currentData(), default_table=self.default_table.currentData(),
                    cytoplasmic=self.cytoplasmic.isChecked(), feature_type=self.ftype.currentText(),
                    use_source=self.use_source.isChecked(), source_path=self.source_path.text().strip(),
                    record_id=self.record_id.text().strip(), organism=self.organism.text().strip(),
                    topology=self.topology.currentText())
