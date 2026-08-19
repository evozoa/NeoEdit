from __future__ import annotations

from PySide6.QtGui import QFont
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QDialog, QFormLayout, QComboBox, QSpinBox, QCheckBox, QPushButton, QHBoxLayout,
                               QVBoxLayout, QLabel, QGroupBox, QTableWidgetItem, QMessageBox, QFileDialog,
                               QSplitter, QPlainTextEdit, QWidget)

from ...analysis import orf_finder as OF
from ...model.alignment import Feature
from .common import codon_table_combo, make_table, NumItem, save_text


class ORFFinderDialog(QDialog):
    """MitoFinder-style ORF finder with alternate genetic codes."""
    orfSelected = Signal(int, int, int)        # row, start(gapped), end(gapped)
    featuresReady = Signal(list)               # list[Feature]

    def __init__(self, model, rows: list[int], parent=None, default_table=1):
        super().__init__(parent)
        self.setWindowTitle("ORF Finder")
        self.resize(900, 600)
        self.model = model
        self.rows = rows
        self.orfs: list[OF.ORF] = []
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
        for b in (self.run_btn, self.annot_btn, self.exp_gff, self.exp_tbl, self.exp_aa, self.exp_nt):
            btns.addWidget(b)
        btns.addStretch()
        lay.addLayout(btns)

        split = QSplitter(Qt.Vertical)
        self.table_w = make_table(["Seq", "Start", "End", "Strand", "Frame", "Length (nt)", "Length (aa)",
                                   "Start codon", "Stop codon", "Partial", "Protein"])
        self.table_w.itemSelectionChanged.connect(self._on_select)
        split.addWidget(self.table_w)
        self.detail = QPlainTextEdit(); self.detail.setReadOnly(True)
        f = self.detail.font(); f.setFamily("DejaVu Sans Mono"); f.setStyleHint(QFont.Monospace); self.detail.setFont(f)
        split.addWidget(self.detail)
        split.setSizes([400, 150])
        lay.addWidget(split, 1)
        self.status = QLabel("")
        lay.addWidget(self.status)

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
        for i, o in enumerate(self.orfs):
            self.table_w.insertRow(i)
            partial = ("5'" if o.partial5 else "") + ("3'" if o.partial3 else "")
            vals = [self.model.rows[o.row].name, NumItem(o.start + 1), NumItem(o.end), "+" if o.strand > 0 else "-",
                    NumItem(o.frame), NumItem(o.length_nt), NumItem(o.length_aa), o.start_codon, o.stop_codon or "(none)",
                    partial or "-", o.aa[:60] + ("…" if len(o.aa) > 60 else "")]
            for j, v in enumerate(vals):
                item = v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v))
                item.setData(Qt.UserRole + 1, i)
                self.table_w.setItem(i, j, item)
        self.table_w.setSortingEnabled(True)
        self.table_w.sortItems(6, Qt.DescendingOrder)
        self.status.setText(f"{len(self.orfs)} ORF(s) found")

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
