"""Restriction map: enzyme sites in one sequence or summarised across an alignment."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QColor, QBrush
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QComboBox, QSpinBox,
                               QCheckBox, QPushButton, QLabel, QLineEdit, QSplitter, QPlainTextEdit,
                               QTableWidgetItem, QMessageBox, QTabWidget, QListWidget, QListWidgetItem)

from ...analysis import restriction as R
from .common import make_table, NumItem, save_text

GREEN = QColor("#d7f5d7"); RED = QColor("#ffdcdc")


class RestrictionDialog(QDialog):
    featuresReady = Signal(list)
    siteSelected = Signal(int, int, int)      # row, start, end (alignment columns)

    def __init__(self, model, view, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Restriction sites")
        self.resize(1000, 700)
        self.setModal(False)
        self.model = model
        self.view = view
        self.hits: list[R.EnzymeHit] = []
        self.summaries: list[R.AlignmentSummary] = []

        lay = QVBoxLayout(self)
        top = QHBoxLayout()

        gsel = QGroupBox("Sequences")
        sf = QFormLayout(gsel)
        self.scope = QComboBox()
        self.scope.addItem("Current sequence", "cur")
        self.scope.addItem("Selected sequences", "sel")
        self.scope.addItem("All sequences (compare)", "all")
        self.seq_combo = QComboBox()
        for i, r in enumerate(model.rows):
            self.seq_combo.addItem(f"{i + 1}. {r.name}", i)
        self.seq_combo.setCurrentIndex(min(view.cur_row, model.nrows - 1))
        self.circular = QCheckBox("Circular molecule (plasmid / mitogenome)")
        self.circular.setChecked(bool(getattr(model, "circular", False)))
        sf.addRow("Search in", self.scope)
        sf.addRow("Sequence", self.seq_combo)
        sf.addRow("", self.circular)
        top.addWidget(gsel, 1)

        genz = QGroupBox("Enzymes")
        ef = QFormLayout(genz)
        self.commercial = QCheckBox("Commercially available only"); self.commercial.setChecked(True)
        self.min_site = QSpinBox(); self.min_site.setRange(3, 12); self.min_site.setValue(6)
        self.max_site = QSpinBox(); self.max_site.setRange(3, 20); self.max_site.setValue(8)
        self.ends = QComboBox(); self.ends.addItem("any ends", None); self.ends.addItem("blunt only", True); self.ends.addItem("sticky only", False)
        self.supplier = QComboBox(); self.supplier.addItem("any supplier", None)
        for sname in R.all_suppliers():
            self.supplier.addItem(sname, sname)
        self.names = QLineEdit(); self.names.setPlaceholderText("or list enzymes: EcoRI, BamHI, HindIII…")
        ef.addRow("", self.commercial)
        ef.addRow("Site length", self.min_site)
        ef.addRow("… to", self.max_site)
        ef.addRow("Ends", self.ends)
        ef.addRow("Supplier", self.supplier)
        ef.addRow("Specific enzymes", self.names)
        top.addWidget(genz, 1)

        gfil = QGroupBox("Report")
        ff = QFormLayout(gfil)
        self.min_cuts = QSpinBox(); self.min_cuts.setRange(0, 100); self.min_cuts.setValue(1)
        self.max_cuts = QSpinBox(); self.max_cuts.setRange(0, 1000); self.max_cuts.setValue(0); self.max_cuts.setSpecialValueText("no limit")
        self.unique_only = QCheckBox("Unique cutters only (exactly 1 site)")
        self.diag_only = QCheckBox("Only enzymes that distinguish sequences")
        self.noncut = QCheckBox("Include non-cutters (0 sites)")
        ff.addRow("Min cuts", self.min_cuts)
        ff.addRow("Max cuts", self.max_cuts)
        ff.addRow("", self.unique_only)
        ff.addRow("", self.diag_only)
        ff.addRow("", self.noncut)
        top.addWidget(gfil, 1)
        lay.addLayout(top)

        row = QHBoxLayout()
        self.run_btn = QPushButton("Find sites"); self.run_btn.setDefault(True); self.run_btn.clicked.connect(self.run)
        self.feat_btn = QPushButton("Add sites as features"); self.feat_btn.clicked.connect(self.add_features)
        self.csv_btn = QPushButton("Export CSV"); self.csv_btn.clicked.connect(self.export_csv)
        self.digest_btn = QPushButton("Digest with selected"); self.digest_btn.clicked.connect(self.digest)
        for b in (self.run_btn, self.feat_btn, self.csv_btn, self.digest_btn):
            row.addWidget(b)
        row.addStretch()
        lay.addLayout(row)

        split = QSplitter(Qt.Vertical)
        self.table = make_table(["Enzyme", "Site", "Cut", "Ends", "Cuts", "Positions", "Suppliers"])
        self.table.itemSelectionChanged.connect(self._show)
        split.addWidget(self.table)
        self.detail = QPlainTextEdit(); self.detail.setReadOnly(True)
        f = self.detail.font(); f.setFamily("Courier New"); f.setStyleHint(QFont.Monospace); self.detail.setFont(f)
        split.addWidget(self.detail)
        split.setSizes([420, 200])
        lay.addWidget(split, 1)
        self.status = QLabel("")
        lay.addWidget(self.status)
        self.scope.currentIndexChanged.connect(lambda: self.seq_combo.setEnabled(self.scope.currentData() == "cur"))

    # ------------------------------------------------------------- helpers
    def _pool(self):
        names = [n for n in self.names.text().replace(";", ",").split(",") if n.strip()]
        return R.enzyme_pool(commercial_only=self.commercial.isChecked(),
                             suppliers=[self.supplier.currentData()] if self.supplier.currentData() else None,
                             min_site=self.min_site.value(), max_site=self.max_site.value(),
                             blunt=self.ends.currentData(), names=names or None)

    def _rows(self):
        sc = self.scope.currentData()
        if sc == "cur":
            return [self.seq_combo.currentData()]
        if sc == "sel":
            return self.view.target_rows()
        return list(range(self.model.nrows))

    # ------------------------------------------------------------- run
    def run(self):
        pool = self._pool()
        rows = self._rows()
        if not rows:
            return
        linear = not self.circular.isChecked()
        min_cuts = 0 if self.noncut.isChecked() else self.min_cuts.value()
        max_cuts = None if self.max_cuts.value() == 0 else self.max_cuts.value()
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        self.hits = []; self.summaries = []
        seqs = [r.seq for r in self.model.rows]
        if len(rows) > 1:
            self.summaries = R.search_alignment(seqs, pool, linear, rows)
            if self.diag_only.isChecked():
                self.summaries = [s for s in self.summaries if s.is_diagnostic(rows)]
            if self.unique_only.isChecked():
                self.summaries = [s for s in self.summaries if all(len(s.cuts_per_row.get(r, [])) <= 1 for r in rows)]
            self._fill_summaries(rows)
        else:
            self.hits = R.search_sequence(seqs[rows[0]], pool, linear, rows[0], min_cuts, max_cuts)
            if self.unique_only.isChecked():
                self.hits = [h for h in self.hits if h.n_cuts == 1]
            self._fill_hits()

    def _fill_hits(self):
        self.table.setColumnCount(7)
        self.table.setHorizontalHeaderLabels(["Enzyme", "Site", "Cut", "Ends", "Cuts", "Positions", "Suppliers"])
        for i, h in enumerate(self.hits):
            self.table.insertRow(i)
            pos = ", ".join(str(p) for p in h.positions[:12]) + ("…" if len(h.positions) > 12 else "")
            vals = [h.enzyme, h.site, h.elucidate, h.overhang, NumItem(h.n_cuts), pos or "—", ", ".join(h.suppliers[:3])]
            for j, v in enumerate(vals):
                it = v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v))
                it.setData(Qt.UserRole + 1, i)
                if j == 4:
                    it.setBackground(QBrush(GREEN if h.n_cuts == 1 else (RED if h.n_cuts == 0 else QColor("#ffffff"))))
                self.table.setItem(i, j, it)
        self.table.setSortingEnabled(True)
        self.status.setText(f"{len(self.hits)} enzyme(s)")

    def _fill_summaries(self, rows):
        names = [self.model.rows[r].name for r in rows]
        self.table.setColumnCount(4 + len(rows))
        self.table.setHorizontalHeaderLabels(["Enzyme", "Site", "Ends", "Total cuts"] + names)
        for i, s in enumerate(self.summaries):
            self.table.insertRow(i)
            total = sum(len(s.cuts_per_row.get(r, [])) for r in rows)
            vals = [s.enzyme, s.site, s.overhang, NumItem(total)]
            for j, v in enumerate(vals):
                it = v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v))
                it.setData(Qt.UserRole + 1, i)
                self.table.setItem(i, j, it)
            for k, r in enumerate(rows):
                n = len(s.cuts_per_row.get(r, []))
                it = NumItem(n)
                it.setBackground(QBrush(GREEN if n else RED))
                it.setData(Qt.UserRole + 1, i)
                self.table.setItem(i, 4 + k, it)
        self.table.setSortingEnabled(True)
        diag = sum(1 for s in self.summaries if s.is_diagnostic(rows))
        self.status.setText(f"{len(self.summaries)} enzyme(s) across {len(rows)} sequences — {diag} distinguish them")

    def _current(self):
        items = self.table.selectedItems()
        if not items:
            return None
        i = items[0].data(Qt.UserRole + 1)
        return (self.summaries or self.hits)[i]

    def _show(self):
        obj = self._current()
        if obj is None:
            return
        if isinstance(obj, R.EnzymeHit):
            gm = R.ungapped_to_columns(self.model.rows[obj.row].seq)
            lines = [f"{obj.enzyme}   {obj.elucidate}   {obj.overhang} ends   {obj.n_cuts} site(s)",
                     f"Suppliers: {', '.join(obj.suppliers) or '—'}", ""]
            seq = "".join(c for c in self.model.rows[obj.row].seq if c not in "-.~").upper()
            for cut in obj.positions:
                st = obj.site_start(cut)
                ctx = seq[max(0, st - 10):st] + "[" + seq[st:st + obj.size] + "]" + seq[st + obj.size:st + obj.size + 10]
                lines.append(f"  site {st + 1:>8,}   cut after {cut - 1:,}   {ctx}")
            self.detail.setPlainText("\n".join(lines))
            if obj.positions and gm:
                st = obj.site_start(obj.positions[0])
                if st < len(gm):
                    self.siteSelected.emit(obj.row, gm[st], gm[min(len(gm) - 1, st + obj.size - 1)] + 1)
        else:
            rows = self._rows()
            lines = [f"{obj.enzyme}   {obj.site}   {obj.overhang} ends", ""]
            for r in rows:
                pos = obj.cuts_per_row.get(r, [])
                lines.append(f"  {self.model.rows[r].name[:30].ljust(30)} {len(pos):>3} site(s)  " +
                             (", ".join(str(p) for p in pos[:10]) if pos else "— none —"))
            self.detail.setPlainText("\n".join(lines))

    # ------------------------------------------------------------- outputs
    def add_features(self):
        hits = self.hits
        if not hits and self.summaries:
            obj = self._current()
            if obj is None:
                return
            hits = []
            for r, pos in obj.cuts_per_row.items():
                if pos:
                    hits.append(R.EnzymeHit(obj.enzyme, obj.site, "", obj.size, obj.overhang, [], pos, 1, r))
        if not hits:
            return
        feats = R.hits_to_features(hits, {i: r.seq for i, r in enumerate(self.model.rows)})
        self.featuresReady.emit(feats)
        self.status.setText(f"Added {len(feats)} site features")

    def export_csv(self):
        if self.summaries:
            rows = self._rows()
            head = "enzyme,site,ends," + ",".join(self.model.rows[r].name for r in rows)
            lines = [head]
            for s in self.summaries:
                lines.append(f"{s.enzyme},{s.site},{s.overhang}," +
                             ",".join(str(len(s.cuts_per_row.get(r, []))) for r in rows))
        elif self.hits:
            lines = ["enzyme,site,cut_pattern,ends,n_cuts,positions,suppliers"]
            for h in self.hits:
                lines.append(f"{h.enzyme},{h.site},{h.elucidate},{h.overhang},{h.n_cuts},"
                             f"\"{' '.join(str(p) for p in h.positions)}\",\"{'; '.join(h.suppliers)}\"")
        else:
            return
        save_text(self, "\n".join(lines) + "\n", "Export restriction map", "CSV (*.csv)")

    def digest(self):
        obj = self._current()
        if obj is None:
            return
        row = obj.row if isinstance(obj, R.EnzymeHit) else self._rows()[0]
        frags = R.digest_fragments(self.model.rows[row].seq, [obj.enzyme], linear=not self.circular.isChecked())
        lines = [f"Digest of {self.model.rows[row].name} with {obj.enzyme} "
                 f"({'circular' if self.circular.isChecked() else 'linear'}): {len(frags)} fragment(s)", ""]
        for i, (a, b, L) in enumerate(sorted(frags, key=lambda f: -f[2]), 1):
            lines.append(f"  {i:>3}. {L:>10,} bp   {a + 1:,}–{b:,}")
        self.detail.setPlainText("\n".join(lines))
