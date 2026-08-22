"""Alignment-aware primer design: conserved (universal) or discriminating (eDNA)."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QColor, QBrush
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout, QGroupBox, QLabel,
                               QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox, QPushButton, QSplitter,
                               QPlainTextEdit, QTableWidgetItem, QMessageBox, QTabWidget, QWidget, QLineEdit)

from ...analysis import primer_design as PD
from ...analysis import primers as P
from ...model.alignment import Feature
from .common import make_table, NumItem, save_text

GREEN = QColor("#d7f5d7"); RED = QColor("#ffdcdc"); GREY = QColor("#eeeeee")


class DesignDialog(QDialog):
    featuresReady = Signal(list)
    regionSelected = Signal(int, int)

    def __init__(self, model, view, parent=None, circular: bool = False):
        super().__init__(parent)
        self.setWindowTitle("Design primers across an alignment")
        self._circular_default = bool(circular)
        self.resize(1100, 720)
        self.setModal(False)
        self.model = model
        self.view = view
        self.evals: list[PD.PairEvaluation] = []
        self.result = None

        lay = QVBoxLayout(self)

        # ---------------- parameters
        top = QHBoxLayout()
        gsel = QGroupBox("Sequences")
        gf = QFormLayout(gsel)
        groups = model.groups()
        self.mode = QComboBox()
        self.mode.addItem("Conserved / universal primers (amplify all)", False)
        self.mode.addItem("Discriminating primers (amplify include, exclude others)", True)
        self.mode.currentIndexChanged.connect(self._mode_changed)
        self.include = QComboBox(); self.exclude = QComboBox(); self.template = QComboBox()
        self.include.addItem("All sequences", "__all__")
        self.exclude.addItem("(none)", "__none__")
        for g in groups:
            n = len(model.group_rows(g))
            self.include.addItem(f"group: {g} ({n})", g)
            self.exclude.addItem(f"group: {g} ({n})", g)
        self.include.addItem("Selected rows", "__sel__")
        self.exclude.addItem("Selected rows", "__sel__")
        for i, r in enumerate(model.rows):
            self.template.addItem(f"{i + 1}. {r.name}", i)
        if groups:
            i = self.include.findData(groups[0])
            if i >= 0:
                self.include.setCurrentIndex(i)
            if len(groups) > 1:
                j = self.exclude.findData(groups[1])
                if j >= 0:
                    self.exclude.setCurrentIndex(j)
        gf.addRow("Design", self.mode)
        gf.addRow("Include (must amplify)", self.include)
        gf.addRow("Exclude (must NOT amplify)", self.exclude)
        gf.addRow("Template (Primer3 designs on)", self.template)
        self.region = QComboBox()
        self.region.addItem("Whole alignment", None)
        self.region.addItem("Current selection", "sel")
        gf.addRow("Restrict to", self.region)
        top.addWidget(gsel, 1)

        gcons = QGroupBox("Conservation & stringency")
        cf = QFormLayout(gcons)
        self.min_cons = QDoubleSpinBox(); self.min_cons.setRange(0.5, 1.0); self.min_cons.setSingleStep(0.01); self.min_cons.setValue(1.0)
        self.min_run = QSpinBox(); self.min_run.setRange(10, 60); self.min_run.setValue(20)
        self.max_mm = QSpinBox(); self.max_mm.setRange(0, 10); self.max_mm.setValue(2)
        self.max_3p = QSpinBox(); self.max_3p.setRange(0, 6); self.max_3p.setValue(0)
        self.win_3p = QSpinBox(); self.win_3p.setRange(1, 12); self.win_3p.setValue(5)
        self.require_last = QCheckBox("3'-terminal mismatch blocks amplification"); self.require_last.setChecked(True)
        self.degenerate = QCheckBox("Also propose IUPAC-degenerate primers")
        self.max_deg = QSpinBox(); self.max_deg.setRange(1, 4096); self.max_deg.setValue(64)
        cf.addRow("Min conservation in include set", self.min_cons)
        cf.addRow("Min conserved run (bp)", self.min_run)
        cf.addRow("Max mismatches / primer", self.max_mm)
        cf.addRow("Max mismatches in 3' window", self.max_3p)
        cf.addRow("3' window size (bp)", self.win_3p)
        cf.addRow("", self.require_last)
        cf.addRow("", self.degenerate)
        cf.addRow("Max degeneracy", self.max_deg)
        top.addWidget(gcons, 1)

        gp3 = QGroupBox("Primer3")
        pf = QFormLayout(gp3)
        self.prod_min = QSpinBox(); self.prod_min.setRange(50, 100000); self.prod_min.setValue(120)
        self.prod_max = QSpinBox(); self.prod_max.setRange(50, 100000); self.prod_max.setValue(400)
        self.opt_tm = QDoubleSpinBox(); self.opt_tm.setRange(40, 80); self.opt_tm.setValue(58)
        self.tm_span = QDoubleSpinBox(); self.tm_span.setRange(1, 15); self.tm_span.setValue(4)
        self.num = QSpinBox(); self.num.setRange(1, 200); self.num.setValue(30)
        pf.addRow("Product size min", self.prod_min)
        pf.addRow("Product size max", self.prod_max)
        self.circ_cb = QCheckBox("Circular molecule: products may span the origin")
        self.circ_cb.setChecked(self._circular_default)
        pf.addRow("", self.circ_cb)
        pf.addRow("Optimal Tm (°C)", self.opt_tm)
        pf.addRow("Tm ± (°C)", self.tm_span)
        pf.addRow("Candidates to evaluate", self.num)
        top.addWidget(gp3, 1)
        lay.addLayout(top)

        row = QHBoxLayout()
        self.run_btn = QPushButton("Design and rank"); self.run_btn.setDefault(True); self.run_btn.clicked.connect(self.run)
        self.feat_btn = QPushButton("Add selected pair as features"); self.feat_btn.clicked.connect(self.add_features)
        self.csv_btn = QPushButton("Export candidates (CSV)"); self.csv_btn.clicked.connect(self.export_csv)
        self.pcr_btn = QPushButton("Export in-silico PCR (CSV)"); self.pcr_btn.clicked.connect(self.export_pcr)
        for b in (self.run_btn, self.feat_btn, self.csv_btn, self.pcr_btn):
            row.addWidget(b)
        row.addStretch()
        lay.addLayout(row)

        split = QSplitter(Qt.Vertical)
        self.cand = make_table(["#", "Score", "Forward (5'→3')", "Reverse (5'→3')", "Product",
                                "Include amplified", "Exclude amplified", "Fwd Tm", "Rev Tm", "Position"])
        self.cand.itemSelectionChanged.connect(self._show_pair)
        split.addWidget(self.cand)
        tabs = QTabWidget()
        self.pcr = make_table(["Sequence", "Group", "Fwd mm", "Fwd 3'", "Rev mm", "Rev 3'", "Amplifies"])
        tabs.addTab(self.pcr, "In-silico PCR")
        self.detail = QPlainTextEdit(); self.detail.setReadOnly(True)
        f = self.detail.font(); f.setFamily("Courier New"); f.setStyleHint(QFont.Monospace); self.detail.setFont(f)
        tabs.addTab(self.detail, "Alignment of primer sites")
        split.addWidget(tabs)
        split.setSizes([300, 320])
        lay.addWidget(split, 1)
        self.status = QLabel("")
        lay.addWidget(self.status)
        self._mode_changed()

    # ------------------------------------------------------------- helpers
    def _mode_changed(self):
        disc = self.mode.currentData()
        self.exclude.setEnabled(bool(disc))

    def _rows_for(self, combo) -> list[int]:
        key = combo.currentData()
        if key == "__all__":
            return list(range(self.model.nrows))
        if key == "__none__":
            return []
        if key == "__sel__":
            return self.view.target_rows()
        return self.model.group_rows(key)

    def _kw(self):
        return dict(max_mm=self.max_mm.value(), max_3p=self.max_3p.value(), require_last=self.require_last.isChecked())

    # ------------------------------------------------------------- run
    def run(self):
        inc = self._rows_for(self.include)
        exc = self._rows_for(self.exclude) if self.mode.currentData() else []
        exc = [r for r in exc if r not in inc]
        if not inc:
            QMessageBox.information(self, "Design", "The include set is empty."); return
        tmpl = self.template.currentData()
        region = None
        if self.region.currentData() == "sel":
            s = self.view.selection()
            if s:
                region = (s[2], s[3] + 1)
        rows = [r.seq for r in self.model.rows]
        names = [r.name for r in self.model.rows]
        try:
            res = PD.design_on_alignment(
                rows, names, template_row=tmpl, include_rows=inc, exclude_rows=exc, region=region,
                discriminating=bool(self.mode.currentData()),
                min_conservation=self.min_cons.value(), min_conserved_run=self.min_run.value(),
                product_range=((self.prod_min.value(), self.prod_max.value()),), num_return=self.num.value(),
                three_prime_window=self.win_3p.value(), degenerate=self.degenerate.isChecked(),
                max_degeneracy=self.max_deg.value(), circular=self.circ_cb.isChecked(),
                primer3_kwargs=dict(opt_tm=self.opt_tm.value(), min_tm=self.opt_tm.value() - self.tm_span.value(),
                                    max_tm=self.opt_tm.value() + self.tm_span.value()),
                **self._kw())
            self.result = res
            self.evals = list(res)
        except Exception as e:
            self.evals = []
            self.cand.setRowCount(0); self.pcr.setRowCount(0)
            self.status.setText(str(e)[:500]); return
        self._fill_candidates()

    def _fill_candidates(self):
        kw = self._kw()
        disc = bool(self.mode.currentData())
        self.cand.setSortingEnabled(False)
        self.cand.setRowCount(0)
        for i, ev in enumerate(self.evals):
            st = ev.stats(**kw)
            self.cand.insertRow(i)
            pos = f"{ev.left_cols[0] + 1:,}–{ev.right_cols[1]:,}" + (" ⟳" if ev.pair.crosses_origin else "")
            vals = [NumItem(i + 1), NumItem(ev.score(disc, **kw), "{:.1f}"),
                    ev.left_seq_deg or ev.pair.left.seq, ev.right_seq_deg or ev.pair.right.seq,
                    NumItem(ev.pair.product_size),
                    f"{st['include_hit']}/{st['include_total']}",
                    f"{st['exclude_hit']}/{st['exclude_total']}" if st["exclude_total"] else "—",
                    NumItem(ev.pair.left.tm, "{:.1f}"), NumItem(ev.pair.right.tm, "{:.1f}"), pos]
            for j, v in enumerate(vals):
                it = v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v))
                it.setData(Qt.UserRole + 1, i)
                if j == 5:
                    it.setBackground(QBrush(GREEN if st["include_hit"] == st["include_total"] else RED))
                if j == 6 and st["exclude_total"]:
                    it.setBackground(QBrush(GREEN if st["exclude_hit"] == 0 else RED))
                self.cand.setItem(i, j, it)
        self.cand.setSortingEnabled(True)
        self.cand.sortItems(1, Qt.DescendingOrder)
        msg = f"{len(self.evals)} candidate pairs evaluated against {self.model.nrows} sequences"
        if self.result is not None:
            msg += f"  |  conserved windows: {self.result.n_conserved_runs} (longest {self.result.longest_run} bp)"
            if not self.result.mask_applied:
                msg = "⚠ " + self.result.notes + "\n" + msg
        self.status.setText(msg)
        self.status.setWordWrap(True)
        self.status.setStyleSheet("color:#a00" if (self.result is not None and not self.result.mask_applied) else "")
        if self.evals:
            self.cand.selectRow(0)

    def _current(self):
        items = self.cand.selectedItems()
        if not items:
            return None
        return self.evals[items[0].data(Qt.UserRole + 1)]

    def _show_pair(self):
        ev = self._current()
        if ev is None:
            return
        kw = self._kw()
        tbl = PD.insilico_table(ev, **kw)
        self.pcr.setSortingEnabled(False)
        self.pcr.setRowCount(0)
        for i, r in enumerate(tbl):
            self.pcr.insertRow(i)
            vals = [r["name"], r["group"], NumItem(r["fwd_mm"] if r["fwd_cov"] else -1),
                    NumItem(r["fwd_3p"] if r["fwd_cov"] else -1), NumItem(r["rev_mm"] if r["rev_cov"] else -1),
                    NumItem(r["rev_3p"] if r["rev_cov"] else -1), "yes" if r["amplifies"] else "no"]
            for j, v in enumerate(vals):
                it = v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v))
                if isinstance(v, NumItem) and v.data(Qt.UserRole) == -1:
                    it.setText("n/c")
                if j == 6:
                    it.setBackground(QBrush(GREEN if r["amplifies"] else RED))
                elif not (r["fwd_cov"] and r["rev_cov"]):
                    it.setBackground(QBrush(GREY))
                self.pcr.setItem(i, j, it)
        self.pcr.setSortingEnabled(True)
        # alignment-of-sites view
        lines = [f"Pair: F {ev.pair.left.seq}  R {ev.pair.right.seq}   product {ev.pair.product_size} bp",
                 f"Columns: forward {ev.cols_text('left')}, reverse {ev.cols_text('right')}"
                 + ("   — product across the origin" if ev.pair.crosses_origin else ""), ""]
        for label, primer, hits in (("FORWARD", ev.left_seq_deg or ev.pair.left.seq, ev.left_hits),
                                    ("REVERSE", ev.right_seq_deg or ev.pair.right.seq, ev.right_hits)):
            lines.append(f"{label}  {primer}   (5'->3'; '.' = match, lower-case = mismatch, 3' end at right)")
            for h in hits:
                if not h.covered:
                    lines.append(f"   {h.name[:26].ljust(26)} {'(not covered)'}")
                    continue
                marks = "".join("." if PD.matches(p, t) else (t.lower() if t else "-")
                                for p, t in zip(primer.upper(), h.site.ljust(len(primer))))
                flag = "  <-- 3' mismatch" if h.last_base_mm else ("  <- 3' region" if h.mm_3prime else "")
                grp = "I" if h.row in ev.include_rows else ("X" if h.row in ev.exclude_rows else " ")
                lines.append(f" {grp} {h.name[:26].ljust(26)} {marks}  {h.mismatches} mm{flag}")
            lines.append("")
        self.detail.setPlainText("\n".join(lines))
        self.regionSelected.emit(ev.left_cols[0], self.model.width if ev.pair.crosses_origin else ev.right_cols[1])

    # ------------------------------------------------------------- outputs
    def add_features(self):
        ev = self._current()
        if ev is None:
            return
        feats = []
        for pieces, strand, nm, seq in ((ev.left_pieces or [ev.left_cols], 1, "F", ev.left_seq_deg or ev.pair.left.seq),
                                        (ev.right_pieces or [ev.right_cols], -1, "R", ev.right_seq_deg or ev.pair.right.seq)):
            for c0, c1 in pieces:
                for r in range(self.model.nrows):
                    feats.append(Feature(r, c0, c1, strand, "primer", f"{nm} {seq}",
                                         "#10b981" if strand > 0 else "#ef4444", data={"seq": seq}))
        self.featuresReady.emit(feats)
        self.status.setText(f"Added primer features to {self.model.nrows} sequences")

    def export_csv(self):
        if not self.evals:
            return
        kw = self._kw(); disc = bool(self.mode.currentData())
        lines = ["rank,score,forward,reverse,product_size,fwd_tm,rev_tm,include_hit,include_total,exclude_hit,exclude_total,fwd_cols,rev_cols"]
        for i, ev in enumerate(self.evals, 1):
            st = ev.stats(**kw)
            lines.append(f"{i},{ev.score(disc, **kw):.2f},{ev.left_seq_deg or ev.pair.left.seq},{ev.right_seq_deg or ev.pair.right.seq},"
                         f"{ev.pair.product_size},{ev.pair.left.tm:.1f},{ev.pair.right.tm:.1f},{st['include_hit']},{st['include_total']},"
                         f"{st['exclude_hit']},{st['exclude_total']},{ev.left_cols[0] + 1}-{ev.left_cols[1]},{ev.right_cols[0] + 1}-{ev.right_cols[1]}")
        save_text(self, "\n".join(lines) + "\n", "Export candidates", "CSV (*.csv)")

    def export_pcr(self):
        ev = self._current()
        if ev is None:
            return
        kw = self._kw()
        lines = ["sequence,group,fwd_mismatches,fwd_3prime,fwd_terminal,rev_mismatches,rev_3prime,rev_terminal,amplifies"]
        for r in PD.insilico_table(ev, **kw):
            lines.append(f"{r['name']},{r['group']},{r['fwd_mm']},{r['fwd_3p']},{r['fwd_last']},"
                         f"{r['rev_mm']},{r['rev_3p']},{r['rev_last']},{'yes' if r['amplifies'] else 'no'}")
        save_text(self, "\n".join(lines) + "\n", "Export in-silico PCR", "CSV (*.csv)")
