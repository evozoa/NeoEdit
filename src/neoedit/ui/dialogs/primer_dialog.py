from __future__ import annotations

from PySide6.QtGui import QFont
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QDialog, QFormLayout, QComboBox, QSpinBox, QDoubleSpinBox, QPushButton, QHBoxLayout,
                               QVBoxLayout, QLabel, QGroupBox, QTableWidgetItem, QMessageBox, QSplitter,
                               QPlainTextEdit, QGridLayout, QLineEdit, QCheckBox)

from ...analysis import primers as P
from ...analysis.orf_finder import gap_map
from .common import make_table, NumItem, save_text


class PrimerDialog(QDialog):
    pairSelected = Signal(int, int, int)     # row, start, end (gapped coords spanning product)
    featuresReady = Signal(list)

    def __init__(self, model, rows: list[int], sel_cols: tuple[int, int] | None, parent=None, circular: bool = False):
        super().__init__(parent)
        self.setWindowTitle("Primer design (Primer3)")
        self.resize(950, 620)
        self.setModal(False)
        self.model = model
        self.rows = rows
        self.pairs: list[P.PrimerPair] = []
        self.template_row = rows[0] if rows else 0

        lay = QVBoxLayout(self)
        params = QGroupBox("Parameters")
        g = QGridLayout(params)
        self.seq_combo = QComboBox()
        for r in rows:
            self.seq_combo.addItem(model.rows[r].name, r)
        g.addWidget(QLabel("Template"), 0, 0); g.addWidget(self.seq_combo, 0, 1, 1, 3)

        # target region (ungapped coords, 1-based in UI)
        gm = gap_map(model.rows[self.template_row].seq) if rows else []
        tgt_start, tgt_len = 0, 0
        if sel_cols and gm:
            # convert gapped selection to ungapped
            ug = [i for i, g_ in enumerate(gm) if sel_cols[0] <= g_ <= sel_cols[1]]
            if ug:
                tgt_start, tgt_len = ug[0], len(ug)
        self.use_target = QCheckBox("Product must span target region"); self.use_target.setChecked(bool(tgt_len))
        self.t_start = QSpinBox(); self.t_start.setRange(1, 10**7); self.t_start.setValue(tgt_start + 1)
        self.t_len = QSpinBox(); self.t_len.setRange(1, 10**7); self.t_len.setValue(max(1, tgt_len))
        g.addWidget(self.use_target, 1, 0); g.addWidget(QLabel("start"), 1, 1); g.addWidget(self.t_start, 1, 2)
        g.addWidget(QLabel("length"), 1, 3); g.addWidget(self.t_len, 1, 4)
        self.circ_cb = QCheckBox("Circular template: primers and products may span the origin")
        self.circ_cb.setChecked(bool(circular))
        self.circ_cb.setToolTip("Follows the window's topology; positions past the origin are shown as 'a-b (across origin)'")
        g.addWidget(self.circ_cb, 1, 5, 1, 2)

        def dsb(lo, hi, v, step=0.5):
            w = QDoubleSpinBox(); w.setRange(lo, hi); w.setValue(v); w.setSingleStep(step); return w

        def sb(lo, hi, v):
            w = QSpinBox(); w.setRange(lo, hi); w.setValue(v); return w

        self.prod_min = sb(20, 100000, 100); self.prod_max = sb(20, 100000, 300)
        self.size_min = sb(10, 40, 18); self.size_opt = sb(10, 40, 20); self.size_max = sb(10, 40, 25)
        self.tm_min = dsb(30, 90, 57); self.tm_opt = dsb(30, 90, 60); self.tm_max = dsb(30, 90, 63)
        self.gc_min = dsb(0, 100, 30, 5); self.gc_max = dsb(0, 100, 70, 5)
        self.num = sb(1, 50, 5)
        self.gc_clamp = sb(0, 5, 0)
        self.poly_x = sb(1, 10, 4)
        self.salt = dsb(0, 500, 50, 5); self.dna = dsb(0, 1000, 50, 10); self.mg = dsb(0, 20, 1.5, 0.1); self.dntp = dsb(0, 10, 0.6, 0.1)

        r = 2
        g.addWidget(QLabel("Product size"), r, 0); g.addWidget(self.prod_min, r, 1); g.addWidget(QLabel("–"), r, 2); g.addWidget(self.prod_max, r, 3)
        r += 1
        g.addWidget(QLabel("Primer size min/opt/max"), r, 0); g.addWidget(self.size_min, r, 1); g.addWidget(self.size_opt, r, 2); g.addWidget(self.size_max, r, 3)
        r += 1
        g.addWidget(QLabel("Tm min/opt/max (°C)"), r, 0); g.addWidget(self.tm_min, r, 1); g.addWidget(self.tm_opt, r, 2); g.addWidget(self.tm_max, r, 3)
        r += 1
        g.addWidget(QLabel("GC% min/max"), r, 0); g.addWidget(self.gc_min, r, 1); g.addWidget(self.gc_max, r, 3)
        r += 1
        g.addWidget(QLabel("Pairs to return"), r, 0); g.addWidget(self.num, r, 1)
        g.addWidget(QLabel("GC clamp"), r, 2); g.addWidget(self.gc_clamp, r, 3)
        g.addWidget(QLabel("Max poly-X"), r, 4); g.addWidget(self.poly_x, r, 5)
        r += 1
        g.addWidget(QLabel("Na+ mM / Mg2+ mM / dNTP mM / primer nM"), r, 0)
        g.addWidget(self.salt, r, 1); g.addWidget(self.mg, r, 2); g.addWidget(self.dntp, r, 3); g.addWidget(self.dna, r, 4)
        lay.addWidget(params)

        btns = QHBoxLayout()
        self.run_btn = QPushButton("Design primers"); self.run_btn.setDefault(True); self.run_btn.clicked.connect(self.run)
        self.feat_btn = QPushButton("Add selected pair as features"); self.feat_btn.clicked.connect(self.add_features)
        self.all_feat_btn = QPushButton("Add all pairs as features"); self.all_feat_btn.clicked.connect(lambda: self.add_features(all_pairs=True))
        self.csv_btn = QPushButton("Export CSV"); self.csv_btn.clicked.connect(self.export_csv)
        self.fa_btn = QPushButton("Export FASTA"); self.fa_btn.clicked.connect(self.export_fasta)
        self.mm_btn = QPushButton("Check against alignment"); self.mm_btn.clicked.connect(self.check_alignment)
        for b in (self.run_btn, self.feat_btn, self.all_feat_btn, self.csv_btn, self.fa_btn, self.mm_btn):
            btns.addWidget(b)
        btns.addStretch()
        lay.addLayout(btns)

        split = QSplitter(Qt.Vertical)
        self.table_w = make_table(["#", "Left primer (5'→3')", "L pos", "L Tm", "L GC%", "Right primer (5'→3')", "R pos", "R Tm", "R GC%",
                                   "Product", "Penalty", "Pair compl. (any/end)"])
        self.table_w.itemSelectionChanged.connect(self._on_select)
        split.addWidget(self.table_w)
        self.detail = QPlainTextEdit(); self.detail.setReadOnly(True)
        f = self.detail.font(); f.setFamily("DejaVu Sans Mono"); f.setStyleHint(QFont.Monospace); self.detail.setFont(f)
        split.addWidget(self.detail)
        split.setSizes([350, 180])
        lay.addWidget(split, 1)
        self.status = QLabel("")
        lay.addWidget(self.status)

    def run(self):
        self.template_row = self.seq_combo.currentData()
        tmpl = self.model.rows[self.template_row].seq
        target = (self.t_start.value() - 1, self.t_len.value()) if self.use_target.isChecked() else None
        try:
            self.pairs = P.design_primers(
                tmpl, target=target, product_range=((self.prod_min.value(), self.prod_max.value()),),
                opt_size=self.size_opt.value(), min_size=self.size_min.value(), max_size=self.size_max.value(),
                opt_tm=self.tm_opt.value(), min_tm=self.tm_min.value(), max_tm=self.tm_max.value(),
                min_gc=self.gc_min.value(), max_gc=self.gc_max.value(), num_return=self.num.value(),
                gc_clamp=self.gc_clamp.value(), max_poly_x=self.poly_x.value(),
                salt_monovalent=self.salt.value(), salt_divalent=self.mg.value(), dntp=self.dntp.value(), dna_conc=self.dna.value(),
                circular=self.circ_cb.isChecked())
        except ValueError as e:
            self.pairs = []
            self.table_w.setRowCount(0)
            self.status.setText(str(e)[:400])
            return
        except Exception as e:  # primer3 errors
            QMessageBox.critical(self, "Primer3 error", str(e))
            return
        self.table_w.setSortingEnabled(False)
        self.table_w.setRowCount(0)
        for i, pr in enumerate(self.pairs):
            self.table_w.insertRow(i)
            L = pr.circular_len
            rend = pr.right.start + pr.right.length
            vals = [NumItem(i + 1), pr.left.seq, NumItem(pr.left.start + 1), NumItem(pr.left.tm, "{:.1f}"), NumItem(pr.left.gc, "{:.0f}"),
                    pr.right.seq, NumItem(rend - L if (L and rend > L) else rend, "{}" + (" ⟳" if (L and rend > L) else "")),
                    NumItem(pr.right.tm, "{:.1f}"), NumItem(pr.right.gc, "{:.0f}"),
                    NumItem(pr.product_size, "{}" + (" ⟳" if pr.crosses_origin else "")), NumItem(pr.penalty, "{:.2f}"), f"{pr.compl_any:.1f}/{pr.compl_end:.1f}"]
            for j, v in enumerate(vals):
                item = v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v))
                item.setData(Qt.UserRole + 1, i)
                self.table_w.setItem(i, j, item)
        self.table_w.setSortingEnabled(True)
        self.table_w.sortItems(0, Qt.AscendingOrder)
        self.status.setText(f"{len(self.pairs)} pair(s)")
        if self.pairs:
            self.table_w.selectRow(0)

    def _current(self):
        items = self.table_w.selectedItems()
        if not items:
            return None, None
        i = items[0].data(Qt.UserRole + 1)
        return i, self.pairs[i]

    def _on_select(self):
        i, pr = self._current()
        if pr is None:
            return
        gm = gap_map(self.model.rows[self.template_row].seq)
        s = gm[pr.left.start] if pr.left.start < len(gm) else 0
        e_idx = pr.right.start + pr.right.length - 1
        if pr.crosses_origin:
            e = len(gm)                        # the product runs off the end and continues at column 1
        else:
            e = gm[e_idx] + 1 if e_idx < len(gm) else len(gm)
        self.pairSelected.emit(self.template_row, s, e)
        L, R = pr.left, pr.right
        CL = pr.circular_len
        self.detail.setPlainText(
            f"Pair {i + 1}   product {pr.product_size} bp{' (across the origin)' if pr.crosses_origin else ''}   penalty {pr.penalty:.2f}\n"
            f"Left : {L.seq}  pos {L.pos_text(CL)}  Tm {L.tm:.1f}  GC {L.gc:.0f}%  "
            f"hairpin {L.hairpin:.1f}  self-any {L.self_any:.1f}  self-end {L.self_end:.1f}  3' stability {L.end_stability:.1f}\n"
            f"Right: {R.seq}  pos {R.pos_text(CL)}  Tm {R.tm:.1f}  GC {R.gc:.0f}%  "
            f"hairpin {R.hairpin:.1f}  self-any {R.self_any:.1f}  self-end {R.self_end:.1f}  3' stability {R.end_stability:.1f}\n"
            f"Heterodimer Tm: {P.heterodimer_tm(L.seq, R.seq):.1f} °C")

    def add_features(self, all_pairs=False):
        if not self.pairs:
            return
        gm = gap_map(self.model.rows[self.template_row].seq)
        feats = []
        if all_pairs:
            for i, pr in enumerate(self.pairs):
                feats += P.pair_to_features(pr, self.template_row, gm, i + 1)
        else:
            i, pr = self._current()
            if pr is None:
                return
            feats = P.pair_to_features(pr, self.template_row, gm, i + 1)
        self.featuresReady.emit(feats)
        self.status.setText(f"Added {len(feats)} primer features")

    def export_csv(self):
        if not self.pairs:
            return
        lines = ["pair,left_seq,left_start,left_len,left_tm,left_gc,right_seq,right_start,right_len,right_tm,right_gc,product_size,penalty"]
        for i, pr in enumerate(self.pairs, 1):
            L, R = pr.left, pr.right
            lines.append(f"{i},{L.seq},{L.start + 1},{L.length},{L.tm:.2f},{L.gc:.1f},{R.seq},{R.start + R.length},{R.length},{R.tm:.2f},{R.gc:.1f},{pr.product_size},{pr.penalty:.3f}")
        save_text(self, "\n".join(lines) + "\n", "Export primers", "CSV (*.csv)")

    def export_fasta(self):
        if not self.pairs:
            return
        name = self.model.rows[self.template_row].name
        out = []
        for i, pr in enumerate(self.pairs, 1):
            out.append(f">{name}_F{i} Tm={pr.left.tm:.1f}\n{pr.left.seq}")
            out.append(f">{name}_R{i} Tm={pr.right.tm:.1f}\n{pr.right.seq}")
        save_text(self, "\n".join(out) + "\n", "Export primers", "FASTA (*.fasta *.fa)")

    def check_alignment(self):
        """Mismatch count of each primer vs every sequence in the alignment, using the
        aligned position of the primer in the template."""
        i, pr = self._current()
        if pr is None:
            return
        gm = gap_map(self.model.rows[self.template_row].seq)
        lines = [f"Pair {i + 1}: mismatches vs each aligned sequence (at template-aligned position)", ""]
        CL = pr.circular_len
        for label, prm, strand in (("Left ", pr.left, 1), ("Right", pr.right, -1)):
            pcs = [(gm[a], gm[b - 1] + 1) for a, b in prm.pieces(CL)]
            lines.append(f"{label} {prm.seq}  aligned cols " + " + ".join(f"{a + 1}-{b}" for a, b in pcs)
                         + (" (across origin)" if len(pcs) > 1 else ""))
            for r in range(self.model.nrows):
                site = "".join(self.model.rows[r].seq[a:b] for a, b in pcs)
                site_ug = "".join(c for c in site if c not in "-.~").upper()
                from Bio.Seq import Seq
                cmp_ = site_ug if strand > 0 else str(Seq(site_ug).reverse_complement())
                mm = sum(1 for a, b in zip(prm.seq.upper(), cmp_) if a != b and b != "N") + abs(len(prm.seq) - len(cmp_))
                mark = "".join("." if a == b else b.lower() if b else "-" for a, b in zip(prm.seq.upper(), cmp_.ljust(len(prm.seq))))
                lines.append(f"   {self.model.rows[r].name[:24].ljust(24)} {mark}  {mm} mm")
            lines.append("")
        self.detail.setPlainText("\n".join(lines))
