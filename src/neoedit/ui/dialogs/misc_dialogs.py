from __future__ import annotations

import numpy as np
from PySide6.QtGui import QFont
from PySide6.QtCore import Qt, Signal, QRectF
from PySide6.QtGui import QPainter, QColor, QPen, QPalette
from PySide6.QtWidgets import (QDialog, QFormLayout, QComboBox, QLineEdit, QCheckBox, QPushButton, QHBoxLayout,
                               QVBoxLayout, QLabel, QDialogButtonBox, QPlainTextEdit, QWidget, QFileDialog,
                               QTableWidgetItem, QGroupBox, QSpinBox, QDoubleSpinBox, QTabWidget, QScrollArea)

from ...analysis import translate as T
from ...analysis import external as EXT
from .common import make_table, NumItem, save_text


# ------------------------------------------------------------------ Find
class FindDialog(QDialog):
    findNext = Signal(str, bool, bool, bool)   # pattern, ignore_gaps, regex, case
    findInNames = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Find")
        self.setModal(False)
        lay = QVBoxLayout(self)
        form = QFormLayout()
        self.pattern = QLineEdit()
        self.ignore_gaps = QCheckBox("Ignore gaps in sequences"); self.ignore_gaps.setChecked(True)
        self.regex = QCheckBox("Regular expression")
        self.case = QCheckBox("Case sensitive")
        form.addRow("Find", self.pattern)
        form.addRow("", self.ignore_gaps); form.addRow("", self.regex); form.addRow("", self.case)
        lay.addLayout(form)
        row = QHBoxLayout()
        b = QPushButton("Find next"); b.setDefault(True)
        b.clicked.connect(lambda: self.findNext.emit(self.pattern.text(), self.ignore_gaps.isChecked(),
                                                     self.regex.isChecked(), self.case.isChecked()))
        n = QPushButton("Find in names"); n.clicked.connect(lambda: self.findInNames.emit(self.pattern.text()))
        c = QPushButton("Close"); c.clicked.connect(self.close)
        row.addWidget(b); row.addWidget(n); row.addWidget(c)
        lay.addLayout(row)


# ------------------------------------------------------------------ Stats
class StatsDialog(QDialog):
    def __init__(self, model, rows, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Sequence statistics")
        self.resize(700, 450)
        lay = QVBoxLayout(self)
        st = model.seq_type
        hdr = ["Name", "Length (aligned)", "Length (ungapped)", "Gaps", "GC%" if st != "protein" else "MW (Da)", "Composition"]
        tbl = make_table(hdr)
        for i, r in enumerate(rows):
            row = model.rows[r]
            ug = row.ungapped()
            tbl.insertRow(i)
            comp = T.composition(row.seq)
            comp_s = " ".join(f"{k}:{v}" for k, v in sorted(comp.items()))
            if st != "protein":
                metric = NumItem(100 * T.gc_content(row.seq), "{:.1f}")
            else:
                mw = T.molecular_weight(row.seq, st)
                metric = NumItem(mw or 0, "{:.0f}")
            vals = [row.name, NumItem(len(row.seq)), NumItem(len(ug)), NumItem(len(row.seq) - len(ug)), metric, comp_s]
            for j, v in enumerate(vals):
                tbl.setItem(i, j, v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v)))
        lay.addWidget(tbl)
        bb = QDialogButtonBox(QDialogButtonBox.Close); bb.rejected.connect(self.reject); lay.addWidget(bb)


# ------------------------------------------------------- Identity matrix
class IdentityDialog(QDialog):
    def __init__(self, model, rows, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Sequence identity matrix")
        self.resize(800, 500)
        lay = QVBoxLayout(self)
        names = [model.rows[r].name for r in rows]
        mat = T.identity_matrix([model.rows[r].seq for r in rows])
        tbl = make_table([""] + names)
        tbl.setSortingEnabled(False)
        for i, n in enumerate(names):
            tbl.insertRow(i)
            tbl.setItem(i, 0, QTableWidgetItem(n))
            for j in range(len(names)):
                it = QTableWidgetItem(f"{mat[i, j]:.3f}")
                v = mat[i, j]
                c = QColor(int(255 - 120 * v), int(255 - 40 * v), 255)
                it.setBackground(c)
                tbl.setItem(i, j + 1, it)
        lay.addWidget(tbl)
        row = QHBoxLayout()
        b = QPushButton("Save as CSV")
        def save():
            lines = ["," + ",".join(names)]
            for i, n in enumerate(names):
                lines.append(n + "," + ",".join(f"{mat[i, j]:.4f}" for j in range(len(names))))
            save_text(self, "\n".join(lines) + "\n", "Save identity matrix", "CSV (*.csv)")
        b.clicked.connect(save)
        row.addWidget(b); row.addStretch()
        bb = QDialogButtonBox(QDialogButtonBox.Close); bb.rejected.connect(self.reject); row.addWidget(bb)
        lay.addLayout(row)


# ------------------------------------------------------- Conservation plot
class PlotWidget(QWidget):
    def __init__(self, values: np.ndarray, title: str, ymax: float, parent=None):
        super().__init__(parent)
        self.values = values
        self.title = title
        self.ymax = ymax or 1.0
        self.setMinimumHeight(220)
        self.setMinimumWidth(max(400, len(values) * 2))

    def paintEvent(self, e):
        p = QPainter(self)
        pal = self.palette()
        p.fillRect(self.rect(), pal.color(QPalette.Base))
        W, H = self.width(), self.height()
        left, bottom, top = 40, 24, 24
        p.setPen(pal.color(QPalette.Text))
        p.drawText(left, 16, self.title)
        p.drawLine(left, top, left, H - bottom)
        p.drawLine(left, H - bottom, W - 4, H - bottom)
        n = len(self.values)
        if n == 0:
            return
        pw = (W - left - 8) / n
        ph = H - top - bottom
        p.setPen(Qt.NoPen)
        for i, v in enumerate(self.values):
            h = ph * float(v) / self.ymax
            t = float(v) / self.ymax
            col = QColor(int(60 + 180 * (1 - t)), int(90 + 120 * t), 200)
            p.setBrush(col)
            p.drawRect(QRectF(left + i * pw, H - bottom - h, max(1.0, pw), h))
        p.setPen(pal.color(QPalette.Text))
        p.drawText(2, top + 6, f"{self.ymax:.1f}")
        p.drawText(2, H - bottom, "0")
        for x in range(0, n, max(1, n // 10)):
            p.drawText(int(left + x * pw), H - 6, str(x + 1))


class PlotDialog(QDialog):
    def __init__(self, model, rows, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Conservation plots")
        self.resize(900, 560)
        seqs = [model.rows[r].seq for r in rows]
        lay = QVBoxLayout(self)
        tabs = QTabWidget()
        ent = T.column_entropy(seqs)
        ident = T.column_identity(seqs)
        for title, vals, ymax in (("Entropy (bits)", ent, max(0.01, float(ent.max()) if len(ent) else 1)),
                                  ("Identity (fraction)", ident, 1.0)):
            sa = QScrollArea(); sa.setWidgetResizable(False)
            sa.setWidget(PlotWidget(vals, title, ymax))
            tabs.addTab(sa, title)
        lay.addWidget(tabs)
        bb = QDialogButtonBox(QDialogButtonBox.Close); bb.rejected.connect(self.reject); lay.addWidget(bb)


# ------------------------------------------------------- MAFFT
class AlignDialog(QDialog):
    def __init__(self, parent=None, settings=None):
        super().__init__(parent)
        self.setWindowTitle("Align with MAFFT")
        self.settings = settings
        lay = QVBoxLayout(self)
        exe = EXT.find_mafft(settings.value("exe/MAFFT") if settings else None)
        status = QLabel(f"MAFFT: {exe}  ({EXT.mafft_version(exe)})" if exe else
                        "MAFFT not found.  " + EXT.mafft_install_hint() + "\nSet its location in Edit > Preferences.")
        status.setWordWrap(True)
        lay.addWidget(status)
        form = QFormLayout()
        self.strategy = QComboBox()
        for label, _ in EXT.MAFFT_STRATEGIES:
            self.strategy.addItem(label)
        self.strategy.setCurrentIndex(int(settings.value("mafft/strategy", 0)) if settings else 0)
        self.threads = QSpinBox(); self.threads.setRange(0, 256); self.threads.setSpecialValueText("all cores")
        self.threads.setValue(int(settings.value("mafft/threads", 0)) if settings else 0)
        self.adjust = QCheckBox("Adjust direction (reverse-complement sequences as needed; nucleotide only)")
        self.adjust.setChecked(bool(settings and settings.value("mafft/adjust", False) in (True, "true")))
        self.keep_order = QCheckBox("Keep input order"); self.keep_order.setChecked(True)
        self.extra = QLineEdit()
        self.scope = QComboBox(); self.scope.addItems(["All sequences", "Selected sequences only"])
        form.addRow("Strategy", self.strategy)
        form.addRow("Threads", self.threads)
        form.addRow("", self.adjust)
        form.addRow("", self.keep_order)
        form.addRow("Extra MAFFT arguments", self.extra)
        form.addRow("Scope", self.scope)
        lay.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Align")
        bb.button(QDialogButtonBox.Ok).setEnabled(bool(exe))
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject); lay.addWidget(bb)

    def accept(self):
        if self.settings:
            self.settings.setValue("mafft/strategy", self.strategy.currentIndex())
            self.settings.setValue("mafft/threads", self.threads.value())
            self.settings.setValue("mafft/adjust", self.adjust.isChecked())
        super().accept()

    def values(self):
        return dict(strategy=self.strategy.currentIndex(), threads=self.threads.value(),
                    adjust_direction=self.adjust.isChecked(), keep_order=self.keep_order.isChecked(),
                    extra_args=self.extra.text().split(), sel_only=self.scope.currentIndex() == 1)


# ------------------------------------------------------- Preferences
class PreferencesDialog(QDialog):
    def __init__(self, settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.settings = settings
        lay = QVBoxLayout(self)
        grp = QGroupBox("MAFFT")
        form = QFormLayout(grp)
        self.edits = {}
        row = QHBoxLayout()
        ed = QLineEdit(settings.value("exe/MAFFT", ""))
        ed.setPlaceholderText("leave blank to search PATH")
        b = QPushButton("…"); b.clicked.connect(lambda _, e=ed: self._browse(e))
        row.addWidget(ed); row.addWidget(b)
        w = QWidget(); w.setLayout(row)
        form.addRow("MAFFT executable", w)
        self.edits["MAFFT"] = ed
        found = EXT.find_mafft(settings.value("exe/MAFFT") or None)
        hint = QLabel((f"Found: {found}  ({EXT.mafft_version(found)})" if found else
                       "Not found. To install: " + EXT.mafft_install_hint()))
        hint.setWordWrap(True)
        form.addRow("", hint)
        lay.addWidget(grp)
        g2 = QGroupBox("Defaults")
        f2 = QFormLayout(g2)
        self.table = QSpinBox(); self.table.setRange(1, 33); self.table.setValue(int(settings.value("default_table", 1)))
        f2.addRow("Default genetic code id", self.table)
        lay.addWidget(g2)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject); lay.addWidget(bb)

    def _browse(self, ed):
        p, _ = QFileDialog.getOpenFileName(self, "Select executable")
        if p:
            ed.setText(p)

    def accept(self):
        for name, ed in self.edits.items():
            self.settings.setValue(f"exe/{name}", ed.text())
        self.settings.setValue("default_table", self.table.value())
        super().accept()


# ------------------------------------------------------- New sequence
class NewSequenceDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("New sequence")
        self.resize(600, 400)
        lay = QVBoxLayout(self)
        form = QFormLayout()
        self.name = QLineEdit("new_seq")
        form.addRow("Name", self.name)
        lay.addLayout(form)
        lay.addWidget(QLabel("Sequence (raw or FASTA; whitespace and digits are stripped):"))
        self.text = QPlainTextEdit()
        f = self.text.font(); f.setFamily("DejaVu Sans Mono"); f.setStyleHint(QFont.Monospace); self.text.setFont(f)
        lay.addWidget(self.text)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject); lay.addWidget(bb)

    def records(self):
        """-> list of (name, seq)"""
        import re
        t = self.text.toPlainText().strip()
        if t.startswith(">"):
            out, name, buf = [], None, []
            for ln in t.splitlines():
                if ln.startswith(">"):
                    if name is not None:
                        out.append((name, "".join(buf)))
                    name, buf = ln[1:].split()[0] if ln[1:].split() else "seq", []
                else:
                    buf.append(re.sub(r"[\s\d]", "", ln))
            if name is not None:
                out.append((name, "".join(buf)))
            return out
        return [(self.name.text() or "new_seq", re.sub(r"[\s\d]", "", t))]


# ------------------------------------------------------------------ Consensus
class ConsensusDialog(QDialog):
    """Build a consensus of all / selected rows: majority (plurality) with a threshold, or an
    IUPAC-degenerate consensus; add it to the alignment as a row, copy it, or save it."""

    def __init__(self, model, selected_rows, threshold=0.5, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Consensus sequence")
        self.resize(720, 480)
        self.model = model
        self.selected = [r for r in selected_rows if 0 <= r < model.nrows]
        lay = QVBoxLayout(self)
        form = QFormLayout()
        self.scope = QComboBox()
        self.scope.addItem(f"All sequences ({model.nrows})", "all")
        if 1 < len(self.selected) < model.nrows:
            self.scope.addItem(f"Selected sequences ({len(self.selected)})", "sel")
            self.scope.setCurrentIndex(1)
        form.addRow("Sequences", self.scope)
        self.method = QComboBox()
        self.method.addItem("Majority / plurality (threshold, else N/X)", "majority")
        if model.is_nucleotide():
            self.method.addItem("IUPAC degenerate (every base above a minimum fraction)", "iupac")
        form.addRow("Method", self.method)
        self.thr = QDoubleSpinBox(); self.thr.setRange(0.0, 1.0); self.thr.setSingleStep(0.05); self.thr.setDecimals(2)
        self.thr.setValue(float(threshold))
        self.thr_label = QLabel()
        row = QHBoxLayout(); row.addWidget(self.thr); row.addWidget(self.thr_label, 1)
        form.addRow("Threshold", row)
        self.ignore_gaps = QCheckBox("Ignore gaps (fractions over residues only; a gap-only column gives '-')")
        self.ignore_gaps.setChecked(True)
        form.addRow("", self.ignore_gaps)
        self.plurality = QCheckBox("Accept the commonest residue when it reaches 50 % even below the threshold")
        self.plurality.setChecked(True)
        form.addRow("", self.plurality)
        self.name = QLineEdit("Consensus")
        form.addRow("Name", self.name)
        lay.addLayout(form)
        self.preview = QPlainTextEdit(); self.preview.setReadOnly(True)
        f = self.preview.font(); f.setFamily("DejaVu Sans Mono"); f.setStyleHint(QFont.Monospace); self.preview.setFont(f)
        lay.addWidget(self.preview, 1)
        self.info = QLabel(""); lay.addWidget(self.info)
        bb = QDialogButtonBox()
        b_add = bb.addButton("Add as sequence", QDialogButtonBox.ActionRole)
        b_copy = bb.addButton("Copy", QDialogButtonBox.ActionRole)
        b_save = bb.addButton("Save FASTA…", QDialogButtonBox.ActionRole)
        bb.addButton(QDialogButtonBox.Close)
        b_add.clicked.connect(self.add_row); b_copy.clicked.connect(self.copy); b_save.clicked.connect(self.save)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        for wdg in (self.scope, self.method):
            wdg.currentIndexChanged.connect(self.recompute)
        self.thr.valueChanged.connect(self.recompute)
        self.ignore_gaps.toggled.connect(self.recompute); self.plurality.toggled.connect(self.recompute)
        self.recompute()

    def rows(self):
        idx = self.selected if self.scope.currentData() == "sel" else range(self.model.nrows)
        return [self.model.rows[i].seq for i in idx]

    def consensus(self) -> str:
        rows = self.rows()
        if self.method.currentData() == "iupac":
            from ...analysis import primer_design as PD
            return PD.degenerate_consensus(rows, threshold=self.thr.value())
        return T.consensus(rows, threshold=self.thr.value(), ignore_gaps=self.ignore_gaps.isChecked(),
                           plurality=self.plurality.isChecked())

    def recompute(self, *_a):
        iupac = self.method.currentData() == "iupac"
        self.thr_label.setText("minimum fraction of a base to be included in the IUPAC code" if iupac
                               else "fraction of residues that must agree")
        self.ignore_gaps.setEnabled(not iupac); self.plurality.setEnabled(not iupac)
        cons = self.consensus()
        self._cons = cons
        self.preview.setPlainText("\n".join(cons[i:i + 100] for i in range(0, len(cons), 100)))
        n_amb = sum(1 for ch in cons if ch not in "ACGTU-" and ch != "X") if self.model.is_nucleotide() else cons.count("X")
        self.info.setText(f"{len(cons):,} columns, {len(self.rows())} sequences; "
                          f"{n_amb:,} ambiguous position(s)" + (" (N = no residue reached the threshold)" if not iupac else ""))

    def add_row(self):
        from ...model import SequenceRow
        self.model.add_row(SequenceRow(self.name.text().strip() or "Consensus", self._cons,
                                       f"consensus ({self.method.currentText().split(' (')[0]}, threshold {self.thr.value():.2f})"))
        self.info.setText(f"Added '{self.name.text().strip() or 'Consensus'}' as the last row.")

    def copy(self):
        from PySide6.QtWidgets import QApplication
        QApplication.clipboard().setText(f">{self.name.text().strip() or 'Consensus'}\n{self._cons}\n")
        self.info.setText("Copied as FASTA.")

    def save(self):
        save_text(self, f">{self.name.text().strip() or 'Consensus'}\n{self._cons}\n", "Save consensus",
                  "FASTA (*.fasta *.fa);;All files (*)", (self.name.text().strip() or "consensus") + ".fasta")


# ------------------------------------------------------------------ Set origin (circular)
class OriginDialog(QDialog):
    """Choose where a circular molecule should start: a position, the start of an annotated
    gene / feature (so it is not split across the origin), or the widest feature-free gap."""
    STANDARD = ("trnf", "trn-f", "trna-phe", "mt-tf", "tRNA-Phe".lower(), "trnp", "phe")

    def __init__(self, length: int, features, cursor_pos: int = 0, parent=None, circular: bool = True):
        """features: list of (label, start, end, strand) in 0-based reference coordinates (unwrapped)."""
        super().__init__(parent)
        self.setWindowTitle("Set origin (rotate circular molecule)")
        self.resize(560, 360)
        self.length = max(1, length)
        self.features = list(features)
        lay = QVBoxLayout(self)
        info = QLabel("Rotate the molecule so the chosen point becomes position 1. Sequences, grid features, "
                      "gene models and ORF tracks rotate together; a gene placed at the origin is no longer split "
                      "across it." + ("" if circular else "<br><b>The molecule is not marked circular</b> — rotating a linear sequence rearranges it."))
        info.setWordWrap(True); lay.addWidget(info)
        form = QFormLayout()
        from PySide6.QtWidgets import QRadioButton, QButtonGroup
        self.r_pos = QRadioButton("Position"); self.r_feat = QRadioButton("Start of feature"); self.r_gap = QRadioButton("Largest gap between features")
        grp = QButtonGroup(self)
        for r in (self.r_pos, self.r_feat, self.r_gap):
            grp.addButton(r)
        self.pos = QSpinBox(); self.pos.setRange(1, self.length); self.pos.setValue(min(self.length, max(1, cursor_pos + 1)))
        self.pos.setSuffix(f"  (1–{self.length:,})")
        form.addRow(self.r_pos, self.pos)
        self.feat = QComboBox()
        std = -1
        for i, (lab, a, b, st) in enumerate(sorted(self.features, key=lambda f: f[1])):
            from ...genome.annotations import fmt_span
            self.feat.addItem(f"{lab}   {fmt_span(a, b, self.length)} ({'+' if st > 0 else '-'})", i)
            if std < 0 and lab.lower().replace("_", "-") in self.STANDARD or lab.lower().startswith(("trnf", "mt-tf", "trna-phe")):
                std = self.feat.count() - 1
        self._sorted = sorted(self.features, key=lambda f: f[1])
        if std >= 0:
            self.feat.setCurrentIndex(std)
        form.addRow(self.r_feat, self.feat)
        form.addRow(self.r_gap, QLabel("origin placed so that no annotated feature crosses it"))
        lay.addLayout(form)
        self.flip = QCheckBox("Reverse-complement as well, so the chosen feature reads forward (minus-strand genes)")
        lay.addWidget(self.flip)
        self.note = QLabel("Vertebrate mitogenomes conventionally start at tRNA-Phe (trnF), followed by 12S rRNA."
                           if std >= 0 else "")
        self.note.setStyleSheet("color: gray"); self.note.setWordWrap(True); lay.addWidget(self.note)
        (self.r_feat if self.features else self.r_pos).setChecked(True)
        self.feat.setEnabled(bool(self.features)); self.r_feat.setEnabled(bool(self.features)); self.r_gap.setEnabled(bool(self.features))
        self.feat.currentIndexChanged.connect(self._feat_changed); self._feat_changed()
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject); lay.addWidget(bb)

    def _feat_changed(self, *_a):
        i = self.feat.currentData()
        if i is None or not self._sorted:
            return
        _lab, _a, _b, st = self._sorted[i]
        self.flip.setChecked(st < 0)

    def values(self) -> tuple[int, bool]:
        """-> (0-based reference position that becomes position 1, reverse-complement?)"""
        from ...genome.annotations import largest_gap_origin
        flip = self.flip.isChecked()
        if self.r_feat.isChecked() and self._sorted:
            _lab, a, b, st = self._sorted[self.feat.currentData()]
            # after rotation + flip the feature occupies [0, len): rotate to its 5' end
            return ((b % self.length) if flip else (a % self.length)), flip
        if self.r_gap.isChecked() and self._sorted:
            o = largest_gap_origin([(a, b) for _l, a, b, _s in self._sorted], self.length)
            return (o if o is not None else 0), flip
        return self.pos.value() - 1, flip
