from __future__ import annotations

import numpy as np
from PySide6.QtGui import QFont
from PySide6.QtCore import Qt, Signal, QRectF
from PySide6.QtGui import QPainter, QColor, QPen, QPalette
from PySide6.QtWidgets import (QDialog, QFormLayout, QComboBox, QLineEdit, QCheckBox, QPushButton, QHBoxLayout,
                               QVBoxLayout, QLabel, QDialogButtonBox, QPlainTextEdit, QWidget, QFileDialog,
                               QTableWidgetItem, QGroupBox, QSpinBox, QTabWidget, QScrollArea)

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
