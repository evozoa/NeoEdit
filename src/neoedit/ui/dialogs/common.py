from __future__ import annotations

from PySide6.QtGui import QFont
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QComboBox, QTableWidget, QTableWidgetItem, QHeaderView, QFileDialog,
                               QMessageBox, QPlainTextEdit, QDialog, QVBoxLayout, QDialogButtonBox)

from ...analysis.translate import codon_tables


def codon_table_combo(default: int = 1) -> QComboBox:
    cb = QComboBox()
    for tid, name in codon_tables():
        cb.addItem(f"{tid}: {name}", tid)
    idx = cb.findData(default)
    cb.setCurrentIndex(max(0, idx))
    cb.setMaxVisibleItems(30)
    return cb


def make_table(headers: list[str]) -> QTableWidget:
    t = QTableWidget(0, len(headers))
    t.setHorizontalHeaderLabels(headers)
    t.setSelectionBehavior(QTableWidget.SelectRows)
    t.setSelectionMode(QTableWidget.SingleSelection)
    t.setEditTriggers(QTableWidget.NoEditTriggers)
    t.setSortingEnabled(True)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
    t.horizontalHeader().setStretchLastSection(True)
    t.verticalHeader().setVisible(False)
    return t


class NumItem(QTableWidgetItem):
    def __init__(self, value, fmt="{}"):
        super().__init__(fmt.format(value))
        self.setData(Qt.UserRole, value)

    def __lt__(self, other):
        a, b = self.data(Qt.UserRole), other.data(Qt.UserRole)
        try:
            return a < b
        except TypeError:
            return super().__lt__(other)


def save_text(parent, text: str, title: str, filt: str, default: str = ""):
    path, _ = QFileDialog.getSaveFileName(parent, title, default, filt)
    if not path:
        return None
    try:
        with open(path, "w") as fh:
            fh.write(text)
    except OSError as e:
        QMessageBox.critical(parent, "Save failed", str(e))
        return None
    return path


class TextDialog(QDialog):
    def __init__(self, parent, title, text, mono=True):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(700, 500)
        lay = QVBoxLayout(self)
        self.edit = QPlainTextEdit()
        self.edit.setPlainText(text)
        self.edit.setReadOnly(True)
        if mono:
            f = self.edit.font(); f.setFamily("DejaVu Sans Mono"); f.setStyleHint(QFont.Monospace); self.edit.setFont(f)
        lay.addWidget(self.edit)
        bb = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Close)
        bb.rejected.connect(self.reject)
        bb.button(QDialogButtonBox.Save).clicked.connect(
            lambda: save_text(self, self.edit.toPlainText(), "Save", "Text files (*.txt);;All files (*)"))
        lay.addWidget(bb)
