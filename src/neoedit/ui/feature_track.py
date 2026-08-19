"""Dock listing annotation features; click to jump."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QTableWidgetItem

from .dialogs.common import make_table, NumItem, save_text


class FeaturePanel(QWidget):
    featureSelected = Signal(object)
    featuresChanged = Signal()

    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.model = model
        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        self.table = make_table(["Seq", "Type", "Label", "Start", "End", "Strand"])
        self.table.itemSelectionChanged.connect(self._sel)
        lay.addWidget(self.table)
        row = QHBoxLayout()
        rm = QPushButton("Remove"); rm.clicked.connect(self.remove_selected)
        clr = QPushButton("Clear all"); clr.clicked.connect(self.clear_all)
        exp = QPushButton("Export GFF3"); exp.clicked.connect(self.export_gff)
        row.addWidget(rm); row.addWidget(clr); row.addWidget(exp); row.addStretch()
        lay.addLayout(row)

    def set_model(self, model):
        self.model = model
        self.refresh()

    def refresh(self):
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        for i, f in enumerate(self.model.features):
            self.table.insertRow(i)
            name = self.model.rows[f.row].name if f.row < self.model.nrows else "?"
            vals = [name, f.type, f.label, NumItem(f.start + 1), NumItem(f.end), "+" if f.strand > 0 else "-"]
            for j, v in enumerate(vals):
                it = v if isinstance(v, QTableWidgetItem) else QTableWidgetItem(str(v))
                it.setData(Qt.UserRole + 1, i)
                self.table.setItem(i, j, it)
        self.table.setSortingEnabled(True)

    def _sel(self):
        items = self.table.selectedItems()
        if items:
            self.featureSelected.emit(self.model.features[items[0].data(Qt.UserRole + 1)])

    def remove_selected(self):
        items = self.table.selectedItems()
        if not items:
            return
        idx = items[0].data(Qt.UserRole + 1)
        del self.model.features[idx]
        self.refresh()
        self.featuresChanged.emit()

    def clear_all(self):
        self.model.features.clear()
        self.refresh()
        self.featuresChanged.emit()

    def export_gff(self):
        lines = ["##gff-version 3"]
        for k, f in enumerate(self.model.features, 1):
            name = self.model.rows[f.row].name if f.row < self.model.nrows else "seq"
            # report ungapped coordinates
            seq = self.model.rows[f.row].seq if f.row < self.model.nrows else ""
            ug_start = sum(1 for c in seq[:f.start] if c not in "-.~") + 1
            ug_end = sum(1 for c in seq[:f.end] if c not in "-.~")
            lines.append("\t".join([name, "neoedit", f.type, str(ug_start), str(ug_end), ".",
                                    "+" if f.strand > 0 else "-", ".", f"ID=f{k};Name={f.label}"]))
        save_text(self, "\n".join(lines) + "\n", "Export features", "GFF3 (*.gff3 *.gff)")
