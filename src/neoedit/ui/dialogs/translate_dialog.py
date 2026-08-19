from __future__ import annotations

from PySide6.QtWidgets import (QDialog, QFormLayout, QComboBox, QDialogButtonBox, QCheckBox, QVBoxLayout, QLabel)

from .common import codon_table_combo


class TranslateDialog(QDialog):
    """Options for translating selected sequences."""

    def __init__(self, parent=None, table=1):
        super().__init__(parent)
        self.setWindowTitle("Translate")
        lay = QVBoxLayout(self)
        form = QFormLayout()
        self.frame = QComboBox()
        self.frame.addItems(["+1", "+2", "+3", "-1", "-2", "-3", "All six frames"])
        self.table = codon_table_combo(table)
        self.output = QComboBox()
        self.output.addItems(["Open in new window", "Show report", "Replace selected sequences (in place)"])
        self.to_stop = QCheckBox("Stop at first stop codon")
        self.keep_align = QCheckBox("Keep alignment (translate codon columns, gaps -> '-')")
        self.keep_align.setChecked(True)
        form.addRow("Frame", self.frame)
        form.addRow("Genetic code", self.table)
        form.addRow("Output", self.output)
        form.addRow("", self.to_stop)
        form.addRow("", self.keep_align)
        lay.addLayout(form)
        lay.addWidget(QLabel("Tip: View > Show translation overlays a live translation under each DNA row."))
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def values(self):
        return dict(frame=self.frame.currentIndex(), table=self.table.currentData(),
                    output=self.output.currentIndex(), to_stop=self.to_stop.isChecked(),
                    keep_align=self.keep_align.isChecked())
