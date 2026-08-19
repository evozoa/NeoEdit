"""Window wrapping the circular genome map, with display options and export."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QCheckBox, QPushButton, QLabel,
                               QFileDialog, QMessageBox, QSpinBox)

from ..circular_view import CircularView


class CircularDialog(QDialog):
    positionClicked = Signal(int)
    geneActivated = Signal(object)

    def __init__(self, parent=None, seqid="", length=0, ann=None, fetch_seq=None, title=""):
        super().__init__(parent)
        self.setWindowTitle(f"Circular map — {title or seqid}")
        self.resize(720, 780)
        self.setModal(False)
        lay = QVBoxLayout(self)
        self.view = CircularView()
        self.view.set_data(seqid, length, ann, fetch_seq, title)
        self.view.positionClicked.connect(self.positionClicked)
        self.view.geneActivated.connect(self.geneActivated)
        lay.addWidget(self.view, 1)
        bar = QHBoxLayout()
        self.cb_gc = QCheckBox("GC content"); self.cb_gc.setChecked(True)
        self.cb_sk = QCheckBox("GC skew"); self.cb_sk.setChecked(True)
        self.cb_lbl = QCheckBox("Gene labels"); self.cb_lbl.setChecked(True)
        for cb, attr in ((self.cb_gc, "show_gc"), (self.cb_sk, "show_skew"), (self.cb_lbl, "show_labels")):
            cb.toggled.connect(lambda on, a=attr: (setattr(self.view, a, on), self.view.update()))
            bar.addWidget(cb)
        bar.addWidget(QLabel("  Export size:"))
        self.size = QSpinBox(); self.size.setRange(400, 8000); self.size.setValue(2000); self.size.setSingleStep(200)
        bar.addWidget(self.size)
        b_svg = QPushButton("Export SVG…"); b_svg.clicked.connect(lambda: self.export("svg"))
        b_png = QPushButton("Export PNG…"); b_png.clicked.connect(lambda: self.export("png"))
        b_close = QPushButton("Close"); b_close.clicked.connect(self.close)
        bar.addStretch(); bar.addWidget(b_svg); bar.addWidget(b_png); bar.addWidget(b_close)
        lay.addLayout(bar)
        lay.addWidget(QLabel("Click a gene or anywhere on the ring to move the alignment there; double-click a gene to open it."))

    def export(self, kind: str):
        filt = "SVG (*.svg)" if kind == "svg" else "PNG (*.png)"
        path, _ = QFileDialog.getSaveFileName(self, "Export circular map", f"{self.view.seqid}.{kind}", filt)
        if not path:
            return
        if not path.lower().endswith("." + kind):
            path += "." + kind
        try:
            self.view.export_image(path, self.size.value())
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e)); return
        self.parent().statusBar().showMessage(f"Saved {path}", 5000) if self.parent() else None
