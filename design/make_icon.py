"""Generate the NeoEdit app icon at all sizes (+ .ico / .icns).
Run: QT_QPA_PLATFORM=offscreen python design/make_icon.py
"""
import os, sys
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QImage, QPainter, QFont, QColor, QPen, QPainterPath
from PySide6.QtCore import Qt, QRectF

NT = {"A": "#008000", "C": "#0000ff", "G": "#000000", "T": "#ff0000"}
NCOL = "#aa00aa"                     # BioEdit's color for the degenerate base N
LETTERS = ["ATGC", "GCCT", "TAGC", "CCTA"]   # filler letters (cells not in the N)
NMASK = {(r, 0) for r in range(4)} | {(r, 3) for r in range(4)} | {(1, 1), (2, 2)}
OUT = os.path.join(os.path.dirname(__file__), "..", "src", "neoedit", "resources", "icons")


def render(S):
    im = QImage(S, S, QImage.Format_ARGB32); im.fill(Qt.transparent)
    p = QPainter(im); p.setRenderHint(QPainter.Antialiasing); p.setRenderHint(QPainter.TextAntialiasing)
    m = S * 6 / 256; rad = S * 40 / 256
    path = QPainterPath(); path.addRoundedRect(QRectF(m, m, S - 2 * m, S - 2 * m), rad, rad)
    p.fillPath(path, QColor("#ffffff")); p.setPen(QPen(QColor("#9a9a9a"), max(1, S * 3 / 256))); p.drawPath(path)
    x0 = y0 = S * 28 / 256; cell = (S - 2 * x0) / 4
    f = QFont("Courier New", max(4, int(cell * 0.62))); f.setBold(True); p.setFont(f); p.setClipPath(path)
    for r in range(4):
        for c in range(4):
            rect = QRectF(x0 + c * cell, y0 + r * cell, cell, cell)
            if (r, c) in NMASK:
                p.fillRect(rect, QColor(NCOL)); p.setPen(QColor("#ffffff")); p.drawText(rect, Qt.AlignCenter, "N")
            else:
                ch = LETTERS[r][c]; p.setPen(QColor(NT[ch])); p.drawText(rect, Qt.AlignCenter, ch)
    p.end(); return im


def main():
    app = QApplication.instance() or QApplication(sys.argv)
    os.makedirs(OUT, exist_ok=True)
    sizes = (16, 24, 32, 48, 64, 128, 256, 512, 1024)
    for S in sizes:
        render(S).save(os.path.join(OUT, f"neoedit_{S}.png"))
    from PIL import Image
    imgs = {S: Image.open(os.path.join(OUT, f"neoedit_{S}.png")) for S in sizes}
    imgs[256].save(os.path.join(OUT, "neoedit.ico"), sizes=[(s, s) for s in (16, 24, 32, 48, 64, 128, 256)])
    imgs[1024].save(os.path.join(OUT, "neoedit.icns"), format="ICNS",
                    sizes=[(s, s) for s in (16, 32, 64, 128, 256, 512, 1024)])
    imgs[256].save(os.path.join(OUT, "neoedit.png"))
    print("icons written to", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
