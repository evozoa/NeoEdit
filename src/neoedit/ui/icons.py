"""Programmatically drawn toolbar icons (no external image files)."""
from __future__ import annotations

from PySide6.QtCore import Qt, QRect, QPointF, QRectF, QSize
from PySide6.QtGui import QIcon, QPixmap, QPainter, QColor, QPen, QFont, QBrush, QPolygonF, QPainterPath

SZ = 22
NT = {"A": "#008000", "C": "#0000ff", "G": "#000000", "T": "#ff0000"}
_cache: dict[str, QIcon] = {}


def _pix():
    pm = QPixmap(SZ, SZ)
    pm.setDevicePixelRatio(1)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    p.setRenderHint(QPainter.TextAntialiasing)
    return pm, p


def _font(size=7, bold=True, family="Courier New"):
    f = QFont(family, size)
    f.setBold(bold)
    return f


def _letters_square(inverse: bool):
    pm, p = _pix()
    cell = 9
    x0, y0 = 2, 2
    p.setPen(QPen(QColor("#888888"), 1))
    p.setBrush(QColor("#ffffff"))
    p.drawRect(x0, y0, 2 * cell, 2 * cell)
    p.setFont(_font(8))
    for i, (ch, col) in enumerate(NT.items()):
        cx, cy = x0 + (i % 2) * cell, y0 + (i // 2) * cell
        r = QRect(cx, cy, cell, cell)
        if inverse:
            p.fillRect(r, QColor(col))
            p.setPen(QColor("#ffffff"))
        else:
            p.setPen(QColor(col))
        p.drawText(r, Qt.AlignCenter, ch)
    p.end()
    return pm


def _pointer():
    pm, p = _pix()
    poly = QPolygonF([QPointF(5, 3), QPointF(5, 17), QPointF(9, 13.5), QPointF(12, 19), QPointF(14.5, 18), QPointF(11.5, 12.5), QPointF(16.5, 12)])
    p.setPen(QPen(QColor("#222"), 1)); p.setBrush(QColor("#fff")); p.drawPolygon(poly)
    p.end(); return pm


def _ibeam_insert():
    pm, p = _pix()
    p.setPen(QPen(QColor("#c02020"), 2))
    p.drawLine(8, 4, 8, 18); p.drawLine(5, 4, 11, 4); p.drawLine(5, 18, 11, 18)
    p.setPen(QColor("#222")); p.setFont(_font(9, family="DejaVu Sans")); p.drawText(QRect(10, 3, 12, 16), Qt.AlignCenter, "+")
    p.end(); return pm


def _overwrite():
    pm, p = _pix()
    p.setPen(QPen(QColor("#c02020"), 1.5)); p.setBrush(QColor("#ffe0e0")); p.drawRect(4, 4, 14, 14)
    p.setPen(QColor("#222")); p.setFont(_font(10)); p.drawText(QRect(4, 4, 14, 14), Qt.AlignCenter, "A")
    p.end(); return pm


def _zoom(plus: bool):
    pm, p = _pix()
    p.setPen(QPen(QColor("#333"), 2)); p.setBrush(Qt.NoBrush); p.drawEllipse(3, 3, 11, 11); p.drawLine(13, 13, 19, 19)
    p.setPen(QPen(QColor("#333"), 1.6)); p.drawLine(6, 8.5, 11, 8.5)
    if plus:
        p.drawLine(8.5, 6, 8.5, 11)
    p.end(); return pm


def _open():
    pm, p = _pix()
    p.setPen(QPen(QColor("#8a6d1f"), 1)); p.setBrush(QColor("#f2c94c"))
    p.drawRoundedRect(2, 6, 18, 12, 1.5, 1.5); p.drawRect(2, 4, 7, 3)
    p.setBrush(QColor("#fbe08a")); p.drawRect(2, 9, 18, 9)
    p.end(); return pm


def _save():
    pm, p = _pix()
    p.setPen(QPen(QColor("#234"), 1)); p.setBrush(QColor("#3b6ea5")); p.drawRoundedRect(3, 3, 16, 16, 1.5, 1.5)
    p.setBrush(QColor("#fff")); p.drawRect(6, 3, 10, 6); p.setBrush(QColor("#cfd8e3")); p.drawRect(6, 13, 10, 6)
    p.setBrush(QColor("#3b6ea5")); p.drawRect(12, 4, 2, 4)
    p.end(); return pm


def _arrow_curve(redo: bool):
    pm, p = _pix()
    path = QPainterPath(); path.moveTo(5, 14); path.cubicTo(5, 6, 16, 6, 16, 12)
    p.setPen(QPen(QColor("#333"), 2)); p.setBrush(Qt.NoBrush); p.drawPath(path)
    head = QPolygonF([QPointF(16, 16), QPointF(12.5, 11), QPointF(19.5, 11)])
    p.setPen(Qt.NoPen); p.setBrush(QColor("#333")); p.drawPolygon(head)
    p.end()
    if not redo:
        pm = pm.transformed(__import__("PySide6.QtGui", fromlist=["QTransform"]).QTransform().scale(-1, 1))
    return pm


def _translation():
    pm, p = _pix()
    p.setFont(_font(6)); p.setPen(QColor("#008000")); p.drawText(QRect(1, 1, 20, 9), Qt.AlignCenter, "ATG")
    p.setPen(QPen(QColor("#666"), 1)); p.drawLine(4, 11, 18, 11)
    p.setPen(QColor("#222")); p.setFont(_font(8)); p.drawText(QRect(1, 11, 20, 10), Qt.AlignCenter, "M")
    p.end(); return pm


def _orf():
    pm, p = _pix()
    p.setPen(QPen(QColor("#555"), 1.5)); p.drawLine(2, 16, 20, 16)
    p.setPen(Qt.NoPen); p.setBrush(QColor("#f59e0b"))
    p.drawPolygon(QPolygonF([QPointF(3, 5), QPointF(14, 5), QPointF(18, 9), QPointF(14, 13), QPointF(3, 13)]))
    p.setPen(QColor("#222")); p.setFont(_font(6)); p.drawText(QRect(3, 5, 12, 8), Qt.AlignCenter, "ORF")
    p.end(); return pm


def _primer():
    pm, p = _pix()
    p.setPen(QPen(QColor("#888"), 1.2)); p.drawLine(2, 11, 20, 11)
    p.setPen(QPen(QColor("#10b981"), 2.2)); p.drawLine(2, 6, 11, 6); p.drawLine(8, 3.5, 11, 6); p.drawLine(8, 8.5, 11, 6)
    p.setPen(QPen(QColor("#ef4444"), 2.2)); p.drawLine(20, 16, 11, 16); p.drawLine(14, 13.5, 11, 16); p.drawLine(14, 18.5, 11, 16)
    p.end(); return pm


def _align():
    pm, p = _pix()
    cols = ["#008000", "#0000ff", "#ff0000", "#000000"]
    for i, (x, w) in enumerate([(3, 14), (6, 11), (4, 13)]):
        for k in range(4):
            p.fillRect(QRectF(x + k * (w / 4), 4 + i * 5, w / 4 - 0.6, 3.6), QColor(cols[(k + i) % 4]))
    p.end(); return pm


def _features():
    pm, p = _pix()
    p.setPen(QPen(QColor("#444"), 1)); p.setBrush(QColor("#f59e0b"))
    p.drawPolygon(QPolygonF([QPointF(3, 4), QPointF(13, 4), QPointF(19, 11), QPointF(13, 18), QPointF(3, 18)]))
    p.setBrush(QColor("#fff")); p.drawEllipse(QPointF(7, 11), 2, 2)
    p.end(); return pm


def _spacing(kind: str):
    pm, p = _pix()
    p.setPen(QPen(QColor("#333"), 1.5))
    if kind in ("row+", "row-"):
        for y in (5, 11, 17):
            p.drawLine(4, y, 18, y)
    else:
        for x in (5, 11, 17):
            p.drawLine(x, 4, x, 18)
    p.setPen(Qt.NoPen); p.setBrush(QColor("#ffffff")); p.drawEllipse(QPointF(16, 16), 5.5, 5.5)
    p.setPen(QPen(QColor("#c02020"), 2)); p.drawLine(12.5, 16, 19.5, 16)
    if kind.endswith("+"):
        p.drawLine(16, 12.5, 16, 19.5)
    p.end(); return pm


_BUILDERS = {
    "normal_view": lambda: _letters_square(False),
    "inverse_view": lambda: _letters_square(True),
    "mode_select": _pointer,
    "mode_insert": _ibeam_insert,
    "mode_overwrite": _overwrite,
    "zoom_in": lambda: _zoom(True),
    "zoom_out": lambda: _zoom(False),
    "open": _open,
    "save": _save,
    "undo": lambda: _arrow_curve(False),
    "redo": lambda: _arrow_curve(True),
    "translation": _translation,
    "orf": _orf,
    "primer": _primer,
    "align": _align,
    "features": _features,
    "row_more": lambda: _spacing("row+"),
    "row_less": lambda: _spacing("row-"),
    "col_more": lambda: _spacing("col+"),
    "col_less": lambda: _spacing("col-"),
}


def icon(name: str) -> QIcon:
    if name not in _cache:
        _cache[name] = QIcon(_BUILDERS[name]())
    return _cache[name]
