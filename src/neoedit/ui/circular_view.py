"""Circular genome map (mitogenome / plasmid): gene arrows, GC rings, focal wedge.

Used as an alternative tier-1 rendering and as a standalone figure exporter.
"""
from __future__ import annotations

import math
from typing import Optional

from PySide6.QtCore import Qt, QRectF, QPointF, Signal, QSize
from PySide6.QtGui import (QPainter, QColor, QPen, QBrush, QFont, QFontMetrics, QPainterPath, QPolygonF,
                           QPalette, QCursor)
from PySide6.QtWidgets import QWidget, QToolTip, QSizePolicy

from ..genome.annotations import Annotation, Gene

TWO_PI = 2 * math.pi


def gc_series(seq: str, n_bins: int) -> tuple[list[float], list[float]]:
    """(GC content, GC skew) per bin. Skew = (G-C)/(G+C)."""
    L = len(seq)
    if not L or n_bins <= 0:
        return [], []
    step = L / n_bins
    gc, skew = [], []
    for i in range(n_bins):
        a, b = int(i * step), int((i + 1) * step)
        s = seq[a:b].upper()
        g, c = s.count("G"), s.count("C")
        n = len(s) - s.count("-") - s.count("N")
        gc.append((g + c) / n if n else 0.0)
        skew.append((g - c) / (g + c) if (g + c) else 0.0)
    return gc, skew


def gene_color(g: Gene) -> QColor:
    bt = (g.biotype or "").lower()
    name = (g.name or "").upper()
    if getattr(g, "cytoplasmic", False):
        return QColor("#d926a9")           # MDP: translated on cytoplasmic ribosomes
    if "trna" in bt or name.startswith("TRN"):
        return QColor("#d98c00")
    if "rrna" in bt or name.startswith(("RRN", "S-RRNA", "L-RRNA", "12S", "16S")):
        return QColor("#8e44ad")
    if "d-loop" in bt or "control" in bt or "D-LOOP" in name:
        return QColor("#7f8c8d")
    if "cds" in bt or "protein" in bt or "mrna" in bt:
        return QColor("#1f5fbf") if g.strand > 0 else QColor("#c0392b")
    return QColor("#16a085")


class CircularView(QWidget):
    """Circular map. Angle 0 = top (12 o'clock), clockwise."""
    positionClicked = Signal(int)          # reference position
    geneActivated = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.length = 1
        self.seqid = ""
        self.ann: Optional[Annotation] = None
        self.fetch_seq = None
        self.focus: tuple[int, int] | None = None
        self.window: tuple[int, int] | None = None
        self.show_gc = True
        self.show_skew = True
        self.show_labels = True
        self.title = ""
        self._gc = None; self._skew = None
        self._hits: list[tuple[QPainterPath, Gene]] = []
        self.setMinimumSize(320, 320)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    # ---------------------------------------------------------------- data
    def set_data(self, seqid: str, length: int, ann: Optional[Annotation], fetch_seq=None, title: str = ""):
        self.seqid = seqid; self.length = max(1, length); self.ann = ann
        self.fetch_seq = fetch_seq; self.title = title or seqid
        self._gc = self._skew = None
        self.update()

    def set_focus(self, s, e):
        self.focus = (s, e); self.update()

    def set_window(self, s, e):
        self.window = (s, e); self.update()

    def _ensure_gc(self, bins: int):
        if self._gc is None and self.fetch_seq:
            seq = self.fetch_seq(0, self.length)
            self._gc, self._skew = gc_series(seq, bins)

    # ---------------------------------------------------------------- geometry
    def _angle(self, pos: float) -> float:
        return TWO_PI * (pos / self.length) - math.pi / 2      # 0 at top, clockwise

    def _pt(self, cx, cy, r, ang) -> QPointF:
        return QPointF(cx + r * math.cos(ang), cy + r * math.sin(ang))

    def _pos_from(self, cx, cy, x, y) -> int:
        ang = math.atan2(y - cy, x - cx) + math.pi / 2
        if ang < 0:
            ang += TWO_PI
        return int(self.length * (ang % TWO_PI) / TWO_PI)

    def _arc_path(self, cx, cy, r_in, r_out, a0, a1, arrow: int = 0) -> QPainterPath:
        """Annular sector; arrow=+1/-1 adds a point at the leading edge."""
        path = QPainterPath()
        if a1 < a0:
            a0, a1 = a1, a0
        span = a1 - a0
        head = min(span * 0.35, 0.05)
        rm = (r_in + r_out) / 2
        if arrow > 0:
            path.moveTo(self._pt(cx, cy, r_in, a0))
            path.arcTo(QRectF(cx - r_in, cy - r_in, 2 * r_in, 2 * r_in), -math.degrees(a0), -math.degrees(span - head))
            path.lineTo(self._pt(cx, cy, rm, a1))
            path.lineTo(self._pt(cx, cy, r_out, a1 - head))
            path.arcTo(QRectF(cx - r_out, cy - r_out, 2 * r_out, 2 * r_out), -math.degrees(a1 - head), math.degrees(span - head))
            path.closeSubpath()
        elif arrow < 0:
            path.moveTo(self._pt(cx, cy, rm, a0))
            path.lineTo(self._pt(cx, cy, r_in, a0 + head))
            path.arcTo(QRectF(cx - r_in, cy - r_in, 2 * r_in, 2 * r_in), -math.degrees(a0 + head), -math.degrees(span - head))
            path.lineTo(self._pt(cx, cy, r_out, a1))
            path.arcTo(QRectF(cx - r_out, cy - r_out, 2 * r_out, 2 * r_out), -math.degrees(a1), math.degrees(span - head))
            path.closeSubpath()
        else:
            path.moveTo(self._pt(cx, cy, r_in, a0))
            path.arcTo(QRectF(cx - r_in, cy - r_in, 2 * r_in, 2 * r_in), -math.degrees(a0), -math.degrees(span))
            path.lineTo(self._pt(cx, cy, r_out, a1))
            path.arcTo(QRectF(cx - r_out, cy - r_out, 2 * r_out, 2 * r_out), -math.degrees(a1), math.degrees(span))
            path.closeSubpath()
        return path

    # ---------------------------------------------------------------- paint
    def render_to(self, p: QPainter, W: int, H: int, for_export: bool = False):
        p.setRenderHint(QPainter.Antialiasing)
        if for_export:
            p.fillRect(0, 0, W, H, QColor("#ffffff"))
        cx, cy = W / 2, H / 2
        R = min(W, H) / 2 - 12
        f = QFont(self.font()); f.setPointSizeF(max(6.0, R * 0.045)); p.setFont(f)
        fm = QFontMetrics(f)
        label_gap = (max(fm.horizontalAdvance(g.name) for genes in self.ann.genes_by_seq.values() for g in genes)
                     if (self.show_labels and self.ann and self.ann.count()) else 0)
        label_gap = min(label_gap + 8, R * 0.30) if self.show_labels else 0
        r_out = R - label_gap
        ring_w = max(8.0, r_out * 0.10)
        r_fwd_out, r_fwd_in = r_out, r_out - ring_w             # + strand outer ring
        r_rev_out, r_rev_in = r_fwd_in - 2, r_fwd_in - 2 - ring_w
        r_gc_out = r_rev_in - 6
        r_gc_in = r_gc_out - ring_w * 0.9
        r_sk_out = r_gc_in - 3
        r_sk_in = r_sk_out - ring_w * 0.9
        self._hits.clear()

        # focus wedge (what the grid is showing)
        if self.focus and self.length:
            a0, a1 = self._angle(self.focus[0]), self._angle(max(self.focus[1], self.focus[0] + self.length * 0.002))
            path = self._arc_path(cx, cy, 0, r_fwd_out + 4, a0, a1)
            p.fillPath(path, QColor(255, 0, 0, 55))
            p.setPen(QPen(QColor("#e00000"), 2.5)); p.setBrush(Qt.NoBrush); p.drawPath(path)

        # backbone
        p.setPen(QPen(QColor("#555"), 1.5)); p.setBrush(Qt.NoBrush)
        rm = (r_fwd_in + r_rev_out) / 2
        p.drawEllipse(QRectF(cx - rm, cy - rm, 2 * rm, 2 * rm))

        # ticks + coordinate labels
        step = 10 ** int(math.floor(math.log10(max(1, self.length / 8))))
        for m in (1, 2, 5, 10):
            if self.length / (m * step) <= 12:
                step = int(m * step); break
        tickf = QFont(f); tickf.setPointSizeF(max(5.5, f.pointSizeF() * 0.85)); p.setFont(tickf)
        tfm = QFontMetrics(tickf)
        for t in range(0, self.length, max(1, step)):
            a = self._angle(t)
            p.setPen(QColor("#666"))
            p.drawLine(self._pt(cx, cy, rm - 4, a), self._pt(cx, cy, rm + 4, a))
            lab = f"{t // 1000} kb" if self.length >= 5000 else str(t)
            if t:
                pt = self._pt(cx, cy, rm - 6 - tfm.height() * 0.6, a)
                wlab = tfm.horizontalAdvance(lab) + 6
                p.setPen(QColor("#8a8a8a"))
                p.drawText(QRectF(pt.x() - wlab / 2, pt.y() - tfm.height() / 2, wlab, tfm.height()),
                           Qt.AlignCenter, lab)
        p.setFont(f)

        # genes
        genes = [g for gs in (self.ann.genes_by_seq.values() if self.ann else []) for g in gs
                 if not self.seqid or g.seqid == self.seqid or True]
        for g in genes:
            a0, a1 = self._angle(g.start), self._angle(max(g.end, g.start + 1))
            r_in, r_o = (r_fwd_in, r_fwd_out) if g.strand > 0 else (r_rev_in, r_rev_out)
            path = self._arc_path(cx, cy, r_in, r_o, a0, a1, arrow=1 if g.strand > 0 else -1)
            col = gene_color(g)
            p.setPen(QPen(col.darker(140), 0.8)); p.setBrush(col)
            p.drawPath(path)
            self._hits.append((path, g))
            if self.show_labels and (a1 - a0) * r_out > fm.height() * 0.6:
                mid = (a0 + a1) / 2
                pt = self._pt(cx, cy, r_out + 4, mid)
                deg = math.degrees(mid)
                p.save()
                p.translate(pt)
                if -90 <= deg <= 90:
                    p.rotate(deg)
                    p.setPen(QColor("#222")); p.drawText(QPointF(2, fm.height() / 3), g.name)
                else:
                    p.rotate(deg + 180)
                    p.setPen(QColor("#222"))
                    p.drawText(QPointF(-fm.horizontalAdvance(g.name) - 2, fm.height() / 3), g.name)
                p.restore()

        # GC content and skew rings
        if (self.show_gc or self.show_skew) and self.fetch_seq:
            bins = max(120, int(TWO_PI * r_gc_out / 2))
            self._ensure_gc(bins)
            if self._gc:
                if self.show_gc:
                    mean = sum(self._gc) / len(self._gc)
                    mid_r = (r_gc_in + r_gc_out) / 2
                    p.setPen(Qt.NoPen)
                    for i, v in enumerate(self._gc):
                        a0 = self._angle(i * self.length / len(self._gc))
                        a1 = self._angle((i + 1) * self.length / len(self._gc))
                        d = (v - mean) * (r_gc_out - r_gc_in) * 3
                        d = max(-(mid_r - r_gc_in), min(mid_r - r_gc_in, d))
                        r0, r1 = (mid_r, mid_r + d) if d >= 0 else (mid_r + d, mid_r)
                        p.setBrush(QColor("#2c3e50"))
                        p.drawPath(self._arc_path(cx, cy, r0, r1, a0, a1))
                    p.setPen(QPen(QColor("#999"), 0.7)); p.setBrush(Qt.NoBrush)
                    p.drawEllipse(QRectF(cx - mid_r, cy - mid_r, 2 * mid_r, 2 * mid_r))
                if self.show_skew:
                    mid_r = (r_sk_in + r_sk_out) / 2
                    p.setPen(Qt.NoPen)
                    for i, v in enumerate(self._skew):
                        a0 = self._angle(i * self.length / len(self._skew))
                        a1 = self._angle((i + 1) * self.length / len(self._skew))
                        d = v * (r_sk_out - r_sk_in)
                        d = max(-(mid_r - r_sk_in), min(mid_r - r_sk_in, d))
                        r0, r1 = (mid_r, mid_r + d) if d >= 0 else (mid_r + d, mid_r)
                        p.setBrush(QColor("#27ae60") if d >= 0 else QColor("#8e44ad"))
                        p.drawPath(self._arc_path(cx, cy, r0, r1, a0, a1))
                    p.setPen(QPen(QColor("#999"), 0.7)); p.setBrush(Qt.NoBrush)
                    p.drawEllipse(QRectF(cx - mid_r, cy - mid_r, 2 * mid_r, 2 * mid_r))

        # centre caption
        p.setPen(QColor("#222"))
        cf = QFont(self.font()); cf.setPointSizeF(max(8.0, R * 0.075)); cf.setBold(True); p.setFont(cf)
        p.drawText(QRectF(cx - R / 2, cy - R * 0.16, R, R * 0.16), Qt.AlignCenter, self.title)
        cf.setBold(False); cf.setPointSizeF(max(7.0, R * 0.055)); p.setFont(cf)
        p.drawText(QRectF(cx - R / 2, cy, R, R * 0.14), Qt.AlignCenter, f"{self.length:,} bp")
        if self.show_gc and self._gc:
            gcm = 100 * sum(self._gc) / len(self._gc)
            p.drawText(QRectF(cx - R / 2, cy + R * 0.13, R, R * 0.14), Qt.AlignCenter, f"GC {gcm:.1f}%")

    def paintEvent(self, ev):
        p = QPainter(self)
        p.fillRect(self.rect(), self.palette().color(QPalette.Base))
        self.render_to(p, self.width(), self.height())
        p.end()

    # ---------------------------------------------------------------- export
    def export_image(self, path: str, size: int = 2000):
        if path.lower().endswith(".svg"):
            from PySide6.QtSvg import QSvgGenerator
            gen = QSvgGenerator()
            gen.setFileName(path); gen.setSize(QSize(size, size))
            gen.setViewBox(QRectF(0, 0, size, size))
            gen.setTitle(self.title)
            p = QPainter(gen); self.render_to(p, size, size, for_export=True); p.end()
        else:
            from PySide6.QtGui import QImage
            img = QImage(size, size, QImage.Format_ARGB32)
            img.fill(Qt.white)
            p = QPainter(img); self.render_to(p, size, size, for_export=True); p.end()
            img.save(path)

    # ---------------------------------------------------------------- mouse
    def _gene_at(self, pos):
        for path, g in reversed(self._hits):
            if path.contains(pos):
                return g
        return None

    def mouseMoveEvent(self, e):
        g = self._gene_at(e.position())
        if g:
            QToolTip.showText(QCursor.pos(), f"<b>{g.name}</b> {g.start + 1:,}-{g.end:,} "
                              f"({'+' if g.strand > 0 else '-'}) {g.biotype}", self)
            self.setCursor(Qt.PointingHandCursor)
        else:
            QToolTip.hideText(); self.setCursor(Qt.ArrowCursor)

    def mousePressEvent(self, e):
        g = self._gene_at(e.position())
        if g:
            self.positionClicked.emit(g.start)
        else:
            self.positionClicked.emit(self._pos_from(self.width() / 2, self.height() / 2,
                                                     e.position().x(), e.position().y()))

    def mouseDoubleClickEvent(self, e):
        g = self._gene_at(e.position())
        if g:
            self.geneActivated.emit(g)
