"""Genome-context panel: chromosome overview (tier 1) + region/gene-model view (tier 2).

Sits above the alignment grid (tier 3). Red boxes: tier-1 box = tier-2 window;
tier-2 box = columns currently visible in the grid.
"""
from __future__ import annotations

import math
from typing import Optional

from PySide6.QtCore import Qt, QRect, QRectF, QPointF, Signal, QSize
from PySide6.QtGui import QPainter, QColor, QPen, QFont, QFontMetrics, QPalette, QBrush, QPolygonF, QCursor
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QComboBox, QLineEdit, QToolButton, QLabel,
                               QToolTip, QCheckBox, QSizePolicy)

from ..genome.annotations import Annotation, Gene, Transcript, SyntenyBlock, pack_lanes, fmt_span
from ..model import colors as C

NT_COL = {"A": "#008000", "C": "#0000ff", "G": "#000000", "T": "#ff0000"}
RED = QColor("#e00000")


def fmt_bp(n: int) -> str:
    return f"{n:,}"


def nice_step(span: int, target_ticks: int = 8) -> int:
    raw = span / max(1, target_ticks)
    mag = 10 ** int(math.floor(math.log10(max(1, raw))))
    for m in (1, 2, 5, 10):
        if m * mag >= raw:
            return int(m * mag)
    return int(10 * mag)


def _hash_color(name: str) -> QColor:
    h = 0
    for ch in name:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFF
    return QColor.fromHsv(h % 360, 140, 200)


# ======================================================================== tier 1
class ChromosomeOverview(QWidget):
    """Whole-contig bar with a draggable red box = the tier-2 window."""
    windowChanged = Signal(int, int)     # new tier-2 window (start, end)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.length = 1
        self.win = (0, 1)
        self.focus = None                 # tier-3 visible region (drawn as a thin marker)
        self.synteny: list[SyntenyBlock] = []
        self.gene_density: list[int] = []
        self.setMinimumHeight(46)
        self.setMaximumHeight(46)
        self.setMouseTracking(True)
        self._drag = None

    def set_length(self, n: int):
        self.length = max(1, n); self.update()

    def set_window(self, s: int, e: int):
        self.win = (max(0, s), min(self.length, e)); self.update()

    def set_focus(self, s, e):
        self.focus = (s, e); self.update()

    def _bar(self) -> QRect:
        return QRect(8, 10, self.width() - 16, 14)

    def _x(self, pos: int) -> float:
        b = self._bar()
        return b.left() + b.width() * pos / self.length

    def _pos(self, x: float) -> int:
        b = self._bar()
        return int(max(0, min(self.length, (x - b.left()) / max(1, b.width()) * self.length)))

    def paintEvent(self, e):
        p = QPainter(self); p.setRenderHint(QPainter.Antialiasing)
        pal = self.palette()
        p.fillRect(self.rect(), pal.color(QPalette.Window))
        b = self._bar()
        p.setPen(QPen(QColor("#777"), 1)); p.setBrush(QColor("#e8e8e8")); p.drawRoundedRect(b, 4, 4)
        # synteny coloring on the bar
        for blk in self.synteny:
            x0, x1 = self._x(blk.qstart), self._x(blk.qend)
            if x1 - x0 < 0.5:
                continue
            col = _hash_color(blk.tname); col.setAlpha(200)
            p.fillRect(QRectF(x0, b.top() + 1, max(1.0, x1 - x0), b.height() - 2), col)
        # gene density (light ticks)
        if self.gene_density:
            mx = max(self.gene_density) or 1
            n = len(self.gene_density)
            for i, d in enumerate(self.gene_density):
                if d:
                    x = b.left() + b.width() * i / n
                    h = (b.height() - 4) * d / mx
                    p.fillRect(QRectF(x, b.bottom() - 1 - h, max(1.0, b.width() / n), h), QColor(0, 0, 0, 60))
        # Ticks. When the tier-2 window covers most of the contig the two rulers show
        # the same thing at the same scale, so only the contig extremes are labelled
        # here and the detailed numbers are left to the region view below.
        p.setPen(QColor("#555")); f = self.font(); f.setPointSize(7); p.setFont(f)
        fm = QFontMetrics(f)
        span = max(1, self.win[1] - self.win[0])
        redundant = span > 0.6 * self.length
        step = nice_step(self.length, 10)
        for t in range(0, self.length + 1, step):
            x = self._x(t)
            p.drawLine(QPointF(x, b.bottom() + 1), QPointF(x, b.bottom() + 4))
            if redundant and t != 0:
                continue
            lab = f"{t / 1e6:g} Mb" if self.length >= 2_000_000 else (f"{t / 1e3:g} kb" if self.length >= 2000 else str(t))
            w = fm.horizontalAdvance(lab)
            xx = min(max(b.left(), x - w / 2), b.right() - w)
            p.drawText(QPointF(xx, b.bottom() + 13), lab)
        if redundant:
            lab = f"{self.length / 1e6:g} Mb" if self.length >= 2_000_000 else f"{self.length / 1e3:g} kb"
            w = fm.horizontalAdvance(lab)
            p.drawText(QPointF(b.right() - w, b.bottom() + 13), lab)
        # focus marker (tier 3) and window box (tier 2)
        if self.focus:
            x0, x1 = self._x(self.focus[0]), self._x(self.focus[1])
            p.fillRect(QRectF(x0, b.top() - 4, max(2.0, x1 - x0), 3), QColor("#1f5fbf"))
        x0, x1 = self._x(self.win[0]), self._x(self.win[1])
        p.setPen(QPen(RED, 2)); p.setBrush(Qt.NoBrush)
        p.drawRect(QRectF(x0, b.top() - 2, max(4.0, x1 - x0), b.height() + 4))
        p.end()

    def mousePressEvent(self, e):
        x = e.position().x()
        x0, x1 = self._x(self.win[0]), self._x(self.win[1])
        if x0 - 4 <= x <= x1 + 4:
            self._drag = ("move", x, self.win)
        else:
            # jump: center window at click
            w = self.win[1] - self.win[0]
            c = self._pos(x)
            s = max(0, min(self.length - w, c - w // 2))
            self.windowChanged.emit(s, s + w)

    def mouseMoveEvent(self, e):
        if self._drag:
            _, x_start, (s0, e0) = self._drag
            d = self._pos(e.position().x()) - self._pos(x_start)
            w = e0 - s0
            s = max(0, min(self.length - w, s0 + d))
            self.windowChanged.emit(s, s + w)
        else:
            self.setCursor(Qt.SizeHorCursor if self._x(self.win[0]) - 4 <= e.position().x() <= self._x(self.win[1]) + 4 else Qt.PointingHandCursor)

    def mouseReleaseEvent(self, e):
        self._drag = None

    def wheelEvent(self, e):
        s, en = self.win
        c = self._pos(e.position().x())
        factor = 0.8 if e.angleDelta().y() > 0 else 1.25
        w = max(100, min(self.length, int((en - s) * factor)))
        ns = int(c - (c - s) * w / max(1, en - s))
        ns = max(0, min(self.length - w, ns))
        self.windowChanged.emit(ns, ns + w)


# ======================================================================== tier 2
MINOR_TYPES = {"misc_feature", "source", "repeat_region", "sts", "variation", "misc_difference",
               "unsure", "assembly_gap", "gap", "primer_bind", "protein_bind", "misc_binding"}


class RegionView(QWidget):
    """Gene models + synteny for the current window; red box = grid's visible columns."""
    windowChanged = Signal(int, int)
    focusRequested = Signal(int, int)       # user clicked: scroll grid to show (start,end)
    geneActivated = Signal(object)           # double-click gene
    insertionClicked = Signal(int, int)      # alignment columns of an insertion
    hoverInfo = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.length = 1
        self.win = (0, 1)
        self.focus = None
        self.ann: Optional[Annotation] = None
        self.seqid = ""
        self.synteny: list[SyntenyBlock] = []
        self.fetch_seq = None               # callable(start,end)->str for base-level drawing
        self.fetch_var = None               # callable(start,end)->list[float] per-ref-position variant freq (or None)
        self.insertions: list[tuple[int, int, int, int]] = []   # (ref_pos, n_cols, n_seqs, first_column)
        self.expanded = False               # show all transcripts
        self.show_minor = False             # misc_feature / repeat_region / variation etc.
        self.orf_tracks: list[dict] = []     # [{"label","color","items":[(s,e,strand,name,tip[,color])]}]
        self.circular = False                # circular molecule: features may continue past the origin (end > length)
        self.setMinimumHeight(150)
        self.setMouseTracking(True)
        self._drag = None
        self._gene_hits: list[tuple[QRectF, Gene, Transcript]] = []
        self._syn_hits: list[tuple[QRectF, SyntenyBlock]] = []
        self._orf_hits: list[tuple[QRectF, str, str]] = []
        self._ins_hits: list[tuple[QRectF, tuple]] = []
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    # --------------------------------------------------------------- state
    def set_window(self, s, e):
        s = max(0, s); e = min(self.length, max(s + 50, e))
        self.win = (s, e); self.update()

    def set_focus(self, s, e):
        self.focus = (s, e); self.update()

    def _x(self, pos) -> float:
        s, e = self.win
        return 40 + (self.width() - 48) * (pos - s) / max(1, e - s)

    def _pos(self, x) -> int:
        s, e = self.win
        return int(s + (x - 40) / max(1, self.width() - 48) * (e - s))

    def bp_per_px(self) -> float:
        return (self.win[1] - self.win[0]) / max(1, self.width() - 48)

    # --------------------------------------------------------------- paint
    def paintEvent(self, ev):
        p = QPainter(self); p.setRenderHint(QPainter.Antialiasing, False)
        pal = self.palette()
        p.fillRect(self.rect(), pal.color(QPalette.Base))
        W, H = self.width(), self.height()
        s, e = self.win
        left, right = 40, W - 8
        self._gene_hits.clear(); self._syn_hits.clear(); self._orf_hits.clear(); self._ins_hits.clear()
        f = self.font(); f.setPointSize(8); p.setFont(f); fm = QFontMetrics(f)

        # ruler
        y = 14
        p.setPen(QColor("#444"))
        p.drawLine(left, y, right, y)
        step = nice_step(e - s, 8)
        t0 = max(step, ((s + 1) // step) * step)
        for t in range(t0, e + 2, step):           # t = 1-based coordinate at a round number
            x = self._x(t - 1)
            if x < left - 1 or x > right + 1:
                continue
            p.drawLine(QPointF(x, y), QPointF(x, y - 5))
            lab = fmt_bp(t)
            p.drawText(QPointF(min(max(left, x - fm.horizontalAdvance(lab) / 2), right - fm.horizontalAdvance(lab)), 9), lab)
        y += 4

        # sequence track when zoomed in enough
        bpp = self.bp_per_px()
        if self.fetch_seq and bpp <= 1 / 6:
            seq = self.fetch_seq(s, e)
            px = 1 / bpp
            mono = QFont("Courier New", max(6, min(12, int(px * 0.9)))); mono.setBold(True)
            p.setFont(mono); fmm = QFontMetrics(mono)
            for i, ch in enumerate(seq):
                x = self._x(s + i)
                p.setPen(QColor(NT_COL.get(ch.upper(), "#888")))
                p.drawText(QPointF(x + (px - fmm.horizontalAdvance(ch)) / 2, y + 12), ch)
            p.setFont(f)
            y += 18
        elif self.fetch_seq and bpp <= 1:
            # GC-ish density strip? keep simple: thin line indicating sequence present
            p.fillRect(QRectF(left, y + 3, right - left, 3), QColor("#ddd")); y += 10

        # insertion carets: sequence the tracks have but the reference lacks
        if self.insertions:
            cy = y + 6
            for (rpos, ncols, nseq, col0) in self.insertions:
                if not (s <= rpos <= e):
                    continue
                x = self._x(rpos)
                w = 3 + min(4.0, ncols / 3.0)
                tri = QPolygonF([QPointF(x - w, cy + 6), QPointF(x + w, cy + 6), QPointF(x, cy - 2)])
                p.setPen(QPen(QColor("#b45309"), 1)); p.setBrush(QColor("#f59e0b"))
                p.drawPolygon(tri)
                self._ins_hits.append((QRectF(x - w - 2, cy - 3, 2 * w + 4, 11), (rpos, ncols, nseq, col0)))
            p.setPen(QColor("#666")); p.drawText(QPointF(2, cy + 6), "ins")
            y = cy + 10

        # variation track (population alignment in the grid)
        if self.fetch_var is not None:
            var = self.fetch_var(s, e)
            if var:
                base_y = y + 14
                p.setPen(QColor("#666")); p.drawText(QPointF(2, base_y - 2), "var")
                p.setPen(Qt.NoPen)
                n = len(var)
                for i, v in enumerate(var):
                    if v <= 0:
                        continue
                    x = self._x(s + i)
                    h = 3 + 9 * min(1.0, v)
                    col = QColor("#d97706") if v < 0.5 else QColor("#dc2626")
                    p.fillRect(QRectF(x, base_y - h, max(1.0, 1 / max(1e-9, self.bp_per_px())), h), col)
                p.setPen(QColor("#bbb")); p.drawLine(QPointF(left, base_y), QPointF(right, base_y))
                y = base_y + 2

        # gene lanes
        genes = self.ann.overlapping(self.seqid, s, e) if (self.ann and self.seqid) else []
        if not self.show_minor:
            genes = [g for g in genes if (g.biotype or "").lower() not in MINOR_TYPES]
        items = []
        L = self.length
        for g in genes:
            txs = g.transcripts if self.expanded else [max(g.transcripts, key=lambda t: (len(t.cds), t.end - t.start))] if g.transcripts else []
            for t in txs:
                label_w = fm.horizontalAdvance(g.name if not self.expanded else f"{g.name} {t.name}") + 6
                # a transcript past the origin is drawn twice: its tail at the right end and its
                # head (shifted by -length) at the left end of the molecule
                offsets = (0, -L) if t.end > L else (0,)
                for off in offsets:
                    a, b = t.start + off, t.end + off
                    if b <= s or a >= e:
                        continue
                    x0, x1 = self._x(max(a, s)), self._x(min(b, e))
                    items.append((int(min(x0, x1)), int(max(x1, x0 + label_w)), (g, t, off)))
        lanes = pack_lanes(items, gap=4)
        lane_h = 26 if not self.expanded else 22
        max_lanes = max(1, (H - y - 40) // lane_h)
        y0 = y + 4
        for li, lane in enumerate(lanes[:max_lanes]):
            yy = y0 + li * lane_h
            for _, _, (g, t, off) in lane:
                self._draw_transcript(p, fm, g, t, yy, lane_h, off)
        if len(lanes) > max_lanes:
            p.setPen(QColor("#a00")); p.drawText(QPointF(left, y0 + max_lanes * lane_h + 10), f"… {len(lanes) - max_lanes} more rows (zoom in or enlarge panel)")
        y = y0 + min(len(lanes), max_lanes) * lane_h + 6
        if not genes and self.ann is None:
            p.setPen(QColor("#888")); p.drawText(QPointF(left, y0 + 14), "No annotation loaded (Genome ▸ Load annotation…)")

        # ORF tracks (below the canonical gene models)
        if self.orf_tracks:
            p.setPen(QPen(QColor("#ccc"), 1))
            p.drawLine(left, y + 1, right, y + 1)
            y += 4
            for track in self.orf_tracks:
                col = QColor(track.get("color", "#7d3c98"))
                vis = []
                for it in track["items"]:
                    for off in ((0, -L) if it[1] > L else (0,)):
                        a, b = it[0] + off, it[1] + off
                        if b > s and a < e:
                            vis.append((a, b) + tuple(it[2:]))
                lanes = pack_lanes([(int(self._x(it[0])), int(self._x(it[1])) + fm.horizontalAdvance(it[3]) + 6, it)
                                    for it in vis], gap=4)
                p.setPen(QColor("#555"))
                p.drawText(QPointF(2, y + 9), track["label"][:6])
                lh = 13
                maxl = max(1, min(len(lanes), (H - y - (34 if self.synteny else 6)) // lh))
                for li, lane in enumerate(lanes[:maxl]):
                    yy = y + li * lh
                    for _a0, _a1, item in lane:
                        a, b, st, nm, tip = item[:5]
                        icol = QColor(item[5]) if len(item) > 5 and item[5] else col
                        x0, x1 = max(left, self._x(a)), min(right, self._x(b))
                        w = max(2.0, x1 - x0)
                        p.setPen(Qt.NoPen); p.setBrush(icol if st > 0 else icol.darker(130))
                        p.drawRect(QRectF(x0, yy + 1, w, 7))
                        # strand tick
                        p.setPen(QPen(icol.darker(160), 1))
                        if st > 0:
                            p.drawLine(QPointF(x1, yy + 1), QPointF(x1 + 3, yy + 4.5)); p.drawLine(QPointF(x1 + 3, yy + 4.5), QPointF(x1, yy + 8))
                        else:
                            p.drawLine(QPointF(x0, yy + 1), QPointF(x0 - 3, yy + 4.5)); p.drawLine(QPointF(x0 - 3, yy + 4.5), QPointF(x0, yy + 8))
                        if w > fm.horizontalAdvance(nm) * 0.6:
                            named = len(item) > 5 and item[5]
                            p.setPen(QColor(icol).darker(180) if named else QColor("#222"))
                            f2 = QFont(f); f2.setBold(bool(named)); p.setFont(f2)
                            p.drawText(QPointF(x0 + 2, yy + 8), nm)
                            p.setFont(f)
                        self._orf_hits.append((QRectF(x0, yy, w, 9), nm, tip))
                y += maxl * lh + 3
                if len(lanes) > maxl:
                    p.setPen(QColor("#a00")); p.drawText(QPointF(left, y + 8), f"… {len(lanes) - maxl} more ORF rows")
                    y += 10

        # synteny track (bottom)
        if self.synteny:
            ty = H - 22
            p.setPen(QColor("#666")); p.drawText(QPointF(2, ty + 10), "syn")
            for blk in self.synteny:
                if blk.qend <= s or blk.qstart >= e:
                    continue
                x0, x1 = max(left, self._x(blk.qstart)), min(right, self._x(blk.qend))
                if x1 - x0 < 1:
                    x1 = x0 + 1
                col = _hash_color(blk.tname)
                r = QRectF(x0, ty, x1 - x0, 12)
                p.fillRect(r, col if blk.strand > 0 else col.darker(140))
                self._syn_hits.append((r, blk))
                lab = f"{blk.tname}:{blk.tstart // 1000}k{'+' if blk.strand > 0 else '-'}"
                if x1 - x0 > fm.horizontalAdvance(lab) + 6:
                    p.setPen(QColor("#000")); p.drawText(QPointF(x0 + 3, ty + 10), lab)

        # origin of a circular molecule: marked at both ends of the linearised view
        if self.circular:
            p.setPen(QPen(QColor("#6d28d9"), 1, Qt.DashLine))
            for pos in (0, L):
                if s <= pos <= e:
                    x = self._x(pos)
                    p.drawLine(QPointF(x, 16), QPointF(x, H - 4))
                    p.setPen(QColor("#6d28d9"))
                    p.drawText(QPointF(min(x + 3, right - 36), H - 6), "origin")
                    p.setPen(QPen(QColor("#6d28d9"), 1, Qt.DashLine))
        # focus box (tier 3 region)
        if self.focus and self.focus[1] > s and self.focus[0] < e:
            fx0, fx1 = self._x(self.focus[0]), self._x(self.focus[1])
            p.setPen(QPen(RED, 2)); p.setBrush(Qt.NoBrush)
            p.drawRect(QRectF(max(left - 1, fx0), 2, max(3.0, min(right + 1, fx1) - max(left - 1, fx0)), H - 4))
        # side labels
        p.setPen(QColor("#666")); p.drawText(QPointF(2, 9), "bp")
        p.end()

    def _draw_transcript(self, p: QPainter, fm, g: Gene, t: Transcript, y: int, lane_h: int, off: int = 0):
        left, right = 40, self.width() - 8
        X = self._x
        self._x = lambda v, _X=X, _o=off: _X(v + _o)          # shift every coordinate of this piece
        try:
            self._draw_transcript_at(p, fm, g, t, y, lane_h, left, right, off)
        finally:
            self._x = X

    def _draw_transcript_at(self, p, fm, g, t, y, lane_h, left, right, off):
        x0, x1 = max(left, self._x(t.start)), min(right, self._x(t.end))
        if x1 < left or x0 > right:
            return
        mid = y + 8
        base = QColor("#1f5fbf") if g.strand > 0 else QColor("#c0392b")
        if g.biotype and "protein" not in g.biotype and g.biotype not in ("gene", "mRNA", "CDS", "bed"):
            base = QColor("#7d3c98") if g.strand > 0 else QColor("#a04000")
        if getattr(g, "cytoplasmic", False):
            base = QColor("#d926a9")        # cytoplasmically translated (MDP)
        if g.low_confidence:
            base.setAlpha(110)          # faded = partial / low-identity lift-over
        p.setPen(QPen(base, 1)); p.drawLine(QPointF(x0, mid), QPointF(x1, mid))
        # strand chevrons along intron line
        if x1 - x0 > 30:
            step = 18
            for cx in range(int(x0) + 9, int(x1) - 6, step):
                if g.strand > 0:
                    p.drawLine(QPointF(cx - 3, mid - 3), QPointF(cx, mid)); p.drawLine(QPointF(cx - 3, mid + 3), QPointF(cx, mid))
                else:
                    p.drawLine(QPointF(cx + 3, mid - 3), QPointF(cx, mid)); p.drawLine(QPointF(cx + 3, mid + 3), QPointF(cx, mid))
        p.setPen(Qt.NoPen)
        # UTR (thin) and exons/CDS (thick)
        cds = t.cds
        for es, ee in t.exons:
            ex0, ex1 = max(left, self._x(es)), min(right, self._x(ee))
            if ex1 < left or ex0 > right:
                continue
            w = max(1.0, ex1 - ex0)
            if cds:
                # draw exon as UTR-height first, then CDS part full height
                p.setBrush(base.lighter(150)); p.drawRect(QRectF(ex0, mid - 3, w, 6))
            else:
                p.setBrush(base); p.drawRect(QRectF(ex0, mid - 5, w, 10))
        for cs, ce in cds:
            cx0, cx1 = max(left, self._x(cs)), min(right, self._x(ce))
            if cx1 < left or cx0 > right:
                continue
            p.setBrush(base); p.drawRect(QRectF(cx0, mid - 6, max(1.0, cx1 - cx0), 12))
        # label
        label = g.name if not self.expanded else f"{g.name} {t.name}"
        p.setPen(QColor("#222") if not g.low_confidence else QColor("#888"))
        lx = max(left, self._x(t.start))
        p.drawText(QPointF(lx, mid + 17 if lane_h > 22 else mid + 15), label)
        hit = QRectF(x0, mid - 7, max(4.0, x1 - x0), 14)
        self._gene_hits.append((hit, g, t))

    # --------------------------------------------------------------- mouse
    def _hit(self, pos):
        for r, g, t in self._gene_hits:
            if r.contains(pos):
                return g, t
        return None

    def mousePressEvent(self, e):
        if e.button() == Qt.LeftButton:
            self._drag = (e.position().x(), self.win, False)

    def mouseMoveEvent(self, e):
        pos = e.position()
        if self._drag and e.buttons() & Qt.LeftButton:
            x0, (s0, e0), moved = self._drag
            dx = pos.x() - x0
            if abs(dx) > 2 or moved:
                d = int(-dx * self.bp_per_px())
                w = e0 - s0
                s = max(0, min(self.length - w, s0 + d))
                self._drag = (x0, (s0, e0), True)
                self.windowChanged.emit(s, s + w)
            return
        h = self._hit(pos)
        if h:
            g, t = h
            self.setCursor(Qt.PointingHandCursor)
            info = (f"<b>{g.name}</b> ({g.id})  {g.seqid}:{fmt_span(g.start, g.end, self.length if g.end > self.length else None)} "
                    f"({'+' if g.strand > 0 else '-'})  {g.biotype}<br>{t.name}: {len(t.exons)} exons"
                    + (f", CDS {sum(b - a for a, b in t.cds):,} bp" if t.cds else "")
                    + "".join(f"<br>{k}: {v}" for k, v in g.attrs.items() if k in ("description", "product"))
                    + (f"<br><i>lift-over coverage {g.attrs.get('coverage')}, identity {g.attrs.get('sequence_ID')}"
                       f"{' — partial' if g.attrs.get('partial_mapping') == 'True' else ''}"
                       f"{' — low identity' if g.attrs.get('low_identity') == 'True' else ''}</i>" if "coverage" in g.attrs else ""))
            QToolTip.showText(QCursor.pos(), info, self)
            self.hoverInfo.emit(f"{g.name}  {g.seqid}:{fmt_span(g.start, g.end, self.length if g.end > self.length else None)} ({'+' if g.strand > 0 else '-'})")
            return
        for r, ins in self._ins_hits:
            if r.contains(pos):
                rpos, ncols, nseq, col0 = ins
                QToolTip.showText(QCursor.pos(),
                                  f"<b>Insertion</b> after reference position {fmt_bp(rpos)}<br>"
                                  f"{ncols} column(s) absent from the reference, present in {nseq} sequence(s)<br>"
                                  f"<i>click to read the inserted bases in the alignment</i>", self)
                self.setCursor(Qt.PointingHandCursor)
                return
        for r, nm, tip in self._orf_hits:
            if r.contains(pos):
                QToolTip.showText(QCursor.pos(), tip, self)
                self.setCursor(Qt.PointingHandCursor)
                return
        for r, blk in self._syn_hits:
            if r.contains(pos):
                QToolTip.showText(QCursor.pos(), f"{blk.qname}:{fmt_bp(blk.qstart + 1)}-{fmt_bp(blk.qend)} {'+' if blk.strand > 0 else '-'} → "
                                  f"{blk.tname}:{fmt_bp(blk.tstart + 1)}-{fmt_bp(blk.tend)}  id {blk.identity:.1%}  mapq {blk.mapq}", self)
                return
        self.setCursor(Qt.OpenHandCursor)
        self.hoverInfo.emit(f"{self.seqid}:{fmt_bp(self._pos(pos.x()) + 1)}")

    def mouseReleaseEvent(self, e):
        if self._drag and not self._drag[2] and e.button() == Qt.LeftButton:
            # click (no drag): scroll the grid here
            for r, ins in self._ins_hits:
                if r.contains(e.position()):
                    self.insertionClicked.emit(ins[3], ins[3] + ins[1])
                    self._drag = None
                    return
            h = self._hit(e.position())
            if h:
                g, t = h
                self.focusRequested.emit(g.start, g.end)
            elif any(r.contains(e.position()) for r, _n, _t in self._orf_hits):
                for r, _n, _t in self._orf_hits:
                    if r.contains(e.position()):
                        a = self._pos(r.left()); b = self._pos(r.right())
                        self.focusRequested.emit(a, max(b, a + 1))
                        break
            else:
                c = self._pos(e.position().x())
                self.focusRequested.emit(c, c)
        self._drag = None

    def mouseDoubleClickEvent(self, e):
        h = self._hit(e.position())
        if h:
            self.geneActivated.emit(h[0])

    def wheelEvent(self, e):
        s, en = self.win
        c = self._pos(e.position().x())
        factor = 0.8 if e.angleDelta().y() > 0 else 1.25
        if e.modifiers() & Qt.ShiftModifier:
            d = int((en - s) * (0.1 if e.angleDelta().y() < 0 else -0.1))
            w = en - s
            ns = max(0, min(self.length - w, s + d))
            self.windowChanged.emit(ns, ns + w); return
        w = max(50, min(self.length, int((en - s) * factor)))
        ns = int(c - (c - s) * w / max(1, en - s))
        ns = max(0, min(self.length - w, ns))
        self.windowChanged.emit(ns, ns + w)


# ======================================================================== panel
class GenomePanel(QWidget):
    regionVisibilityChanged = Signal(bool)  # region (gene-model) view shown / hidden
    contigSelected = Signal(str)            # user picked a contig
    focusRequested = Signal(int, int)       # scroll grid to region
    geneActivated = Signal(object)
    insertionClicked = Signal(int, int)     # alignment columns (not reference coords)
    openRegionRequested = Signal(int, int)  # "open region in new editor window"

    def __init__(self, parent=None):
        super().__init__(parent)
        self.length = 1
        self.seqid = ""
        self.ann: Optional[Annotation] = None
        lay = QVBoxLayout(self); lay.setContentsMargins(2, 2, 2, 0); lay.setSpacing(2)
        bar = QHBoxLayout(); bar.setSpacing(4)
        self.contig_combo = QComboBox(); self.contig_combo.setMinimumWidth(160)
        self.contig_combo.currentIndexChanged.connect(self._contig_changed)
        self.region_edit = QLineEdit(); self.region_edit.setPlaceholderText("region, e.g. 1,200,000-1,250,000  or  gene name")
        self.region_edit.returnPressed.connect(self._goto_text)
        self.zoom_in = QToolButton(); self.zoom_in.setText("+"); self.zoom_in.setToolTip("Zoom in")
        self.zoom_out = QToolButton(); self.zoom_out.setText("−"); self.zoom_out.setToolTip("Zoom out")
        self.zoom_in.clicked.connect(lambda: self._zoom(0.5)); self.zoom_out.clicked.connect(lambda: self._zoom(2.0))
        self.expand_cb = QCheckBox("All transcripts")
        self.expand_cb.toggled.connect(self._toggle_expanded)
        self.minor_cb = QCheckBox("Minor features")
        self.minor_cb.setToolTip("Also show misc_feature, repeat_region, variation and similar "
                                 "non-genic annotations (e.g. the rCRS placeholder at 3107)")
        self.minor_cb.toggled.connect(self._toggle_minor)
        self.open_btn = QToolButton(); self.open_btn.setText("Open region in editor"); self.open_btn.setToolTip("Open the red-box region as a new alignment window with gene features")
        self.open_btn.clicked.connect(lambda: self.openRegionRequested.emit(*self.region.win))
        self.region_btn = QToolButton(); self.region_btn.setCheckable(True); self.region_btn.setChecked(True)
        self.region_btn.setText("Gene models ▾"); self.region_btn.setToolTip("Show / hide the region view (gene models, ORF tracks, synteny); the chromosome overview stays")
        self.region_btn.toggled.connect(self.set_region_visible)
        self.info = QLabel(""); self.info.setMinimumWidth(200)
        bar.addWidget(QLabel("Contig:")); bar.addWidget(self.contig_combo); bar.addWidget(self.region_edit, 1)
        bar.addWidget(self.zoom_in); bar.addWidget(self.zoom_out); bar.addWidget(self.expand_cb)
        bar.addWidget(self.minor_cb); bar.addWidget(self.open_btn); bar.addWidget(self.region_btn)
        bar.addWidget(self.info)
        lay.addLayout(bar)
        self.overview = ChromosomeOverview()
        self.region = RegionView()
        lay.addWidget(self.overview)
        lay.addWidget(self.region, 1)
        self.overview.windowChanged.connect(self.set_window)
        self.region.windowChanged.connect(self.set_window)
        self.region.focusRequested.connect(self.focusRequested)
        self.region.geneActivated.connect(self.geneActivated)
        self.region.hoverInfo.connect(self.info.setText)
        self.region.insertionClicked.connect(self.insertionClicked)

    # ------------------------------------------------------------ API
    def set_contigs(self, contigs: list[tuple[str, int]], current: str | None = None):
        self.contig_combo.blockSignals(True)
        self.contig_combo.clear()
        for n, L in contigs:
            self.contig_combo.addItem(f"{n}  ({L / 1e6:.2f} Mb)" if L >= 1e6 else f"{n}  ({L:,} bp)", n)
        if current:
            i = self.contig_combo.findData(current)
            if i >= 0:
                self.contig_combo.setCurrentIndex(i)
        self.contig_combo.blockSignals(False)

    def set_contig(self, seqid: str, length: int, fetch_seq=None):
        self.seqid = seqid; self.length = max(1, length)
        self.overview.set_length(self.length)
        self.region.length = self.length; self.region.seqid = seqid; self.region.fetch_seq = fetch_seq
        self.overview.synteny = []; self.region.synteny = []
        self._update_density()
        self.set_window(0, min(self.length, 200_000))

    def set_annotation(self, ann: Optional[Annotation]):
        self.ann = ann; self.region.ann = ann
        self._update_density()
        self.region.update()

    def set_insertions(self, ins: list[tuple[int, int, int, int]]):
        self.region.insertions = ins
        self.region.update()

    def set_circular(self, on: bool):
        self.region.circular = bool(on)
        self.region.update()

    def set_region_visible(self, on: bool):
        """Hide tier 2 (gene models / ORF tracks / synteny) and keep only the bar + chromosome overview."""
        on = bool(on)
        if self.region.isVisible() == on and self.region_btn.isChecked() == on:
            return
        self.region.setVisible(on)
        self.region_btn.blockSignals(True); self.region_btn.setChecked(on); self.region_btn.blockSignals(False)
        self.region_btn.setText("Gene models ▾" if on else "Gene models ▸")
        for w in (self.expand_cb, self.minor_cb, self.open_btn):
            w.setEnabled(on)
        self.regionVisibilityChanged.emit(on)

    def region_visible(self) -> bool:
        return self.region.isVisibleTo(self)

    def preferred_height(self) -> int:
        return 260 if self.region_visible() else self.overview.height() + 34

    def add_orf_track(self, label: str, items: list[tuple[int, int, int, str, str]], color: str = "#7d3c98",
                      replace: bool = True):
        """items: (start, end, strand, name, tooltip) in reference coordinates."""
        if replace:
            self.region.orf_tracks = [t for t in self.region.orf_tracks if t["label"] != label]
        self.region.orf_tracks.append({"label": label, "color": color, "items": items})
        self.region.update()

    def clear_orf_tracks(self):
        self.region.orf_tracks = []
        self.region.update()

    def orf_track_labels(self) -> list[str]:
        return [t["label"] for t in self.region.orf_tracks]

    def set_synteny(self, blocks: list[SyntenyBlock]):
        blocks = [b for b in blocks if b.qname == self.seqid]
        self.overview.synteny = blocks; self.region.synteny = blocks
        self.overview.update(); self.region.update()

    def set_focus(self, s: int, e: int):
        """Grid is showing columns [s,e) – draw red box; keep window containing it if it was."""
        self.overview.set_focus(s, e); self.region.set_focus(s, e)

    def set_window(self, s: int, e: int):
        s = max(0, s); e = min(self.length, max(s + 50, e))
        self.overview.set_window(s, e); self.region.set_window(s, e)
        self.region_edit.setText(f"{fmt_bp(s + 1)}-{fmt_bp(e)}")

    def window(self):
        return self.region.win

    def ensure_window_contains(self, s: int, e: int):
        ws, we = self.region.win
        if s >= ws and e <= we:
            return
        w = max(we - ws, e - s + 50)
        ns = max(0, min(self.length - w, (s + e) // 2 - w // 2))
        self.set_window(ns, ns + w)

    # ------------------------------------------------------------ internals
    def _update_density(self):
        n = 200
        dens = [0] * n
        if self.ann and self.seqid:
            for g in self.ann.genes_by_seq.get(self.seqid, []):
                i = min(n - 1, int(g.start / self.length * n)); dens[i] += 1
        self.overview.gene_density = dens if any(dens) else []
        self.overview.update()

    def _contig_changed(self, i):
        sid = self.contig_combo.itemData(i)
        if sid:
            self.contigSelected.emit(sid)

    def _zoom(self, factor):
        s, e = self.region.win
        c = (s + e) // 2
        w = max(50, min(self.length, int((e - s) * factor)))
        ns = max(0, min(self.length - w, c - w // 2))
        self.set_window(ns, ns + w)

    def _toggle_expanded(self, on):
        self.region.expanded = on; self.region.update()

    def _toggle_minor(self, on):
        self.region.show_minor = on; self.region.update()

    def _goto_text(self):
        txt = self.region_edit.text().strip().replace(",", "")
        import re
        m = re.match(r"^(?:[^:]+:)?\s*(\d+)\s*[-–:]\s*(\d+)$", txt)
        if m:
            s, e = int(m.group(1)) - 1, int(m.group(2))
            if e > s:
                self.set_window(s, e); self.focusRequested.emit(s, e)
            return
        if txt.isdigit():
            c = int(txt) - 1
            w = self.region.win[1] - self.region.win[0]
            self.set_window(c - w // 2, c - w // 2 + w); self.focusRequested.emit(c, c)
            return
        if self.ann:
            hits = [g for g in self.ann.find(txt) if g.seqid == self.seqid] or self.ann.find(txt)
            if hits:
                g = hits[0]
                if g.seqid != self.seqid:
                    self.info.setText(f"{g.name} is on {g.seqid}")
                    i = self.contig_combo.findData(g.seqid)
                    if i >= 0:
                        self.contig_combo.setCurrentIndex(i)
                    return
                pad = max(2000, len(g) // 2)
                self.set_window(g.start - pad, g.end + pad)
                self.focusRequested.emit(g.start, min(g.end, self.length))
                self.info.setText(f"{g.name}: {fmt_span(g.start, g.end, self.length if g.end > self.length else None)}")
                return
        self.info.setText("Not found")
