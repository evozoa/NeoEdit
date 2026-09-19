"""A basic tree viewer: rectangular phylogram / cladogram with support values and a scale bar;
click a branch to select its clade, reroot, ladderize, rotate, export SVG/PNG, and select the
clade's sequences in the alignment."""
from __future__ import annotations

import math
import os
import re

from PySide6.QtCore import Qt, QRectF, QPointF, QSize, Signal
from PySide6.QtGui import (QPainter, QColor, QPen, QFont, QFontMetricsF, QPalette, QAction, QKeySequence,
                           QGuiApplication)
from PySide6.QtWidgets import (QMainWindow, QWidget, QScrollArea, QToolBar, QLabel, QSpinBox, QFileDialog,
                               QMessageBox, QMenu, QToolTip, QSizePolicy)

from ..analysis import tree as T

MARGIN = 16
SCALE_BAR_H = 34
_ACCESSION_PREFIX = re.compile(r"^[A-Z]{1,6}$")


def spaced(label: str) -> str:
    """Tree label with "_" shown as spaces, except inside accessions: "NC_012920.1_Homo_sapiens"
    -> "NC_012920.1 Homo sapiens"."""
    parts = label.split("_")
    out = [parts[0]]
    for prev, part in zip(parts, parts[1:]):
        glue = "_" if _ACCESSION_PREFIX.match(prev) and part[:1].isdigit() else " "
        out.append(glue + part)
    return "".join(out)


def nice_length(target: float) -> float:
    """The 1, 2 or 5 x 10^k closest to `target` (on a log scale)."""
    if target <= 0:
        return 0.0
    k = math.floor(math.log10(target))
    cands = [m * 10.0 ** k for m in (1, 2, 5, 10)]
    return min(cands, key=lambda c: abs(math.log(c / target)))


class TreeCanvas(QWidget):
    nodeSelected = Signal(object)          # Node or None
    contextRequested = Signal(object, QPointF)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.root: T.Node | None = None
        self.names: dict[str, str] = {}
        self.full_names = False
        self.underscores_as_spaces = True
        self.show_support = True
        self.cladogram = False
        self.font_pt = 10
        self.row_factor = 1.0          # vertical zoom
        self.width_factor = 1.0        # horizontal zoom
        self.selected: T.Node | None = None
        self._layout: T.Layout | None = None
        self._hits: list[tuple[QRectF, T.Node]] = []
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    # ------------------------------------------------------------ data
    def set_tree(self, root: T.Node, names: dict[str, str] | None = None):
        self.root = root
        if names is not None:
            self.names = names
        self.selected = None
        self.relayout()

    def relayout(self):
        self._layout = T.layout(self.root, self.cladogram) if self.root else None
        self.updateGeometry()
        self.adjustSize()
        self.update()

    def tip_text(self, node: T.Node) -> str:
        s = node.name
        if self.full_names and s in self.names:
            return self.names[s]
        return spaced(s) if self.underscores_as_spaces else s

    def _font(self, scale: float = 1.0) -> QFont:
        f = QFont(self.font())
        f.setPointSizeF(self.font_pt * scale)
        return f

    def row_height(self) -> float:
        return QFontMetricsF(self._font()).height() * 1.25 * self.row_factor

    def label_width(self) -> float:
        if not self.root:
            return 0.0
        fm = QFontMetricsF(self._font())
        return max((fm.horizontalAdvance(self.tip_text(t)) for t in self.root.tips()), default=0.0)

    def sizeHint(self) -> QSize:
        if not self._layout:
            return QSize(400, 300)
        h = MARGIN * 2 + SCALE_BAR_H + max(1, self._layout.ntips - 1) * self.row_height() + self.row_height()
        vw = self.parent().width() if self.parent() else 600
        w = max(vw, 300) * self.width_factor
        return QSize(int(w), int(h))

    def minimumSizeHint(self) -> QSize:
        return self.sizeHint() if self._layout else QSize(200, 150)

    # ------------------------------------------------------------ drawing
    def paintEvent(self, ev):
        p = QPainter(self)
        p.fillRect(self.rect(), self.palette().color(QPalette.Base))
        self.render_to(p, self.width(), self.height())
        p.end()

    def render_to(self, p: QPainter, w: float, h: float, for_export: bool = False):
        self._hits = []
        L = self._layout
        if not L or not self.root:
            return
        p.setRenderHint(QPainter.Antialiasing, True)
        pal = self.palette()
        ink = QColor("black") if for_export else pal.color(QPalette.Text)
        faint = QColor("#666666") if for_export else pal.color(QPalette.PlaceholderText)
        hi = QColor("#1f6fd1") if for_export else pal.color(QPalette.Highlight)
        font, small = self._font(), self._font(0.8)
        fm, fms = QFontMetricsF(font), QFontMetricsF(small)
        labelw = self.label_width()
        left, top = MARGIN, MARGIN + fm.height() / 2
        right = labelw + MARGIN + 8
        plot_w = max(40.0, w - left - right)
        xs = plot_w / L.width if L.width > 0 else 0.0
        rows = max(1, L.ntips - 1)
        row_h = max(1.0, (h - top - MARGIN - SCALE_BAR_H - fm.height() / 2) / rows)

        def X(n):
            return left + L.x[id(n)] * xs

        def Y(n):
            return top + L.y[id(n)] * row_h

        sel = set(id(n) for n in self.selected.walk()) if self.selected is not None else set()
        pen = QPen(ink, 1.3)
        pen_sel = QPen(hi, 2.4)
        # branches
        for n in self.root.walk():
            x, y = X(n), Y(n)
            p.setPen(pen_sel if id(n) in sel else pen)
            if n.parent is not None:
                px = X(n.parent)
                p.drawLine(QPointF(px, y), QPointF(x, y))
                self._hits.append((QRectF(px, y - 5, max(x - px, 6), 10), n))
            if n.children:
                ys = [Y(c) for c in n.children]
                p.setPen(pen_sel if id(n) in sel else pen)
                p.drawLine(QPointF(x, min(ys)), QPointF(x, max(ys)))
        # tip labels
        p.setFont(font)
        for t in self.root.tips():
            text = self.tip_text(t)
            x, y = X(t) + 5, Y(t)
            rect = QRectF(x, y - fm.height() / 2, fm.horizontalAdvance(text) + 2, fm.height())
            if id(t) in sel:
                p.fillRect(rect, QColor(hi.red(), hi.green(), hi.blue(), 50))
            p.setPen(hi if id(t) in sel else ink)
            p.drawText(QPointF(x, y + fm.ascent() - fm.height() / 2), text)
            self._hits.append((rect, t))
        # support values, just above the branch, ending at the node
        if self.show_support:
            p.setFont(small)
            p.setPen(faint)
            for n in self.root.walk():
                if n.children and n.parent is not None and n.label:
                    tw = fms.horizontalAdvance(n.label)
                    p.drawText(QPointF(X(n) - tw - 2, Y(n) - 3), n.label)
        # scale bar
        if L.has_lengths and xs > 0:
            bar = nice_length(plot_w / xs / 5)
            y = h - MARGIN - SCALE_BAR_H / 2
            p.setPen(QPen(ink, 1.3))
            p.drawLine(QPointF(left, y), QPointF(left + bar * xs, y))
            for bx in (left, left + bar * xs):
                p.drawLine(QPointF(bx, y - 3), QPointF(bx, y + 3))
            p.setFont(small)
            p.drawText(QPointF(left, y + fms.ascent() + 4), f"{bar:g}")

    # ------------------------------------------------------------ export
    def export_image(self, path: str, scale: float = 2.0):
        w, h = max(self.width(), 600), max(self.height(), 200)
        if path.lower().endswith(".svg"):
            from PySide6.QtSvg import QSvgGenerator
            gen = QSvgGenerator()
            gen.setFileName(path)
            gen.setSize(QSize(w, h))
            gen.setViewBox(QRectF(0, 0, w, h))
            gen.setTitle(os.path.basename(path))
            p = QPainter(gen)
            self.render_to(p, w, h, for_export=True)
            p.end()
        else:
            from PySide6.QtGui import QImage
            img = QImage(int(w * scale), int(h * scale), QImage.Format_ARGB32)
            img.fill(Qt.white)
            p = QPainter(img)
            p.scale(scale, scale)
            self.render_to(p, w, h, for_export=True)
            p.end()
            if not img.save(path):
                raise OSError(f"could not write {path}")

    # ------------------------------------------------------------ mouse
    def node_at(self, pos: QPointF) -> T.Node | None:
        for rect, n in reversed(self._hits):
            if rect.contains(pos):
                return n
        return None

    def mousePressEvent(self, e):
        n = self.node_at(e.position())
        if e.button() == Qt.RightButton:
            if n is not None:
                self.select(n)
            self.contextRequested.emit(n, e.globalPosition())
            return
        self.select(n)

    def select(self, n):
        self.selected = n
        self.nodeSelected.emit(n)
        self.update()

    def mouseMoveEvent(self, e):
        n = self.node_at(e.position())
        if n is None:
            QToolTip.hideText()
            self.setCursor(Qt.ArrowCursor)
            return
        self.setCursor(Qt.PointingHandCursor)
        bits = []
        if n.is_tip:
            bits.append(f"<b>{self.names.get(n.name, n.name)}</b>")
        else:
            bits.append(f"clade of {len(n.tips())} sequences")
            if n.label:
                bits.append(f"support {n.label}")
        if n.length is not None:
            bits.append(f"branch length {n.length:.6g}")
        QToolTip.showText(e.globalPosition().toPoint(), "<br>".join(bits), self)

    def wheelEvent(self, e):
        if e.modifiers() & Qt.ControlModifier:
            self.zoom(1.25 if e.angleDelta().y() > 0 else 0.8)
            e.accept()
            return
        super().wheelEvent(e)

    def zoom(self, factor: float):
        self.row_factor = min(20.0, max(0.2, self.row_factor * factor))
        self.width_factor = min(20.0, max(1.0, self.width_factor * factor))
        self.relayout()

    def fit(self):
        self.row_factor = self.width_factor = 1.0
        self.relayout()


class TreeWindow(QMainWindow):
    """Top-level tree viewer. `select_names(list of tip labels)` selects them in the alignment."""

    def __init__(self, parent=None, root: T.Node | None = None, title: str = "", path: str = "",
                 names: dict[str, str] | None = None, select_names=None):
        super().__init__(parent)
        self.setWindowFlag(Qt.Window, True)
        self.setAttribute(Qt.WA_DeleteOnClose, False)
        self.resize(900, 700)
        self.select_names = select_names
        self.path = path
        self.original: T.Node | None = None
        self.canvas = TreeCanvas()
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.canvas)
        self.scroll.setBackgroundRole(QPalette.Base)
        self.setCentralWidget(self.scroll)
        self.canvas.nodeSelected.connect(self._selection_changed)
        self.canvas.contextRequested.connect(self._context_menu)
        self._build_toolbar()
        if root is not None:
            self.set_tree(root, title or os.path.basename(path), names)

    # ------------------------------------------------------------ setup
    def _act(self, text, slot, shortcut=None, checkable=False, tip=None):
        a = QAction(text, self)
        a.triggered.connect(slot)
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        a.setCheckable(checkable)
        if tip:
            a.setToolTip(tip); a.setStatusTip(tip)
        return a

    def _build_toolbar(self):
        tb = QToolBar("Tree")
        tb.setMovable(False)
        tb.setToolButtonStyle(Qt.ToolButtonTextOnly)
        self.addToolBar(tb)
        self.a_open = self._act("Open…", self.open_dialog, "Ctrl+O", tip="Open a Newick or NEXUS tree")
        self.a_save = self._act("Save tree…", self.save_tree, "Ctrl+S", tip="Save the tree as shown (rooting, order) as Newick")
        self.a_svg = self._act("Export SVG…", lambda: self.export("svg"))
        self.a_png = self._act("Export PNG…", lambda: self.export("png"))
        self.a_zin = self._act("Zoom in", lambda: self.canvas.zoom(1.25), "Ctrl++")
        self.a_zout = self._act("Zoom out", lambda: self.canvas.zoom(0.8), "Ctrl+-")
        self.a_fit = self._act("Fit", self.canvas.fit, "Ctrl+0", tip="Fit the tree to the window")
        self.a_clado = self._act("Cladogram", self._toggle("cladogram", relayout=True), checkable=True,
                                 tip="Ignore branch lengths and line the tips up")
        self.a_support = self._act("Support values", self._toggle("show_support"), checkable=True)
        self.a_support.setChecked(True)
        self.a_full = self._act("Full names", self._toggle("full_names", relayout=True), checkable=True,
                                tip="Show the full sequence names saved with the run (names.tsv) instead of the tree labels")
        self.a_mid = self._act("Midpoint root", self.midpoint_root,
                               tip="Root the tree halfway along its longest tip-to-tip path")
        self.a_root = self._act("Root here", self.root_selected, tip="Root the tree on the selected branch")
        self.a_rotate = self._act("Rotate", self.rotate_selected, tip="Swap the order of the selected node's branches")
        self.a_ladder = self._act("Ladderize", self.ladderize, tip="Order branches by clade size")
        self.a_reset = self._act("Original", self.reset, tip="Go back to the tree as it was loaded")
        self.a_select = self._act("Select in alignment", self.select_in_alignment,
                                  tip="Select the sequences of the selected clade in the alignment")
        self.a_close = self._act("Close", self.close, "Ctrl+W")
        mb = self.menuBar()
        for title, acts in (("&File", (self.a_open, self.a_save, None, self.a_svg, self.a_png, None, self.a_close)),
                            ("&View", (self.a_zin, self.a_zout, self.a_fit, None, self.a_clado, self.a_support,
                                       self.a_full)),
                            ("&Tree", (self.a_mid, self.a_root, self.a_rotate, self.a_ladder, self.a_reset, None,
                                       self.a_select))):
            m = mb.addMenu(title)
            for a in acts:
                m.addAction(a) if a else m.addSeparator()
        for a in (self.a_zin, self.a_zout, self.a_fit, None, self.a_mid, self.a_root, self.a_rotate,
                  self.a_ladder, None, self.a_select):
            tb.addAction(a) if a else tb.addSeparator()
        tb.addSeparator()
        tb.addWidget(QLabel(" Font "))
        self.font_size = QSpinBox()
        self.font_size.setRange(6, 36)
        self.font_size.setValue(self.canvas.font_pt)
        self.font_size.valueChanged.connect(self._font_changed)
        tb.addWidget(self.font_size)

    def _toggle(self, attr, relayout=False):
        def fn(on):
            setattr(self.canvas, attr, on)
            self.canvas.relayout() if relayout else self.canvas.update()
        return fn

    def _font_changed(self, v):
        self.canvas.font_pt = v
        self.canvas.relayout()

    # ------------------------------------------------------------ tree
    def set_tree(self, root: T.Node, title: str = "", names: dict[str, str] | None = None):
        self.original = T.parse_newick(T.to_newick(root))
        names = names or {}
        self.canvas.set_tree(root, names)
        # full names (from the run's names.tsv) are exact, so they are the default when known
        self.a_full.setEnabled(bool(names))
        self.a_full.setChecked(bool(names))
        self.canvas.full_names = bool(names)
        self.canvas.relayout()
        self.setWindowTitle(f"Tree — {title}" if title else "Tree")
        self._selection_changed(None)

    def load(self, path: str):
        root = T.read_tree(path)
        nf = T.names_file_for(path)
        names = T.read_names(nf) if nf else {}
        self.path = path
        self.set_tree(root, os.path.basename(path), names)

    def open_dialog(self):
        d = os.path.dirname(self.path) if self.path else ""
        path, _ = QFileDialog.getOpenFileName(self, "Open tree", d, TREE_FILTER)
        if path:
            try:
                self.load(path)
            except (T.TreeError, OSError) as e:
                QMessageBox.warning(self, "Open tree", f"{path}\n\n{e}")

    def _replace(self, root: T.Node):
        self.canvas.set_tree(root)
        self._selection_changed(None)

    def midpoint_root(self):
        self._replace(T.midpoint_root(self.canvas.root))
        self.statusBar().showMessage("Rooted at the midpoint", 4000)

    def root_selected(self):
        n = self.canvas.selected
        if n is None or n.parent is None:
            return
        self._replace(T.reroot(self.canvas.root, n))
        self.statusBar().showMessage("Rerooted on the selected branch", 4000)

    def rotate_selected(self):
        n = self.canvas.selected
        if n is not None and n.children:
            T.rotate(n)
            self.canvas.relayout()

    def ladderize(self):
        T.ladderize(self.canvas.root)
        self.canvas.relayout()

    def reset(self):
        if self.original is not None:
            self._replace(T.parse_newick(T.to_newick(self.original)))

    def _selection_changed(self, n):
        has = n is not None
        self.a_root.setEnabled(has and n.parent is not None)
        self.a_rotate.setEnabled(has and bool(n.children))
        self.a_select.setEnabled(has and self.select_names is not None)
        root = self.canvas.root
        if root is None:
            self.statusBar().showMessage("")
        elif n is None:
            self.statusBar().showMessage(f"{len(root.tips())} sequences.  Click a branch to select its clade; "
                                         "right-click for more.")
        elif n.is_tip:
            self.statusBar().showMessage(self.canvas.names.get(n.name, n.name))
        else:
            sup = f", support {n.label}" if n.label else ""
            self.statusBar().showMessage(f"Clade of {len(n.tips())} sequences{sup}")

    def _context_menu(self, n, gpos):
        m = QMenu(self)
        if n is not None:
            m.addAction(self.a_root)
            m.addAction(self.a_rotate)
            m.addAction(self.a_select)
            copy = m.addAction("Copy names")
            copy.triggered.connect(lambda: QGuiApplication.clipboard().setText(
                "\n".join(self.canvas.names.get(t.name, t.name) for t in n.tips())))
            m.addSeparator()
        m.addAction(self.a_mid)
        m.addAction(self.a_ladder)
        m.addAction(self.a_reset)
        m.popup(gpos.toPoint())

    def select_in_alignment(self):
        n = self.canvas.selected
        if n is None or self.select_names is None:
            return
        labels = [t.name for t in n.tips()]
        found = self.select_names(labels, self.canvas.names)
        self.statusBar().showMessage(f"Selected {found} of {len(labels)} sequences in the alignment", 5000)

    # ------------------------------------------------------------ files
    def save_tree(self):
        base = os.path.splitext(self.path)[0] + ".rerooted.nwk" if self.path else "tree.nwk"
        path, _ = QFileDialog.getSaveFileName(self, "Save tree", base, "Newick (*.nwk *.newick *.tre)")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(T.to_newick(self.canvas.root) + "\n")
        except OSError as e:
            QMessageBox.critical(self, "Save tree", str(e)); return
        self.statusBar().showMessage(f"Saved {path}", 5000)

    def export(self, kind: str):
        base = os.path.splitext(self.path)[0] if self.path else "tree"
        filt = "SVG (*.svg)" if kind == "svg" else "PNG (*.png)"
        path, _ = QFileDialog.getSaveFileName(self, "Export tree", f"{base}.{kind}", filt)
        if not path:
            return
        if not path.lower().endswith("." + kind):
            path += "." + kind
        try:
            self.canvas.export_image(path)
        except Exception as e:           # noqa: BLE001
            QMessageBox.critical(self, "Export failed", str(e)); return
        self.statusBar().showMessage(f"Saved {path}", 5000)

    def resizeEvent(self, e):
        super().resizeEvent(e)
        self.canvas.relayout()


TREE_EXTENSIONS = (".nwk", ".newick", ".tre", ".tree", ".treefile", ".contree")
TREE_FILTER = ("Trees (*.nwk *.newick *.tre *.tree *.treefile *.contree *.nex *.nexus *.txt);;"
               "All files (*)")
