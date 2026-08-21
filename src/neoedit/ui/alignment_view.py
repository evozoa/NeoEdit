"""Custom-painted alignment grid (BioEdit-style) on a QAbstractScrollArea."""
from __future__ import annotations

from PySide6.QtCore import Qt, QRect, QPoint, QPointF, Signal, QTimer
from PySide6.QtGui import (QPainter, QColor, QFont, QFontMetrics, QPen, QBrush, QKeyEvent, QPolygonF,
                           QMouseEvent, QWheelEvent, QPalette)
from PySide6.QtWidgets import QAbstractScrollArea, QApplication, QMenu, QInputDialog, QToolButton

from ..model.alignment import AlignmentModel, GAP_CHARS, Feature
from ..model import colors as C
from ..analysis import translate as T

GAPSET = set(GAP_CHARS)


class _ColMap:
    """Dict-backed stand-in for a reference string: supports len() and [c]."""
    __slots__ = ("d", "n")

    def __init__(self, d, n):
        self.d, self.n = d, n

    def __len__(self):
        return self.n

    def __getitem__(self, c):
        return self.d.get(c, "")


class AlignmentView(QAbstractScrollArea):
    cursorChanged = Signal(int, int)          # row, col
    selectionChanged = Signal()
    modeChanged = Signal(str)
    featureActivated = Signal(object)         # Feature
    contextMenuWanted = Signal(QPoint)        # viewport position

    MODES = ("slide", "edit", "grab")          # BioEdit: Select/Slide, Edit, Grab & Drag
    MODE_LABELS = {"slide": "Select / Slide", "edit": "Edit", "grab": "Grab & Drag"}
    # BioEdit-style right-click behaviours
    RIGHT_CLICK_ACTIONS = (
        ("ins_sel", "Insert gap in selected sequence"),
        ("del_sel", "Delete gap in selected sequence"),
        ("ins_other", "Insert gap in all unselected sequences"),
        ("del_other", "Delete gap in all unselected sequences"),
        ("menu", "Show context menu"),
    )

    def __init__(self, model: AlignmentModel, parent=None):
        super().__init__(parent)
        self.model = model
        self.model.add_listener(self._on_model)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setMouseTracking(True)
        self.viewport().setCursor(Qt.IBeamCursor)

        self.font_size = 11
        self.row_pad = 0          # extra vertical pixels per row (line spacing)
        self.col_pad = 0          # extra horizontal pixels per character
        self.font_family = "Courier New"
        self.text_weight = "regular"   # regular | semibold (synthetic) | bold
        self.crisp_text = True         # no anti-aliasing -> BioEdit/GDI-like crisp glyphs
        self._font = QFont(self.font_family)
        self._font.setStyleHint(QFont.Monospace)   # fall back to any monospace if missing
        self._font.setFixedPitch(True)
        self._apply_font()

        self.name_w = 140
        self.ruler_h = 22
        self.cur_row = 0
        self.cur_col = 0
        self.anchor = None            # (row, col) selection anchor or None
        self.sel_rows: set[int] = set()   # whole-row selection (name click)
        self.mode = "slide"
        self.edit_submode = "insert"       # insert | overwrite (Edit mode only)
        self.slide_downstream_default = False   # BioEdit toggle: sliding moves whole downstream sequence by default
        self.right_click_action = "ins_sel"

        self.scheme_name = None
        self.scheme = {}
        self.color_mode = "scheme"    # scheme | identity | none
        self.color_target = "text"    # "text" = BioEdit normal view (colored letters); "background" = inverse view
        self.identity_ref = "consensus"   # consensus | first
        self.identity_color = "#000000"   # BioEdit shades identities black (white letters)
        self.dots_for_identity = False
        self.show_reference = True      # pinned reference row above the grid (the row the map/gene models follow)
        self.show_translation = False
        self.trans_frame = 0            # 0-2 = +1..+3, 3-5 = -1..-3
        self.trans_table = 1
        self.trans_from_features = True # translate annotated ORFs/CDS in their own frame and code
        self.trans_fill = "regions"     # "regions" = only inside features, "all" = whole row in the chosen frame
        self.translation_provider = None  # callable(row, c0, c1) -> [(start, end, strand, table, name)]
        self.show_features = True
        self.shade_features = False     # tint residues under a feature (off: only a marker bar)
        self.feature_alpha = 45
        self.shade_threshold = 0.5
        self._consensus_cache = None
        self._drag_block = None       # (rows, start, end, last_col, kind)  kind: "crunch" | "downstream"
        self.feature_provider = None  # optional callable(row, c0, c1) -> list[Feature] (e.g. genome annotation)
        self._dragging_sel = False
        self._press_pos = None

        self._refresh_timer = QTimer(self)
        self._refresh_timer.setSingleShot(True)
        self._refresh_timer.setInterval(0)
        self._refresh_timer.timeout.connect(self._refresh)
        self._auto_scheme()
        self._update_scrollbars()
        self._make_corner_buttons()

    # ------------------------------------------------------------- basics
    def set_model(self, model: AlignmentModel):
        self.model.remove_listener(self._on_model)
        self.model = model
        self.model.add_listener(self._on_model)
        self.cur_row = self.cur_col = 0
        self.anchor = None
        self.sel_rows.clear()
        self._auto_scheme()
        self._refresh()

    def _auto_scheme(self):
        schemes = C.schemes_for(self.model.seq_type)
        if self.scheme_name not in schemes:
            self.scheme_name = next(iter(schemes))
        self.scheme = schemes[self.scheme_name]

    def set_scheme(self, name):
        self.scheme_name = name
        self._auto_scheme()
        self.viewport().update()

    def set_font_family(self, family: str):
        self.font_family = family
        self._font = QFont(family)
        self._font.setStyleHint(QFont.Monospace)
        self._font.setFixedPitch(True)
        self._apply_font()
        self._refresh()

    def set_text_style(self, weight=None, crisp=None):
        if weight is not None:
            self.text_weight = weight
        if crisp is not None:
            self.crisp_text = crisp
        self._apply_font()
        self._refresh()

    def _apply_font(self):
        self._font.setPointSize(self.font_size)
        self._font.setBold(self.text_weight == "bold")
        self._font.setStyleStrategy(QFont.NoAntialias if self.crisp_text else QFont.PreferAntialias)
        fm = QFontMetrics(self._font)
        self.cell_w = fm.horizontalAdvance("W") + self.col_pad
        self.cell_h = fm.height() + self.row_pad
        self._fm = fm
        # glyph offsets to centre a character in its cell
        self._tx = max(0, (self.cell_w - fm.horizontalAdvance("W")) // 2)
        self._ty = (self.cell_h + fm.ascent() - fm.descent()) // 2

    def set_spacing(self, row_pad=None, col_pad=None):
        if row_pad is not None:
            self.row_pad = max(0, min(40, row_pad))
        if col_pad is not None:
            self.col_pad = max(0, min(30, col_pad))
        self._apply_font()
        self._refresh()

    def _make_corner_buttons(self):
        """BioEdit-style spacing toggles in the top-left corner of the grid."""
        self._corner_btns = []
        specs = [("\u25b2", "Increase line spacing", lambda: self.set_spacing(row_pad=self.row_pad + 1)),
                 ("\u25bc", "Decrease line spacing", lambda: self.set_spacing(row_pad=self.row_pad - 1)),
                 ("\u25b6", "Increase character spacing", lambda: self.set_spacing(col_pad=self.col_pad + 1)),
                 ("\u25c0", "Decrease character spacing", lambda: self.set_spacing(col_pad=self.col_pad - 1))]
        for txt, tip, fn in specs:
            b = QToolButton(self.viewport())
            b.setText(txt); b.setToolTip(tip); b.setAutoRaise(True)
            b.setFocusPolicy(Qt.NoFocus)
            f = b.font(); f.setPointSize(6); b.setFont(f)
            b.clicked.connect(fn)
            b.show()
            self._corner_btns.append(b)
        self._place_corner_buttons()

    def _place_corner_buttons(self):
        if not getattr(self, "_corner_btns", None):
            return
        size = max(14, min(18, self.ruler_h - 2))
        x = 2
        for b in self._corner_btns:
            b.setFixedSize(size, size)
            b.move(x, 1)
            x += size + 1

    def zoom(self, delta: int):
        self.font_size = max(5, min(32, self.font_size + delta))
        self._apply_font()
        self._refresh()

    def _on_model(self, what):
        self._consensus_cache = None
        if self.cur_row >= self.model.nrows:
            self.cur_row = max(0, self.model.nrows - 1)
        self._refresh_timer.start()

    def _refresh(self):
        self._auto_scheme()
        self._update_name_width()
        self._update_scrollbars()
        self.viewport().update()

    def _update_name_width(self):
        w = max((self._fm.horizontalAdvance(r.name) for r in self.model.rows), default=60)
        self.name_w = max(90, min(260, w + 16))

    @property
    def row_h(self) -> int:
        return self.cell_h * (2 if (self.show_translation and self.model.is_nucleotide()) else 1)

    @property
    def header_h(self) -> int:
        return self.ruler_h + (self.cell_h if self.pinned_row() >= 0 else 0)

    def pinned_row(self) -> int:
        """Index of the reference row shown in the header strip, or -1 when the strip is hidden
        (toggled off, or fewer than two sequences — a lone row needs no copy of itself)."""
        m = self.model
        if not self.show_reference or m.nrows < 2:
            return -1
        r = m.ref_row
        return r if 0 <= r < m.nrows else 0

    def _update_scrollbars(self):
        vp = self.viewport()
        cols = self.model.width + 1
        rows = self.model.nrows
        self.horizontalScrollBar().setRange(0, max(0, cols - self._visible_cols()))
        self.horizontalScrollBar().setPageStep(max(1, self._visible_cols()))
        self.verticalScrollBar().setRange(0, max(0, rows - self._visible_rows()))
        self.verticalScrollBar().setPageStep(max(1, self._visible_rows()))

    def _visible_cols(self):
        return max(1, (self.viewport().width() - self.name_w) // self.cell_w)

    def _visible_rows(self):
        return max(1, (self.viewport().height() - self.header_h) // self.row_h)

    def resizeEvent(self, e):
        self._update_scrollbars()
        super().resizeEvent(e)

    def scrollContentsBy(self, dx, dy):
        self.viewport().update()

    # ------------------------------------------------------- coordinates
    def cell_at(self, pos: QPoint):
        """-> (row, col) for a viewport pixel; col may exceed width; row may be -1 (header)."""
        x = pos.x() - self.name_w
        y = pos.y() - self.header_h
        col = x // self.cell_w + self.horizontalScrollBar().value()
        row = y // self.row_h + self.verticalScrollBar().value()
        if y < 0:
            row = -1
        return int(row), int(max(0, col))

    def cell_rect(self, row, col) -> QRect:
        x = self.name_w + (col - self.horizontalScrollBar().value()) * self.cell_w
        y = self.header_h + (row - self.verticalScrollBar().value()) * self.row_h
        return QRect(x, y, self.cell_w, self.cell_h)

    def ensure_visible(self, row=None, col=None):
        row = self.cur_row if row is None else row
        col = self.cur_col if col is None else col
        hs, vs = self.horizontalScrollBar(), self.verticalScrollBar()
        if col < hs.value():
            hs.setValue(col)
        elif col >= hs.value() + self._visible_cols():
            hs.setValue(col - self._visible_cols() + 1)
        if row < vs.value():
            vs.setValue(row)
        elif row >= vs.value() + self._visible_rows():
            vs.setValue(row - self._visible_rows() + 1)

    # ------------------------------------------------------------ selection
    def selection(self):
        """-> (row0, row1, col0, col1) inclusive, or None."""
        if self.sel_rows:
            return (min(self.sel_rows), max(self.sel_rows), 0, max(0, self.model.width - 1))
        if self.anchor is None:
            return None
        r0, c0 = self.anchor
        r1, c1 = self.cur_row, self.cur_col
        return (min(r0, r1), max(r0, r1), min(c0, c1), max(c0, c1))

    def selected_rows(self) -> list[int]:
        if self.sel_rows:
            return sorted(self.sel_rows)
        s = self.selection()
        if s:
            return list(range(s[0], s[1] + 1))
        return [self.cur_row] if self.model.nrows else []

    def target_rows(self) -> list[int]:
        """Rows an operation should apply to: selection if any, else current row."""
        return self.selected_rows()

    def select_all(self):
        self.sel_rows = set(range(self.model.nrows))
        self.anchor = None
        self.selectionChanged.emit()
        self.viewport().update()

    def clear_selection(self):
        self.anchor = None
        self.sel_rows.clear()
        self.selectionChanged.emit()
        self.viewport().update()

    def select_region(self, r0, r1, c0, c1):
        self.sel_rows.clear()
        self.anchor = (r0, c0)
        self.cur_row, self.cur_col = r1, c1
        self.ensure_visible()
        self.selectionChanged.emit()
        self.cursorChanged.emit(self.cur_row, self.cur_col)
        self.viewport().update()

    def selected_text(self) -> str:
        s = self.selection()
        if not s:
            return ""
        r0, r1, c0, c1 = s
        out = []
        for i in range(r0, r1 + 1):
            seq = self.model.rows[i].seq
            out.append(f">{self.model.rows[i].name}\n{seq[c0:c1 + 1]}")
        return "\n".join(out) + "\n"

    def set_cursor(self, row, col, extend=False):
        row = max(0, min(self.model.nrows - 1, row)) if self.model.nrows else 0
        col = max(0, col)
        if extend:
            if self.anchor is None:
                self.anchor = (self.cur_row, self.cur_col)
        else:
            self.anchor = None
            self.sel_rows.clear()
        self.cur_row, self.cur_col = row, col
        self.ensure_visible()
        self.cursorChanged.emit(row, col)
        self.selectionChanged.emit()
        self.viewport().update()

    def set_mode(self, mode: str):
        self.mode = mode
        self.modeChanged.emit(mode)
        self.viewport().update()

    def set_edit_submode(self, sub: str):
        self.edit_submode = sub
        self.modeChanged.emit(self.mode)
        self.viewport().update()

    @property
    def typing(self) -> str | None:
        """'insert' / 'overwrite' when typing edits residues (Edit mode), else None."""
        return self.edit_submode if self.mode == "edit" else None

    def mode_text(self) -> str:
        t = self.MODE_LABELS[self.mode]
        if self.mode == "edit":
            t += f" ({self.edit_submode})"
        return t

    def _begin_drag(self, rows, start, end, col, shift):
        downstream = shift != self.slide_downstream_default   # shift reverses the default
        self._drag_block = (rows, start, end, col, "downstream" if downstream else "crunch")
        self.model.begin_batch("Move downstream" if downstream else "Slide block")
        self.viewport().setCursor(Qt.ClosedHandCursor)

    # ------------------------------------------------------------ consensus
    def consensus_at(self, c: int) -> str:
        """Consensus residue for one column, lazily cached (cache cleared on any edit)."""
        cache = self._consensus_cache
        if cache is None:
            cache = self._consensus_cache = {}
        ch = cache.get(c)
        if ch is None:
            counts = {}
            n = 0
            for r in self.model.rows:
                sq = r.seq
                if c < len(sq):
                    x = sq[c]
                    if x not in GAPSET:
                        x = x.upper()
                        counts[x] = counts.get(x, 0) + 1
                        n += 1
            if not n:
                ch = "-"
            else:
                best, bn = max(counts.items(), key=lambda kv: kv[1])
                frac = bn / n
                ch = best if (frac >= self.shade_threshold or frac >= 0.5) else ("N" if self.model.is_nucleotide() else "X")
            cache[c] = ch
        return ch

    def consensus(self) -> str:
        return "".join(self.consensus_at(c) for c in range(self.model.width))

    # ------------------------------------------------------------- painting
    def paintEvent(self, ev):
        p = QPainter(self.viewport())
        p.setFont(self._font)
        pal = self.palette()
        bg = pal.color(QPalette.Base)
        fg = pal.color(QPalette.Text)
        dark = bg.lightness() < 128
        p.fillRect(self.viewport().rect(), bg)
        m = self.model
        hs, vs = self.horizontalScrollBar().value(), self.verticalScrollBar().value()
        ncols = self._visible_cols() + 1
        nrows = self._visible_rows() + 1
        c0, c1 = hs, min(m.width, hs + ncols)
        r0, r1 = vs, min(m.nrows, vs + nrows)
        sel = self.selection()
        need_cons = (self.color_mode == "identity" and self.identity_ref == "consensus") or self.dots_for_identity
        if need_cons:
            cons_vis = {c: self.consensus_at(c) for c in range(c0, c1)}
        else:
            cons_vis = {}
        if self.identity_ref == "consensus":
            ref_seq = _ColMap(cons_vis, m.width)
        else:
            ref_seq = m.rows[0].seq if m.rows else ""
        grid_left = self.name_w
        W = self.viewport().width()
        H = self.viewport().height()

        # --- ruler
        p.setPen(fg)
        ruler_y = 0
        p.fillRect(0, 0, W, self.header_h, pal.color(QPalette.Window))
        for c in range(c0, c1 + 1):
            x = grid_left + (c - hs) * self.cell_w
            n = c + 1
            if n % 10 == 0:
                p.drawLine(x + self.cell_w // 2, self.ruler_h - 6, x + self.cell_w // 2, self.ruler_h - 1)
                tw = self._fm.horizontalAdvance(str(n))
                p.drawText(x + self.cell_w // 2 - tw // 2, self.ruler_h - 8, str(n))
            elif n % 5 == 0:
                p.drawLine(x + self.cell_w // 2, self.ruler_h - 4, x + self.cell_w // 2, self.ruler_h - 1)
        # --- pinned reference row (stays put while the grid scrolls; same colors as in the grid)
        pr = self.pinned_row()
        if pr >= 0:
            y = self.ruler_h
            row = m.rows[pr]
            p.setPen(Qt.NoPen); p.setBrush(QColor("#e00000"))
            p.drawPolygon(QPolygonF([QPointF(5, y + 3), QPointF(5, y + self.cell_h - 3), QPointF(10, y + self.cell_h / 2)]))
            p.setPen(fg)
            p.drawText(13, y + self._ty, self._fm.elidedText(row.name, Qt.ElideRight, grid_left - 17))
            seq = row.seq
            for c in range(c0, c1):
                ch = seq[c] if c < len(seq) else ""
                x = grid_left + (c - hs) * self.cell_w
                self._draw_cell(p, x, y, ch, self._cell_color(ch, c, ref_seq), fg, False,
                                force_bg=(self.color_mode == "identity"))
        p.setPen(QPen(pal.color(QPalette.Mid)))
        p.drawLine(0, self.header_h - 1, W, self.header_h - 1)
        p.drawLine(grid_left - 1, 0, grid_left - 1, H)

        # --- rows
        selcol = QColor(40, 90, 200, 110) if not dark else QColor(120, 170, 255, 110)
        for r in range(r0, r1):
            row = m.rows[r]
            y = self.header_h + (r - vs) * self.row_h
            # name
            row_selected = r in self.sel_rows or (sel and sel[0] <= r <= sel[1] and self.sel_rows)
            p.fillRect(0, y, grid_left - 1, self.row_h,
                       pal.color(QPalette.Highlight) if r in self.sel_rows else pal.color(QPalette.Window))
            name_x = 4
            if r == m.ref_row and m.nrows > 1:
                # pinned reference: small triangle marker
                p.setPen(Qt.NoPen); p.setBrush(QColor("#e00000"))
                p.drawPolygon(QPolygonF([QPointF(name_x + 1, y + 3), QPointF(name_x + 1, y + self.cell_h - 3),
                                         QPointF(name_x + 6, y + self.cell_h / 2)]))
                name_x += 9
            if row.group:
                gcol = QColor(m.group_color(row.group))
                p.fillRect(0, y, 5, self.row_h, gcol)
                name_x = 8
            p.setPen(pal.color(QPalette.HighlightedText) if r in self.sel_rows else fg)
            name = self._fm.elidedText(row.name, Qt.ElideRight, grid_left - name_x - 4)
            p.drawText(name_x, y + self._ty, name)
            # residues
            seq = row.seq
            for c in range(c0, c1):
                ch = seq[c] if c < len(seq) else ""
                x = grid_left + (c - hs) * self.cell_w
                bgc = self._cell_color(ch, c, ref_seq)
                disp = ch
                if self.dots_for_identity and ch and c < len(ref_seq) and r != (0 if self.identity_ref == "first" else -1) \
                        and ch.upper() == ref_seq[c].upper() and ch not in GAPSET:
                    disp = "."
                in_sel = bool(sel and sel[0] <= r <= sel[1] and sel[2] <= c <= sel[3])
                self._draw_cell(p, x, y, disp, bgc, fg, in_sel, selcol, force_bg=(self.color_mode == "identity"))
            # translation line beneath the residues
            if self.show_translation and m.is_nucleotide():
                self._draw_translation(p, r, seq, y, c0, c1, hs, grid_left, fg, dark)
        # --- features overlay
        feats = list(m.features) if m.features else []
        if self.show_features and self.feature_provider is not None:
            for r in range(r0, r1):
                try:
                    feats.extend(self.feature_provider(r, c0, c1))
                except Exception:
                    pass
        if self.show_features and feats:
            for f in feats:
                if not (r0 <= f.row < r1):
                    continue
                if f.end <= c0 or f.start >= c1:
                    continue
                y = self.header_h + (f.row - vs) * self.row_h
                xs = grid_left + (max(f.start, c0) - hs) * self.cell_w
                xe = grid_left + (min(f.end, c1) - hs) * self.cell_w
                col = QColor(f.color)
                if self.shade_features:
                    col.setAlpha(self.feature_alpha)
                    p.fillRect(xs, y, xe - xs, self.cell_h, col)
                # marker bar under the row: keeps residues legible even when a
                # feature (e.g. a gene model) spans the whole visible width
                bar_h = 2
                p.fillRect(xs, y + self.cell_h - bar_h, max(1, xe - xs), bar_h, QColor(f.color))
                p.setPen(QPen(QColor(f.color), 1))
                # arrow head
                if f.strand > 0:
                    p.drawLine(xe - 4, y + self.cell_h - 6, xe - 1, y + self.cell_h - 2)
                else:
                    p.drawLine(xs + 4, y + self.cell_h - 6, xs + 1, y + self.cell_h - 2)
        # --- cursor
        if m.nrows and r0 <= self.cur_row < r1 and c0 <= self.cur_col <= c1:
            rc = self.cell_rect(self.cur_row, self.cur_col)
            pen = QPen(QColor("#ff3030") if self.mode == "edit" else pal.color(QPalette.Highlight), 2)
            p.setPen(pen)
            if self.typing == "insert":
                p.drawLine(rc.left(), rc.top(), rc.left(), rc.bottom())
            else:
                p.drawRect(rc.adjusted(0, 0, -1, -1))
        p.end()

    AA_REGION_BG = "#e8e8e8"          # neutral shading for residues inside an annotated ORF
    AA_REGION_BG_DARK = "#3a3a3a"

    def _translation_regions(self, row, c0, c1):
        if not self.trans_from_features or self.translation_provider is None:
            return []
        try:
            return self.translation_provider(row, c0, c1)
        except Exception:
            return []

    def _draw_translation(self, p: QPainter, row, seq, y, c0, c1, hs, grid_left, fg, dark):
        """Amino acids under the nucleotides.

        Residues inside an annotated ORF/CDS are translated in that feature's own
        frame with its own genetic code and shaded a neutral grey; everywhere else
        the chosen frame/table is used (or nothing, in 'regions' mode)."""
        ay = y + self.cell_h
        muted = QColor("#8a8a8a") if not dark else QColor("#9a9a9a")
        regions = self._translation_regions(row, c0, c1)
        covered: set[int] = set()
        bg = QColor(self.AA_REGION_BG_DARK if dark else self.AA_REGION_BG)
        for (rs, re_, strand, table, name) in regions:
            x0 = grid_left + (max(rs, c0) - hs) * self.cell_w
            x1 = grid_left + (min(re_, c1) - hs) * self.cell_w
            if x1 > x0:
                p.fillRect(int(x0), ay, int(x1 - x0), self.cell_h, bg)
            aa = T.translate_region(seq, rs, re_, table, strand)
            for k, ch in enumerate(aa):
                # column of the middle base of this codon, in feature orientation
                cc = (rs + k * 3 + 1) if strand > 0 else (re_ - 1 - (k * 3 + 1))
                covered.add(cc)
                if not (c0 <= cc < c1):
                    continue
                x = grid_left + (cc - hs) * self.cell_w
                p.setPen(QColor("#c02020") if ch == "*" else (QColor("#222") if not dark else QColor("#ddd")))
                self._text(p, x + self._tx, ay + self._ty, ch)
        if self.trans_fill == "all" or not regions:
            aa = T.translate_aligned(seq, self.trans_table, self.trans_frame % 3) \
                if self.trans_frame < 3 else T.translate_aligned_reverse(seq, self.trans_table, self.trans_frame - 3)
            p.setPen(muted)
            for k, ch in enumerate(aa):
                cc = (self.trans_frame % 3) + k * 3 + 1 if self.trans_frame < 3 else len(seq) - 1 - ((self.trans_frame - 3) + k * 3 + 1)
                if cc in covered or not (c0 <= cc < c1):
                    continue
                x = grid_left + (cc - hs) * self.cell_w
                self._text(p, x + self._tx, ay + self._ty, ch)
        p.setPen(fg)

    def _cell_color(self, ch, c, ref_seq):
        if not ch or ch in GAPSET:
            return None
        if self.color_mode == "none":
            return None
        if self.color_mode == "identity":
            if c < len(ref_seq) and ch.upper() == ref_seq[c].upper():
                return self.identity_color
            return None
        return self.scheme.get(ch.upper())

    def _text(self, p: QPainter, x, y, ch):
        p.drawText(QPointF(x, y), ch)
        if self.text_weight == "semibold":
            p.drawText(QPointF(x + (1.0 if self.crisp_text else 0.6), y), ch)

    @staticmethod
    def _invert(c: QColor) -> QColor:
        return QColor(255 - c.red(), 255 - c.green(), 255 - c.blue())

    def _draw_cell(self, p: QPainter, x, y, ch, bgc, fg, in_sel, selcol=None, force_bg=False):
        """Paint one cell. Selected cells are drawn with inverted colors (BioEdit-style):
        a gap becomes a white dash on black; in inverse view G/C/A/T backgrounds become
        white/yellow/fuchsia/cyan with black letters."""
        base = self.palette().color(QPalette.Base)
        if bgc and (self.color_target == "background" or force_bg):
            bg, pen = QColor(bgc), QColor(C.text_for(bgc))
        elif bgc:
            bg, pen = base, QColor(bgc)
        else:
            bg, pen = base, QColor(fg)
        if in_sel:
            bg, pen = self._invert(bg), self._invert(pen)
        if in_sel or bg != base:
            p.fillRect(x, y, self.cell_w, self.cell_h, bg)
        p.setPen(pen)
        if ch:
            self._text(p, x + self._tx, y + self._ty, ch)

    # ------------------------------------------------------------ mouse
    def mousePressEvent(self, e: QMouseEvent):
        self.setFocus()
        pos = e.position().toPoint()
        if e.button() == Qt.RightButton:
            self._right_click(pos, e.modifiers())
            return
        row, col = self.cell_at(pos)
        self._press_pos = pos
        if row < 0 and pos.y() >= self.ruler_h and self.pinned_row() >= 0:
            # pinned reference strip: jump to the reference row (name) or to that cell (grid)
            pr = self.pinned_row()
            if pos.x() < self.name_w:
                self.sel_rows = {pr}; self.anchor = None; self.cur_row = pr
                self.ensure_visible(pr, None)
                self.cursorChanged.emit(self.cur_row, self.cur_col); self.selectionChanged.emit()
                self.viewport().update()
            else:
                self.select_region(pr, pr, col, col)
                self.ensure_visible(pr, col)
            return
        if pos.x() < self.name_w:
            # name panel: row selection
            if row < 0 or row >= self.model.nrows:
                return
            if e.modifiers() & Qt.ControlModifier:
                self.sel_rows ^= {row}
            elif e.modifiers() & Qt.ShiftModifier and self.sel_rows:
                a = min(self.sel_rows | {row}); b = max(self.sel_rows | {row})
                self.sel_rows = set(range(a, b + 1))
            else:
                self.sel_rows = {row}
            self.anchor = None
            self.cur_row = row
            self.cursorChanged.emit(self.cur_row, self.cur_col)
            self.selectionChanged.emit()
            self.viewport().update()
            return
        if row < 0:
            # ruler click: select column
            if self.model.nrows:
                self.select_region(0, self.model.nrows - 1, col, col)
            return
        if row >= self.model.nrows:
            return
        shift = bool(e.modifiers() & Qt.ShiftModifier)
        sel = self.selection()
        inside = sel and not self.sel_rows and sel[0] <= row <= sel[1] and sel[2] <= col <= sel[3]
        if self.mode == "slide" and inside:
            # drag the selected block (Shift: move whole downstream sequence instead of crunching gaps)
            self._begin_drag(list(range(sel[0], sel[1] + 1)), sel[2], sel[3] + 1, col, shift)
            return
        if self.mode == "grab":
            seq = self.model.rows[row].seq
            if col < len(seq) and seq[col] not in GAPSET:
                # grab a single residue and drag it (Shift: everything downstream of it)
                self.anchor = None; self.sel_rows.clear()
                self.cur_row, self.cur_col = row, col
                self._begin_drag([row], col, col + 1, col, shift)
                return
        # selection / cursor placement (extend with Shift only when not dragging blocks)
        self.set_cursor(row, col, extend=shift and self.mode != "grab")
        self._dragging_sel = True

    def _right_click(self, pos, mods):
        m = self.model
        row, col = self.cell_at(pos)
        act = self.right_click_action
        if act == "menu" or (mods & Qt.ShiftModifier) or pos.x() < self.name_w or row < 0 or row >= m.nrows:
            self.contextMenuWanted.emit(pos)
            return
        # "selected" = the clicked row, or the selected rows if the click falls inside a row selection
        sel = self.selection()
        if self.sel_rows and row in self.sel_rows:
            target = sorted(self.sel_rows)
        elif sel and not self.sel_rows and sel[0] <= row <= sel[1] and sel[0] != sel[1]:
            target = list(range(sel[0], sel[1] + 1))
        else:
            target = [row]
        others = [r for r in range(m.nrows) if r not in target]
        if act == "ins_sel":
            m.insert_gap_columns(col, 1, target)
        elif act == "del_sel":
            if not m.delete_gap_columns(col, 1, target):
                QApplication.beep()
        elif act == "ins_other":
            if others:
                m.insert_gap_columns(col, 1, others)
        elif act == "del_other":
            if others and not m.delete_gap_columns(col, 1, others):
                QApplication.beep()
        self.cur_row, self.cur_col = row, col
        self.cursorChanged.emit(row, col)
        self.viewport().update()

    def mouseMoveEvent(self, e: QMouseEvent):
        pos = e.position().toPoint()
        row, col = self.cell_at(pos)
        if self._drag_block is not None and e.buttons() & Qt.LeftButton:
            rows, start, end, last, kind = self._drag_block
            delta = col - last
            if delta:
                ok = (self.model.move_downstream(rows, start, delta) if kind == "downstream"
                      else self.model.block_shift(rows, start, end, delta))
                if ok:
                    start += delta; end += delta
                    if self.mode == "grab" and end - start == 1:
                        self.anchor = None
                        self.cur_row, self.cur_col = rows[0], start
                        self.cursorChanged.emit(self.cur_row, self.cur_col)
                    else:
                        self.anchor = (rows[0], start)
                        self.cur_row, self.cur_col = rows[-1], end - 1
                    self._drag_block = (rows, start, end, col, kind)
                    self.ensure_visible(self.cur_row, col)
                    self.selectionChanged.emit()
                    self.viewport().update()
            return
        if self._dragging_sel and e.buttons() & Qt.LeftButton:
            if row < 0:
                row = 0
            if pos.x() < self.name_w:
                col = self.horizontalScrollBar().value()
            row = min(max(0, row), self.model.nrows - 1)
            if self.anchor is None:
                self.anchor = (self.cur_row, self.cur_col)
            self.cur_row, self.cur_col = row, col
            self.ensure_visible()
            self.cursorChanged.emit(row, col)
            self.selectionChanged.emit()
            self.viewport().update()
            return
        if 0 <= row < self.model.nrows and pos.x() < self.name_w:
            g = self.model.rows[row].group
            d = self.model.rows[row].description
            tip = "\n".join(x for x in (self.model.rows[row].name,
                                        "pinned reference" if row == self.model.ref_row and self.model.nrows > 1 else "",
                                        f"group: {g}" if g else "", d) if x)
            self.setToolTip(tip)
            return
        # hover: feature tooltip
        if row >= 0 and row < self.model.nrows and pos.x() >= self.name_w:
            for f in self._features_at(row, col):
                if f.row == row and f.start <= col < f.end:
                    self.setToolTip(f"{f.type}: {f.label}\n{f.start + 1}-{f.end} ({'+' if f.strand > 0 else '-'})")
                    return
        self.setToolTip("")

    def mouseReleaseEvent(self, e: QMouseEvent):
        if self._drag_block is not None:
            self._drag_block = None
            self.model.end_batch()
            self.viewport().setCursor(Qt.IBeamCursor)
        self._dragging_sel = False
        # click without drag inside the selection clears it and places cursor
        if self._press_pos is not None and e.position().toPoint() == self._press_pos and not self.sel_rows:
            row, col = self.cell_at(self._press_pos)
            if 0 <= row < self.model.nrows and self._press_pos.x() >= self.name_w:
                if self.anchor is not None and not (e.modifiers() & Qt.ShiftModifier):
                    s = self.selection()
                    if s and (s[0] != s[1] or s[2] != s[3]):
                        self.set_cursor(row, col)
        self._press_pos = None

    def mouseDoubleClickEvent(self, e: QMouseEvent):
        pos = e.position().toPoint()
        row, col = self.cell_at(pos)
        if pos.x() < self.name_w and 0 <= row < self.model.nrows:
            name, ok = QInputDialog.getText(self, "Rename sequence", "Name:", text=self.model.rows[row].name)
            if ok and name:
                self.model.rename(row, name)
            return
        # double-click a feature selects it
        for f in self._features_at(row, col):
            if f.row == row and f.start <= col < f.end:
                self.select_region(row, row, f.start, f.end - 1)
                self.featureActivated.emit(f)
                return

    def _features_at(self, row, col):
        out = [f for f in self.model.features if f.row == row and f.start <= col < f.end]
        if self.feature_provider is not None:
            try:
                out += [f for f in self.feature_provider(row, col, col + 1) if f.start <= col < f.end]
            except Exception:
                pass
        return out

    def wheelEvent(self, e: QWheelEvent):
        if e.modifiers() & Qt.ControlModifier:
            self.zoom(1 if e.angleDelta().y() > 0 else -1)
            return
        if e.modifiers() & Qt.ShiftModifier:
            hs = self.horizontalScrollBar()
            hs.setValue(hs.value() - e.angleDelta().y() // 40)
            return
        super().wheelEvent(e)

    # ------------------------------------------------------------ keyboard
    def keyPressEvent(self, e: QKeyEvent):
        m = self.model
        if not m.nrows:
            return
        k = e.key()
        mods = e.modifiers()
        shift = bool(mods & Qt.ShiftModifier)
        ctrl = bool(mods & Qt.ControlModifier)
        r, c = self.cur_row, self.cur_col
        rows = self.target_rows()
        sel = self.selection()
        multi = sel is not None and (sel[0] != sel[1] or bool(self.sel_rows))

        if k == Qt.Key_Left:
            self.set_cursor(r, max(0, c - (10 if ctrl else 1)), shift); return
        if k == Qt.Key_Right:
            self.set_cursor(r, c + (10 if ctrl else 1), shift); return
        if k == Qt.Key_Up:
            self.set_cursor(r - 1, c, shift); return
        if k == Qt.Key_Down:
            self.set_cursor(r + 1, c, shift); return
        if k == Qt.Key_Home:
            self.set_cursor(0 if ctrl else r, 0, shift); return
        if k == Qt.Key_End:
            self.set_cursor(m.nrows - 1 if ctrl else r, max(0, len(m.rows[r].seq) - 1) if not ctrl else m.width - 1, shift); return
        if k == Qt.Key_PageDown:
            self.set_cursor(r + self._visible_rows(), c, shift); return
        if k == Qt.Key_PageUp:
            self.set_cursor(r - self._visible_rows(), c, shift); return
        if k == Qt.Key_Escape:
            self.clear_selection(); return

        if k in (Qt.Key_Space, Qt.Key_Minus) or (k == Qt.Key_Period and self.typing is None):
            # insert gap(s) at cursor; Ctrl -> whole column in all rows; in selection -> all selected rows
            if ctrl:
                m.insert_gap_columns(c, 1)
            elif multi:
                m.insert_gap_columns(sel[2] if not self.sel_rows else c, 1, rows)
                if not self.sel_rows:
                    # keep the block selected (it moved one column right)
                    self.anchor = (sel[0], sel[2] + 1)
                    self.cur_row, self.cur_col = sel[1], sel[3] + 1
                    self.selectionChanged.emit()
            else:
                m.insert_gaps(r, c, 1)
            if not multi:
                self.set_cursor(r, c + 1)
            else:
                self.viewport().update()
            return
        if k == Qt.Key_Delete:
            if ctrl:
                m.delete_gap_columns(c, 1)
            elif multi and not self.sel_rows:
                # delete selected gap columns in selected rows
                m.begin_batch("Delete gaps")
                for _ in range(sel[3] - sel[2] + 1):
                    if not m.delete_gap_columns(sel[2], 1, rows):
                        break
                m.end_batch()
                self.clear_selection()
            elif self.typing is None:
                if not m.delete_gaps(r, c, 1):
                    QApplication.beep()
            else:
                m.delete_range(r, c, 1)
            self.viewport().update()
            return
        if k == Qt.Key_Backspace:
            if c == 0:
                return
            if ctrl:
                if m.delete_gap_columns(c - 1, 1):
                    self.set_cursor(r, c - 1)
            elif multi:
                if m.delete_gap_columns(sel[2] - 1 if not self.sel_rows else c - 1, 1, rows):
                    self.anchor = (sel[0], sel[2] - 1)
                    self.set_cursor(sel[1], sel[3] - 1, extend=True)
            elif self.typing is None:
                if m.delete_gaps(r, c - 1, 1):
                    self.set_cursor(r, c - 1)
                else:
                    QApplication.beep()
            else:
                m.delete_range(r, c - 1, 1)
                self.set_cursor(r, c - 1)
            return
        txt = e.text()
        if k == Qt.Key_Insert and self.mode == "edit":
            self.set_edit_submode("overwrite" if self.edit_submode == "insert" else "insert")
            return
        if txt and txt.isprintable() and not ctrl and self.typing is not None:
            ch = txt.upper() if not shift else txt
            if self.typing == "overwrite":
                m.overwrite(r, c, ch)
            else:
                m.insert_text(r, c, ch)
            self.set_cursor(r, c + 1)
            return
        if txt and txt.isalpha() and self.typing is None and not ctrl:
            # typing in select mode: jump to next occurrence? keep BioEdit-like: beep
            QApplication.beep()
            return
        super().keyPressEvent(e)
