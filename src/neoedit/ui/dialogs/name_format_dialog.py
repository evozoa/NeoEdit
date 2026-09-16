"""How sequences imported from NCBI are named: drag fields (accession, genus, species, isolate,
...) from "Available fields" into "Use in name"; the name is built from them in that order."""
from __future__ import annotations

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QListWidget,
                               QListWidgetItem, QAbstractItemView, QPushButton, QComboBox, QCheckBox,
                               QDialogButtonBox, QGroupBox)

from ...remote import ncbi as N

KEY_ROLE = Qt.UserRole
SETTING_FIELDS = "remote/ncbi_name_fields"
SETTING_SEP = "remote/ncbi_name_sep"
SETTING_SPACES = "remote/ncbi_name_spaces"

# Shown in the preview when no search results are at hand (real NCBI records)
EXAMPLES = [
    N.Summary("1", "MT455673.1", "Pimephales promelas voucher USNM:FISH:429789 cytochrome oxidase subunit 1 "
              "(COI) gene, partial cds; mitochondrial", 655, "Pimephales promelas", taxid="90988",
              genome="mitochondrion",
              source={"isolate": "SERCFISH0714", "specimen_voucher": "USNM:FISH:429789",
                      "country": "USA: Maryland, Anne Arundel County, Patuxent River, Galloway Branch",
                      "lat_lon": "38.804 N 76.693 W", "collection_date": "20-May-2013",
                      "collected_by": "Robert Aguilar"}),
    N.Summary("2", "NC_073531.1", "Rhynchocyon stuhlmanni voucher RMCA26091 mitochondrion, complete genome",
              16899, "Rhynchocyon stuhlmanni", taxid="3036954", genome="mitochondrion",
              source={"specimen_voucher": "RMCA26091"}),
    N.Summary("3", "KX353935.1", "Prunus necrotic ringspot virus isolate SHN-40 coat protein gene, complete cds",
              675, "Prunus necrotic ringspot virus", taxid="37733",
              source={"isolate": "SHN-40", "host": "nectarine", "country": "Iran", "segment": "RNA3"}),
]


def name_settings(settings) -> tuple[list[str], str, bool]:
    """(fields, separator, spaces) from QSettings; unknown keys are dropped."""
    raw = settings.value(SETTING_FIELDS, "") or ""
    fields = [k for k in str(raw).split(",") if k in N.NAME_LABELS] or list(N.DEFAULT_NAME_FIELDS)
    sep = settings.value(SETTING_SEP, " ")
    sep = sep if sep in dict((v, k) for k, v in N.NAME_SEPARATORS) else " "
    spaces = settings.value(SETTING_SPACES, False) in (True, "true")
    return fields, sep, spaces


def is_default(fields: list[str], sep: str) -> bool:
    """NCBI's own "accession definition-line" header: nothing to rename."""
    return fields == N.DEFAULT_NAME_FIELDS and sep == " "


def describe(fields: list[str], sep: str) -> str:
    if is_default(fields, sep):
        return "NCBI default (accession + definition line)"
    joiner = {" ": " · ", "_": " _ ", "|": " | ", "-": " - "}.get(sep, " ")
    return joiner.join(N.NAME_LABELS[k] for k in fields)


class _FieldList(QListWidget):
    """A list whose items can be dragged to the other list or reordered within it.

    Drops are handled here rather than by QListWidget, whose stock drop can overwrite the item
    under the cursor: the drop reports (source list, dragged keys, insert row) through
    `keysDropped` and is accepted as a copy, so the source never removes anything itself."""
    keysDropped = Signal(object, list, int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragDropMode(QAbstractItemView.DragDrop)
        self.setDefaultDropAction(Qt.MoveAction)
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.setDropIndicatorShown(True)
        self.setAlternatingRowColors(True)

    def keys(self) -> list[str]:
        return [self.item(i).data(KEY_ROLE) for i in range(self.count())]

    def dragEnterEvent(self, e):
        if isinstance(e.source(), _FieldList):
            e.setDropAction(Qt.MoveAction)
            e.accept()
        else:
            e.ignore()

    def dragMoveEvent(self, e):
        super().dragMoveEvent(e)
        if isinstance(e.source(), _FieldList):
            e.setDropAction(Qt.MoveAction)
            e.accept()
        else:
            e.ignore()

    def drop_row(self, pos) -> int:
        """Insert position for a drop at viewport point `pos`: before or after the item under it."""
        idx = self.indexAt(pos)
        if not idx.isValid():
            return self.count()
        return idx.row() + (1 if pos.y() > self.visualRect(idx).center().y() else 0)

    def dropEvent(self, e):
        src = e.source()
        if not isinstance(src, _FieldList):
            e.ignore()
            return
        items = sorted(src.selectedItems(), key=src.row)
        keys = [it.data(KEY_ROLE) for it in items]
        row = self.drop_row(e.position().toPoint())
        e.setDropAction(Qt.CopyAction)
        e.accept()
        if keys:
            self.keysDropped.emit(src, keys, row)


def _item(key: str) -> QListWidgetItem:
    it = QListWidgetItem(N.NAME_LABELS[key])
    it.setData(KEY_ROLE, key)
    it.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsDragEnabled)
    return it


class NameFormatDialog(QDialog):
    def __init__(self, parent, settings, samples: list[N.Summary] | None = None):
        super().__init__(parent)
        self.setWindowTitle("NCBI sequence names")
        self.resize(700, 560)
        self.settings = settings
        self.samples = (samples or [])[:3] or EXAMPLES
        fields, sep, spaces = name_settings(settings)

        lay = QVBoxLayout(self)
        intro = QLabel("Choose what the names of imported NCBI sequences are made of. Drag fields from "
                       "<b>Available fields</b> to <b>Use in name</b>, and drag them up or down to set "
                       "their order. Fields a record doesn't have are left out.")
        intro.setWordWrap(True)
        lay.addWidget(intro)

        grid = QGridLayout()
        grid.addWidget(QLabel("Available fields"), 0, 0)
        grid.addWidget(QLabel("Use in name (in this order)"), 0, 2)
        self.available = _FieldList()
        self.chosen = _FieldList()
        for k in fields:
            self.chosen.addItem(_item(k))
        grid.addWidget(self.available, 1, 0)
        grid.addWidget(self.chosen, 1, 2)
        mid = QVBoxLayout()
        mid.addStretch(1)
        self.b_add = QPushButton("Add →"); self.b_add.clicked.connect(self.add_selected)
        self.b_remove = QPushButton("← Remove"); self.b_remove.clicked.connect(self.remove_selected)
        self.b_up = QPushButton("Move up"); self.b_up.clicked.connect(lambda: self.move_selected(-1))
        self.b_down = QPushButton("Move down"); self.b_down.clicked.connect(lambda: self.move_selected(1))
        for b in (self.b_add, self.b_remove, self.b_up, self.b_down):
            mid.addWidget(b)
        mid.addStretch(1)
        grid.addLayout(mid, 1, 1)
        lay.addLayout(grid, 1)
        self.available.itemDoubleClicked.connect(lambda _it: self.add_selected())
        self.chosen.itemDoubleClicked.connect(lambda _it: self.remove_selected())

        opts = QHBoxLayout()
        self.sep = QComboBox()
        for label, value in N.NAME_SEPARATORS:
            self.sep.addItem(label, value)
        self.sep.setCurrentIndex(max(0, self.sep.findData(sep)))
        self.spaces = QCheckBox("Also use it for spaces inside fields")
        self.spaces.setToolTip('e.g. "Pimephales promelas" becomes "Pimephales_promelas"')
        self.spaces.setChecked(spaces)
        opts.addWidget(QLabel("Separator")); opts.addWidget(self.sep); opts.addWidget(self.spaces); opts.addStretch(1)
        lay.addLayout(opts)

        grp = QGroupBox("Preview")
        gl = QVBoxLayout(grp)
        self.preview = QLabel("")
        f = QFont("DejaVu Sans Mono"); f.setStyleHint(QFont.Monospace); self.preview.setFont(f)
        self.preview.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.preview.setWordWrap(True)
        gl.addWidget(self.preview)
        lay.addWidget(grp)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel | QDialogButtonBox.RestoreDefaults)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        bb.button(QDialogButtonBox.RestoreDefaults).setText("NCBI default")
        bb.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self.restore_default)
        lay.addWidget(bb)

        # a drop is applied once the drag has finished, so no list changes under a live drag
        for lst in (self.available, self.chosen):
            lst.keysDropped.connect(
                lambda src, keys, row, dst=lst: QTimer.singleShot(0, lambda: self.place_keys(src, dst, keys, row)))
        self.sep.currentIndexChanged.connect(self._update_preview)
        self.spaces.toggled.connect(self._update_preview)
        self._refresh()

    # ------------------------------------------------------------ list upkeep
    def place_keys(self, src, dst, keys: list[str], row: int):
        """Apply a drop: `keys` dragged from list `src` to position `row` of list `dst`."""
        order = self.chosen.keys()
        if dst is self.chosen:
            if src is self.chosen:
                row -= sum(1 for k in keys if order.index(k) < row)    # rows vacated above the drop point
            rest = [k for k in order if k not in keys]
            row = max(0, min(row, len(rest)))
            order = rest[:row] + keys + rest[row:]
        elif src is self.chosen:
            order = [k for k in order if k not in keys]
        else:
            return                          # within "available": its order is fixed
        self._set_chosen(order, select=keys if dst is self.chosen else ())

    def _set_chosen(self, keys, select=()):
        self.chosen.clear()
        for k in keys:
            it = _item(k)
            self.chosen.addItem(it)
            it.setSelected(k in select)
        self._refresh()

    def _refresh(self):
        """Available = every field not in use, in the canonical order."""
        used = set(self.chosen.keys())
        want = [k for k, _ in N.NAME_FIELDS if k not in used]
        if self.available.keys() != want:
            sel = {it.data(KEY_ROLE) for it in self.available.selectedItems()}
            self.available.clear()
            for k in want:
                it = _item(k)
                self.available.addItem(it)
                it.setSelected(k in sel)
        self._update_preview()

    def _selected(self, lst) -> list[str]:
        return [it.data(KEY_ROLE) for it in sorted(lst.selectedItems(), key=lst.row)]

    def add_selected(self):
        keys = self._selected(self.available)
        if keys:
            self.place_keys(self.available, self.chosen, keys, self.chosen.count())

    def remove_selected(self):
        keys = self._selected(self.chosen)
        if keys:
            self.place_keys(self.chosen, self.available, keys, 0)

    def move_selected(self, step: int):
        order = self.chosen.keys()
        keys = self._selected(self.chosen)
        rows = [order.index(k) for k in keys]
        if not rows or (step < 0 and rows[0] == 0) or (step > 0 and rows[-1] == len(order) - 1):
            return
        for r in (rows if step < 0 else reversed(rows)):
            order[r], order[r + step] = order[r + step], order[r]
        self._set_chosen(order, select=keys)

    def restore_default(self):
        self.sep.setCurrentIndex(self.sep.findData(" "))
        self.spaces.setChecked(False)
        self._set_chosen(N.DEFAULT_NAME_FIELDS)

    # ------------------------------------------------------------ result
    def values(self) -> tuple[list[str], str, bool]:
        return self.chosen.keys() or list(N.DEFAULT_NAME_FIELDS), self.sep.currentData(), self.spaces.isChecked()

    def _update_preview(self, *_):
        fields, sep, spaces = self.values()
        lines = [N.format_name(N.name_fields(s), fields, sep, spaces) for s in self.samples]
        if not self.chosen.count():
            lines.insert(0, "(nothing chosen: NCBI's own names are used)")
        self.preview.setText("\n".join(lines))

    def accept(self):
        fields, sep, spaces = self.values()
        self.settings.setValue(SETTING_FIELDS, ",".join(fields))
        self.settings.setValue(SETTING_SEP, sep)
        self.settings.setValue(SETTING_SPACES, spaces)
        super().accept()
