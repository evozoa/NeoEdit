"""Import sequences from NCBI (Entrez) or Ensembl (REST, pinned to a release).

The dialog is modeless. Network work runs in a worker thread; the result is a file in
the downloads folder which is then opened (new alignment) or added to the current one
through the main window's normal file paths."""
from __future__ import annotations

import os

from PySide6.QtCore import Qt, QThread, Signal, QStandardPaths
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QTabWidget, QWidget,
                               QComboBox, QLineEdit, QPlainTextEdit, QPushButton, QLabel, QRadioButton, QSpinBox,
                               QCheckBox, QTableWidget, QTableWidgetItem, QHeaderView, QProgressBar, QFileDialog,
                               QMessageBox, QButtonGroup, QDialogButtonBox, QAbstractItemView, QSizePolicy)

from ...remote import RemoteError, ncbi as N, ensembl as E
from .common import NumItem


class _Worker(QThread):
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, fn, parent=None):
        super().__init__(parent)
        self._fn = fn

    def run(self):
        try:
            self.done.emit(self._fn())
        except RemoteError as e:
            self.failed.emit(str(e))
        except Exception as e:          # noqa: BLE001 - show anything to the user rather than dying silently
            self.failed.emit(f"{type(e).__name__}: {e}")


def default_download_dir() -> str:
    base = QStandardPaths.writableLocation(QStandardPaths.DownloadLocation) or os.path.expanduser("~")
    return os.path.join(base, "NeoEdit")


class ImportDialog(QDialog):
    """`on_open(path)` replaces the alignment, `on_add(path)` appends to it; both are main-window callbacks."""

    def __init__(self, parent, settings, on_open, on_add, autoload: bool = True):
        super().__init__(parent)
        self.setWindowTitle("Import from NCBI / Ensembl")
        self.setModal(False)
        self.resize(760, 640)
        self.settings = settings
        self.on_open, self.on_add = on_open, on_add
        self.autoload = autoload
        self._workers: list[_Worker] = []
        self._species_cache: dict[str, list[E.Species]] = {}
        self._ens_client: E.EnsemblClient | None = None

        lay = QVBoxLayout(self)
        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_ncbi(), "NCBI")
        self.tabs.addTab(self._build_ensembl(), f"Ensembl {E.DEFAULT_RELEASE}")
        self.tabs.currentChanged.connect(self._tab_changed)
        lay.addWidget(self.tabs, 1)

        # destination
        dest = QGroupBox("Destination")
        dl = QFormLayout(dest)
        row = QHBoxLayout()
        self.r_add = QRadioButton("Add to the current alignment")
        self.r_open = QRadioButton("Replace it (open as a new alignment)")
        self.r_add.setChecked(True)
        row.addWidget(self.r_add); row.addWidget(self.r_open); row.addStretch(1)
        dl.addRow("Import into", row)
        row2 = QHBoxLayout()
        self.dir_edit = QLineEdit(settings.value("remote/download_dir", default_download_dir()))
        b = QPushButton("…"); b.setFixedWidth(32); b.clicked.connect(self._browse_dir)
        row2.addWidget(self.dir_edit); row2.addWidget(b)
        dl.addRow("Save downloads in", row2)
        lay.addWidget(dest)

        self.status = QLabel(""); self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(Qt.TextSelectableByMouse)
        lay.addWidget(self.status)
        self.progress = QProgressBar(); self.progress.setRange(0, 0); self.progress.hide()
        lay.addWidget(self.progress)
        bb = QDialogButtonBox()
        self.b_fetch = bb.addButton("Fetch && import", QDialogButtonBox.AcceptRole)
        self.b_fetch.setDefault(True)
        bb.addButton(QDialogButtonBox.Close)
        bb.accepted.connect(self.fetch)
        bb.rejected.connect(self.close)
        lay.addWidget(bb)

    def set_default_destination(self, has_data: bool):
        """Imports add to what is loaded; only an empty editor defaults to 'open as new'."""
        (self.r_add if has_data else self.r_open).setChecked(True)

    # ------------------------------------------------------------ NCBI tab
    def _build_ncbi(self) -> QWidget:
        w = QWidget(); lay = QVBoxLayout(w)
        form = QFormLayout()
        self.n_db = QComboBox()
        for lbl, val in N.DATABASES:
            self.n_db.addItem(lbl, val)
        self.n_fmt = QComboBox()
        for lbl, val in N.FORMATS:
            self.n_fmt.addItem(lbl, val)
        form.addRow("Database", self.n_db)
        form.addRow("Format", self.n_fmt)
        self.n_ids = QPlainTextEdit()
        self.n_ids.setPlaceholderText("Accessions or GI numbers, e.g.  NC_012920.1  MN996528.1  KX353935\n"
                                      "(separate with spaces, commas or new lines)")
        self.n_ids.setFixedHeight(64)
        form.addRow("Accessions", self.n_ids)
        rng = QHBoxLayout()
        self.n_from = QSpinBox(); self.n_from.setRange(0, 2_000_000_000); self.n_from.setSpecialValueText("start")
        self.n_to = QSpinBox(); self.n_to.setRange(0, 2_000_000_000); self.n_to.setSpecialValueText("end")
        self.n_strand = QComboBox(); self.n_strand.addItem("plus strand", 1); self.n_strand.addItem("minus strand", 2)
        for x in (QLabel("from"), self.n_from, QLabel("to"), self.n_to, self.n_strand):
            rng.addWidget(x)
        rng.addStretch(1)
        form.addRow("Sub-range (single accession)", rng)
        lay.addLayout(form)

        grp = QGroupBox("Search Entrez")
        gl = QVBoxLayout(grp)
        srow = QHBoxLayout()
        self.n_term = QLineEdit()
        self.n_term.setPlaceholderText('e.g.  Acantharchus pomotis[Organism] AND mitochondrion[filter] AND "complete genome"')
        self.n_term.returnPressed.connect(self.ncbi_search)
        self.b_search = QPushButton("Search"); self.b_search.clicked.connect(self.ncbi_search)
        self.n_max = QSpinBox(); self.n_max.setRange(1, 1000); self.n_max.setValue(100); self.n_max.setPrefix("max ")
        srow.addWidget(self.n_term, 1); srow.addWidget(self.n_max); srow.addWidget(self.b_search)
        gl.addLayout(srow)
        self.n_table = QTableWidget(0, 5)
        self.n_table.setHorizontalHeaderLabels(["", "Accession", "Length", "Organism", "Title"])
        self.n_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.n_table.horizontalHeader().setStretchLastSection(True)
        self.n_table.verticalHeader().setVisible(False)
        self.n_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.n_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.n_table.setSortingEnabled(True)
        gl.addWidget(self.n_table, 1)
        brow = QHBoxLayout()
        b_all = QPushButton("Check all"); b_all.clicked.connect(lambda: self._check_all(True))
        b_none = QPushButton("Uncheck all"); b_none.clicked.connect(lambda: self._check_all(False))
        self.n_hint = QLabel("Checked results are imported together with the accessions typed above.")
        brow.addWidget(b_all); brow.addWidget(b_none); brow.addWidget(self.n_hint, 1)
        gl.addLayout(brow)
        lay.addWidget(grp, 1)

        ident = QHBoxLayout()
        self.n_email = QLineEdit(self.settings.value("remote/ncbi_email", ""))
        self.n_email.setPlaceholderText("e-mail (NCBI asks clients to identify themselves)")
        self.n_key = QLineEdit(self.settings.value("remote/ncbi_api_key", ""))
        self.n_key.setPlaceholderText("API key (optional, lifts the 3 requests/s limit)")
        self.n_key.setEchoMode(QLineEdit.Password)
        ident.addWidget(QLabel("Identity")); ident.addWidget(self.n_email, 1); ident.addWidget(self.n_key, 1)
        lay.addLayout(ident)
        return w

    def _check_all(self, on: bool):
        for r in range(self.n_table.rowCount()):
            self.n_table.item(r, 0).setCheckState(Qt.Checked if on else Qt.Unchecked)

    def _ncbi_client(self) -> N.NCBIClient:
        self.settings.setValue("remote/ncbi_email", self.n_email.text().strip())
        self.settings.setValue("remote/ncbi_api_key", self.n_key.text().strip())
        return N.NCBIClient(self.n_email.text(), self.n_key.text())

    def ncbi_search(self):
        term = self.n_term.text().strip()
        if not term:
            return
        db, retmax = self.n_db.currentData(), self.n_max.value()
        client = self._ncbi_client()

        def job():
            count, uids = client.search(db, term, retmax)
            return count, client.summaries(db, uids) if uids else []
        self._run(job, self._ncbi_search_done, f"Searching NCBI {db}…")

    def _ncbi_search_done(self, res):
        count, sums = res
        t = self.n_table
        t.setSortingEnabled(False); t.setRowCount(0)
        for s in sums:
            r = t.rowCount(); t.insertRow(r)
            chk = QTableWidgetItem(); chk.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            chk.setCheckState(Qt.Unchecked); chk.setData(Qt.UserRole, s.accession)
            t.setItem(r, 0, chk)
            t.setItem(r, 1, QTableWidgetItem(s.accession))
            t.setItem(r, 2, NumItem(s.length, "{:,}"))
            t.setItem(r, 3, QTableWidgetItem(s.organism))
            t.setItem(r, 4, QTableWidgetItem(s.title))
        t.setSortingEnabled(True)
        shown = len(sums)
        self.status.setText(f"{count:,} hit(s); showing {shown}. Tick the records to import." if count else "No hits.")

    def _checked_accessions(self) -> list[str]:
        out = []
        for r in range(self.n_table.rowCount()):
            it = self.n_table.item(r, 0)
            if it and it.checkState() == Qt.Checked:
                out.append(it.data(Qt.UserRole))
        return out

    def _ncbi_job(self):
        ids = N.parse_ids(self.n_ids.toPlainText())
        for a in self._checked_accessions():
            if a not in ids:
                ids.append(a)
        if not ids:
            raise RemoteError("Enter at least one accession, or search and tick some records.")
        db, fmt = self.n_db.currentData(), self.n_fmt.currentData()
        kw = {}
        if len(ids) == 1:
            if self.n_from.value():
                kw["seq_start"] = self.n_from.value()
            if self.n_to.value():
                kw["seq_stop"] = self.n_to.value()
            if self.n_strand.currentData() == 2:
                kw["strand"] = 2
        client = self._ncbi_client()
        out_dir = self._out_dir()

        def job():
            path, text = client.download(db, ids, fmt, out_dir, **kw)
            n = N.count_records(text, fmt)
            return path, f"NCBI {db}: {n} record(s) → {path}"
        return job, f"Fetching {len(ids)} record(s) from NCBI {db}…"

    # ------------------------------------------------------------ Ensembl tab
    def _build_ensembl(self) -> QWidget:
        w = QWidget(); lay = QVBoxLayout(w)
        form = QFormLayout()
        self.e_server = QComboBox(); self.e_server.setEditable(True)
        self.e_server.addItem(f"Release {E.DEFAULT_RELEASE} (archive server, falls back to live)", "")
        self.e_server.addItem("Live server (current release) — rest.ensembl.org", E.LIVE_SERVER)
        self.e_server.addItem("Human GRCh37 — grch37.rest.ensembl.org", E.GRCH37_SERVER)
        saved = self.settings.value("remote/ensembl_server", "")
        idx = self.e_server.findData(saved)
        if idx >= 0:
            self.e_server.setCurrentIndex(idx)
        elif saved:
            self.e_server.setEditText(saved)
        self.e_server.currentIndexChanged.connect(self._server_changed)
        self.e_server.lineEdit().editingFinished.connect(self._server_changed)
        form.addRow("Server", self.e_server)
        self.e_release = QLabel("release: not checked yet"); self.e_release.setStyleSheet("color: gray")
        form.addRow("", self.e_release)
        srow = QHBoxLayout()
        self.e_div = QComboBox()
        for lbl, val in E.DIVISIONS:
            self.e_div.addItem(lbl, val)
        self.e_div.currentIndexChanged.connect(lambda _i: self._load_species())
        self.e_species = QComboBox(); self.e_species.setEditable(True)
        self.e_species.setInsertPolicy(QComboBox.NoInsert)
        self.e_species.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.e_species.setEditText(self.settings.value("remote/ensembl_species", "homo_sapiens"))
        srow.addWidget(self.e_div); srow.addWidget(self.e_species, 1)
        form.addRow("Species", srow)
        lay.addLayout(form)

        qg = QGroupBox("What to fetch")
        ql = QFormLayout(qg)
        self.e_r_gene = QRadioButton("Gene symbol / stable ID"); self.e_r_region = QRadioButton("Region")
        self.e_r_gene.setChecked(True)
        grp = QButtonGroup(self); grp.addButton(self.e_r_gene); grp.addButton(self.e_r_region)
        self.e_gene = QLineEdit(); self.e_gene.setPlaceholderText("BRCA1, MT-CO1, ENSG00000012048, ENST00000357654, ENSDARG…")
        self.e_gene.textEdited.connect(lambda _t: self.e_r_gene.setChecked(True))
        self.e_gene.returnPressed.connect(self.ensembl_lookup)
        grow = QHBoxLayout(); grow.addWidget(self.e_gene, 1)
        self.b_lookup = QPushButton("Look up"); self.b_lookup.clicked.connect(self.ensembl_lookup); grow.addWidget(self.b_lookup)
        ql.addRow(self.e_r_gene, grow)
        self.e_region = QLineEdit(); self.e_region.setPlaceholderText("17:43044295-43170245   or   MT:1-16569:-1   (1-based, optional strand)")
        self.e_region.textEdited.connect(lambda _t: self.e_r_region.setChecked(True))
        ql.addRow(self.e_r_region, self.e_region)
        self.e_info = QLabel(""); self.e_info.setWordWrap(True); self.e_info.setStyleSheet("color: gray")
        ql.addRow("", self.e_info)
        lay.addWidget(qg)

        sg = QGroupBox("Sequence")
        sl = QFormLayout(sg)
        self.e_type = QComboBox()
        for lbl, val in E.SEQ_TYPES:
            self.e_type.addItem(lbl, val)
        self.e_type.currentIndexChanged.connect(self._type_changed)
        sl.addRow("Type", self.e_type)
        frow = QHBoxLayout()
        self.e_f5 = QSpinBox(); self.e_f5.setRange(0, 1_000_000); self.e_f5.setSuffix(" bp upstream")
        self.e_f3 = QSpinBox(); self.e_f3.setRange(0, 1_000_000); self.e_f3.setSuffix(" bp downstream")
        frow.addWidget(self.e_f5); frow.addWidget(self.e_f3); frow.addStretch(1)
        sl.addRow("Flanks (gene)", frow)
        self.e_annot = QCheckBox("Include gene models (genes, transcripts, CDS) as features"); self.e_annot.setChecked(True)
        self.e_orient = QCheckBox("Orient to the gene's strand (reverse-complement minus-strand genes)"); self.e_orient.setChecked(True)
        self.e_mask = QCheckBox("Soft-mask repeats (lowercase)")
        for cb in (self.e_annot, self.e_orient, self.e_mask):
            sl.addRow("", cb)
        lay.addWidget(sg)
        lay.addStretch(1)
        return w

    def _type_changed(self, _i=None):
        genomic = self.e_type.currentData() == "genomic"
        for x in (self.e_f5, self.e_f3, self.e_annot, self.e_orient, self.e_mask):
            x.setEnabled(genomic)

    def _tab_changed(self, i: int):
        if i == 1 and self.autoload:
            self._ensure_server()
            if not self._species_cache:
                self._load_species()

    def _server_changed(self, *_a):
        self._ens_client = None
        self.e_release.setText("release: not checked yet")
        self.settings.setValue("remote/ensembl_server", self._server_value())
        if self.autoload and self.tabs.currentIndex() == 1:
            self._ensure_server()

    def _server_value(self) -> str:
        data = self.e_server.currentData()
        text = self.e_server.currentText().strip()
        if data is not None and self.e_server.itemText(self.e_server.currentIndex()) == text:
            return data or ""
        return text if text.startswith("http") else ""

    def _ensembl_client(self) -> E.EnsemblClient:
        if self._ens_client is None:
            srv = self._server_value()
            self._ens_client = E.EnsemblClient(release=E.DEFAULT_RELEASE, server=srv or None)
        return self._ens_client

    def _ensure_server(self):
        client = self._ensembl_client()
        if client.release is not None:
            return
        self._run(lambda: client.release_note(), self._release_done, "Contacting Ensembl…", quiet=True)

    def _release_done(self, note: str):
        self.e_release.setText(note)
        bad = "not available" in note
        self.e_release.setStyleSheet("color: #b45309; font-weight: bold" if bad else "color: gray")

    def _load_species(self):
        div = self.e_div.currentData()
        if div in self._species_cache:
            self._fill_species(self._species_cache[div]); return
        client = self._ensembl_client()
        self._run(lambda: (div, client.species(div)), self._species_done, f"Loading {self.e_div.currentText()} species…", quiet=True)

    def _species_done(self, res):
        div, species = res
        self._species_cache[div] = species
        self._fill_species(species)
        if self.e_release.text().startswith("release: not"):
            self._release_done(self._ensembl_client().release_note())

    def _fill_species(self, species: list[E.Species]):
        cur = self.e_species.currentText().strip()
        self.e_species.blockSignals(True)
        self.e_species.clear()
        for s in species:
            self.e_species.addItem(s.label(), s.name)
        self.e_species.blockSignals(False)
        idx = self.e_species.findData(cur)
        if idx >= 0:
            self.e_species.setCurrentIndex(idx)
        else:
            self.e_species.setEditText(cur)

    def _species_value(self) -> str:
        cb = self.e_species
        text = cb.currentText().strip()
        idx = cb.currentIndex()
        if idx >= 0 and cb.itemText(idx) == text and cb.itemData(idx):
            return cb.itemData(idx)
        # typed by hand: accept 'Human', 'homo sapiens', production names
        for i in range(cb.count()):
            if cb.itemText(i).lower().startswith(text.lower() + " ") or cb.itemData(i) == E.normalize_species(text):
                return cb.itemData(i)
        return E.normalize_species(text)

    def ensembl_lookup(self):
        q = self.e_gene.text().strip()
        if not q:
            return
        self.e_r_gene.setChecked(True)
        client = self._ensembl_client(); sp = self._species_value()
        self._run(lambda: client.lookup(q, sp), self._lookup_done, f"Looking up {q}…")

    def _lookup_done(self, lk: E.Lookup):
        L = lk.end - lk.start + 1
        self.e_info.setText(f"{lk.display_name or lk.id} — {lk.object_type} {lk.id}, {lk.biotype}; "
                            f"{lk.species} {lk.assembly_name} {lk.region()} ({L:,} bp)"
                            + (f"\n{lk.description}" if lk.description else ""))
        self.status.setText("")

    def _ensembl_job(self):
        client = self._ensembl_client()
        sp = self._species_value()
        self.settings.setValue("remote/ensembl_species", sp)
        seq_type = self.e_type.currentData()
        out_dir = self._out_dir()
        by_gene = self.e_r_gene.isChecked()
        q = self.e_gene.text().strip() if by_gene else self.e_region.text().strip()
        if not q:
            raise RemoteError("Enter a gene / ID or a region first.")
        f5, f3 = self.e_f5.value(), self.e_f3.value()
        annotate, orient, mask = self.e_annot.isChecked(), self.e_orient.isChecked(), self.e_mask.isChecked()
        if not by_gene:
            chrom, s, e, strand = E.parse_region(q)
            if seq_type != "genomic":
                raise RemoteError("cDNA / CDS / protein need a gene or transcript ID; choose 'Genomic' for a region.")

            def job():
                rec = client.fetch_genomic(sp, chrom, s, e, strand, annotate=annotate, mask="soft" if mask else None)
                path = E.save_record(rec, out_dir)
                return path, f"{client.release_note()}: {chrom}:{s}-{e} ({len(rec.seq):,} bp, {len(rec.features) - 1} features) → {path}"
            return job, f"Fetching {chrom}:{s}-{e} from Ensembl…"

        def job():
            lk = client.lookup(q, sp)
            species = lk.species or sp
            if seq_type == "genomic":
                s, e = lk.start, lk.end
                if lk.strand >= 0:
                    s, e = s - f5, e + f3
                else:
                    s, e = s - f3, e + f5
                s = max(1, s)
                strand = lk.strand if orient else 1
                rec = client.fetch_genomic(species, lk.seq_region_name, s, e, strand, annotate=annotate,
                                           mask="soft" if mask else None, gene=lk)
                path = E.save_record(rec, out_dir)
                return path, (f"{client.release_note()}: {lk.display_name or lk.id} {lk.seq_region_name}:{s}-{e}"
                              f"{'(-)' if strand < 0 else '(+)'} ({len(rec.seq):,} bp, {len(rec.features) - 1} features) → {path}")
            recs = client.fetch_sequences(lk.id, seq_type, lk.display_name, species)
            path = E.save_records(recs, out_dir, f"{lk.display_name or lk.id}_{seq_type}")
            return path, f"{client.release_note()}: {len(recs)} {seq_type} sequence(s) for {lk.display_name or lk.id} → {path}"
        return job, f"Fetching {q} from Ensembl…"

    # ------------------------------------------------------------ shared
    def _out_dir(self) -> str:
        d = self.dir_edit.text().strip() or default_download_dir()
        self.settings.setValue("remote/download_dir", d)
        return d

    def _browse_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Folder for downloaded records", self.dir_edit.text())
        if d:
            self.dir_edit.setText(d)

    def shutdown(self, timeout_ms: int = 3000):
        """Let running workers finish (bounded) so Qt never destroys a live QThread."""
        for w in list(self._workers):
            w.wait(timeout_ms)

    def closeEvent(self, e):
        self.shutdown()
        super().closeEvent(e)

    def _run(self, fn, on_done, busy_text: str, quiet: bool = False):
        w = _Worker(fn, self)
        self._workers.append(w)

        def finished():
            self._workers.remove(w)
            if not self._workers:
                self.progress.hide(); self.b_fetch.setEnabled(True); self.b_search.setEnabled(True)
            w.deleteLater()

        def ok(res):
            finished(); on_done(res)

        def failed(msg):
            finished()
            self.status.setText(msg)
            if not quiet:
                QMessageBox.warning(self, "Import", msg)
        w.done.connect(ok); w.failed.connect(failed)
        if not quiet:
            self.status.setText(busy_text)
        self.progress.show(); self.b_fetch.setEnabled(False); self.b_search.setEnabled(False)
        w.start()

    def fetch(self):
        try:
            job, busy = self._ncbi_job() if self.tabs.currentIndex() == 0 else self._ensembl_job()
        except RemoteError as e:
            QMessageBox.warning(self, "Import", str(e)); return
        self._run(job, self._fetched, busy)

    def _fetched(self, res):
        path, msg = res
        self.status.setText(msg)
        try:
            if self.r_add.isChecked():
                self.on_add(path)
            else:
                self.on_open(path)
        except Exception as e:   # noqa: BLE001
            QMessageBox.critical(self, "Import", f"Downloaded to {path}, but it could not be opened:\n{e}")
