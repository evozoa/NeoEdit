"""Tree-building dialogs (Analysis > Phylogeny): maximum likelihood with IQ-TREE 3, and
neighbor joining computed by NeoEdit itself.

Both are modeless and run in the background (IQ-TREE in a QProcess with its log streamed into
the dialog, NJ in a worker thread), so the alignment stays editable meanwhile: the data are
copied when a run starts. NeoEdit doesn't draw trees; each run writes a Newick file into its
own folder, to open in FigTree / iTOL."""
from __future__ import annotations

import os

from PySide6.QtCore import Qt, QProcess, QProcessEnvironment, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFont, QGuiApplication, QTextCursor
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QComboBox, QLineEdit, QLabel,
                               QPlainTextEdit, QPushButton, QSpinBox, QCheckBox, QFileDialog, QMessageBox,
                               QWidget, QProgressBar)

from ...analysis import nj as NJ
from ...analysis import phylo as PH
from .import_dialog import default_download_dir

ITOL_UPLOAD = "https://itol.embl.de/upload.cgi"


def _wrap(layout) -> QWidget:
    layout.setContentsMargins(0, 0, 0, 0)
    w = QWidget(); w.setLayout(layout)
    return w


class TreeDialogBase(QDialog):
    """Scope (sequences/columns), outgroup, output folder, a log, and the result buttons.
    Subclasses add their options in `build_options(form)` and implement start/stop/running."""

    TITLE = ""
    PROGRAM = ""                 # name used in messages
    SETTINGS = ""                # QSettings group
    FOLDER_SUFFIX = ""           # <stem><suffix>/ per run
    MIN_SUPPORT_SEQS = 4

    def __init__(self, parent, settings, model, selected_rows: list[int], col_range: tuple[int, int] | None):
        super().__init__(parent)
        self.setWindowTitle(self.TITLE)
        self.setModal(False)
        self.resize(760, 720)
        self.settings = settings
        self.model = model
        self.selected_rows = selected_rows
        self.col_range = col_range
        self.out_folder: str | None = None
        self.tree = None                      # result with .newick and .tree_path
        self.seq_type = "protein" if model.seq_type == "protein" else "dna"

        lay = QVBoxLayout(self)
        self.status_label = QLabel(""); self.status_label.setWordWrap(True)
        lay.addWidget(self.status_label)

        form = QFormLayout()
        self.rows_scope = QComboBox()
        self.rows_scope.addItem(f"All {model.nrows} sequences")
        if 3 <= len(selected_rows) < model.nrows:
            self.rows_scope.addItem(f"Selected sequences only ({len(selected_rows)})")
        form.addRow("Sequences", self.rows_scope)
        self.cols_scope = QComboBox()
        self.cols_scope.addItem(f"All {model.width} columns")
        if col_range:
            self.cols_scope.addItem(f"Selected columns only ({col_range[0] + 1}–{col_range[1]})")
        form.addRow("Columns", self.cols_scope)
        form.addRow("Data type", QLabel(("Protein" if self.seq_type == "protein" else "DNA / RNA")
                                        + "   (change with Sequence > Set sequence type)"))

        self.build_options(form)

        self.outgroup = QComboBox()
        self.outgroup.addItem("(none — unrooted tree)", -1)
        for i, r in enumerate(model.rows):
            self.outgroup.addItem(r.name, i)
        form.addRow("Outgroup", self.outgroup)

        base = self.settings.value(f"{self.SETTINGS}/out_dir", "") or (
            os.path.dirname(model.path) if model.path else default_download_dir())
        self.out_dir = QLineEdit(base)
        b = QPushButton("…"); b.clicked.connect(self._browse)
        orow = QHBoxLayout(); orow.addWidget(self.out_dir, 1); orow.addWidget(b)
        form.addRow("Save results in", _wrap(orow))
        stem = os.path.splitext(os.path.basename(model.path))[0] if model.path else "alignment"
        self.stem = PH.safe_stem(stem)
        form.addRow("", QLabel(f"A new folder {self.stem}{self.FOLDER_SUFFIX} is created there for each run."))
        lay.addLayout(form)

        self.log = QPlainTextEdit(); self.log.setReadOnly(True)
        f = QFont("DejaVu Sans Mono"); f.setStyleHint(QFont.Monospace); self.log.setFont(f)
        self.log.setMaximumBlockCount(20000)
        lay.addWidget(self.log, 1)
        self.busy = QProgressBar(); self.busy.setRange(0, 0); self.busy.hide()
        lay.addWidget(self.busy)
        self.summary = QLabel(""); self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        lay.addWidget(self.summary)

        btns = QHBoxLayout()
        self.b_run = QPushButton("Run"); self.b_run.setDefault(True); self.b_run.clicked.connect(self.start)
        self.b_stop = QPushButton("Stop"); self.b_stop.clicked.connect(self.stop); self.b_stop.setEnabled(False)
        self.b_folder = QPushButton("Open folder"); self.b_folder.clicked.connect(self.open_folder)
        self.b_copy = QPushButton("Copy tree"); self.b_copy.clicked.connect(self.copy_tree)
        self.b_copy.setToolTip("Copy the tree (Newick, real sequence names) to the clipboard")
        self.b_itol = QPushButton("View in iTOL…"); self.b_itol.clicked.connect(self.open_itol)
        self.b_itol.setToolTip("Copy the tree and open iTOL's upload page in the browser; paste the tree there")
        for w in (self.b_folder, self.b_copy, self.b_itol):
            w.setEnabled(False)
        close = QPushButton("Close"); close.clicked.connect(self.close)
        for w in (self.b_run, self.b_stop):
            btns.addWidget(w)
        btns.addStretch(1)
        for w in (self.b_folder, self.b_copy, self.b_itol, close):
            btns.addWidget(w)
        lay.addLayout(btns)

    # ------------------------------------------------------------ for subclasses
    def build_options(self, form: QFormLayout):
        pass

    def running(self) -> bool:
        return False

    def start(self):
        raise NotImplementedError

    def stop(self):
        pass

    def wait_stopped(self):
        pass

    def resolve_scope(self, needs_support: bool):
        """-> (row indexes, column range or None, outgroup position within the rows or None),
        or None after telling the user what is wrong."""
        idx = self.selected_rows if self.rows_scope.currentIndex() == 1 else list(range(self.model.nrows))
        if needs_support and len(idx) < self.MIN_SUPPORT_SEQS:
            QMessageBox.warning(self, self.PROGRAM, "Branch support needs at least four sequences.")
            return None
        og = self.outgroup.currentData()
        if og is not None and og >= 0 and og not in idx:
            QMessageBox.warning(self, self.PROGRAM, "The outgroup is not among the sequences being analysed.")
            return None
        cols = self.col_range if self.cols_scope.currentIndex() == 1 else None
        return idx, cols, (idx.index(og) if og is not None and og >= 0 else None)

    def new_output_folder(self) -> str:
        return PH.unique_dir(os.path.join(self.out_dir.text() or default_download_dir(),
                                          self.stem + self.FOLDER_SUFFIX))

    def set_running(self, on: bool):
        self.b_run.setEnabled(not on and self.can_run())
        self.b_stop.setEnabled(on)
        self.busy.setVisible(on)
        self.b_folder.setEnabled(self.out_folder is not None and not on)
        for w in (self.b_copy, self.b_itol):
            w.setEnabled(self.tree is not None and not on)

    def can_run(self) -> bool:
        return True

    # ------------------------------------------------------------ shared behaviour
    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, f"Save {self.PROGRAM} results in", self.out_dir.text())
        if d:
            self.out_dir.setText(d)

    def closeEvent(self, ev):
        if self.running():
            if QMessageBox.question(self, self.PROGRAM, f"{self.PROGRAM} is still running. Stop it?") != QMessageBox.Yes:
                ev.ignore(); return
            self.stop()
            self.wait_stopped()
        super().closeEvent(ev)

    def _status(self, msg: str):
        parent = self.parent()
        if parent is not None and hasattr(parent, "statusBar"):
            parent.statusBar().showMessage(msg, 8000)

    def open_folder(self):
        if self.out_folder:
            QDesktopServices.openUrl(QUrl.fromLocalFile(self.out_folder))

    def copy_tree(self):
        if self.tree:
            QGuiApplication.clipboard().setText(self.tree.newick)
            self._status("Tree copied to the clipboard (Newick)")

    def open_itol(self):
        if self.tree:
            self.copy_tree()
            QDesktopServices.openUrl(QUrl(ITOL_UPLOAD))


# ------------------------------------------------------------------ IQ-TREE
class IQTreeDialog(TreeDialogBase):
    TITLE = "Maximum-likelihood tree (IQ-TREE)"
    PROGRAM = "IQ-TREE"
    SETTINGS = "iqtree"
    FOLDER_SUFFIX = "_iqtree"

    def __init__(self, parent, settings, model, selected_rows, col_range):
        self.proc: QProcess | None = None
        self.run: PH.IQTreeRun | None = None
        self._stopping = False
        self.exe = PH.find_iqtree(settings.value("exe/IQ-TREE") or None)
        self.problem = None
        super().__init__(parent, settings, model, selected_rows, col_range)
        if self.exe:
            ver = PH.iqtree_version(self.exe)
            self.problem = PH.version_problem(ver)
            status = f"IQ-TREE {ver}: {self.exe}" + (f"\n{self.problem}." if self.problem else "")
        else:
            status = ("IQ-TREE not found.  " + PH.iqtree_install_hint()
                      + "\nSet its location in Edit > Preferences.")
        self.status_label.setText(status)
        self.b_run.setEnabled(self.can_run())

    def can_run(self) -> bool:
        return bool(self.exe) and not self.problem

    def build_options(self, form):
        s = self.settings
        self.model_box = QComboBox()
        for label, m in PH.MODEL_PRESETS[self.seq_type]:
            self.model_box.addItem(label, m)
        self.model_box.addItem("Other (type an IQ-TREE model below)", "")
        self.model_edit = QLineEdit()
        self.model_edit.setPlaceholderText("e.g. TIM2+F+R3  or  MFP+MERGE")
        saved = s.value(f"iqtree/model_{self.seq_type}", "MFP")
        i = self.model_box.findData(saved)
        if i < 0:
            i = self.model_box.count() - 1
            self.model_edit.setText(saved)
        self.model_box.setCurrentIndex(i)
        self.model_box.currentIndexChanged.connect(self._model_changed)
        self._model_changed()
        form.addRow("Model", self.model_box)
        form.addRow("", self.model_edit)

        self.boot = QComboBox()
        for label, *_ in PH.BOOTSTRAP_METHODS:
            self.boot.addItem(label)
        self.reps = QSpinBox(); self.reps.setRange(1, 100000)
        self.boot.setCurrentIndex(int(s.value("iqtree/bootstrap", 1)))
        self.reps.setValue(int(s.value("iqtree/replicates", 1000)))
        self.boot.currentIndexChanged.connect(self._boot_changed)
        self._boot_changed(keep=True)
        brow = QHBoxLayout(); brow.addWidget(self.boot, 1); brow.addWidget(QLabel("replicates")); brow.addWidget(self.reps)
        form.addRow("Branch support", _wrap(brow))
        self.alrt = QCheckBox("Also SH-aLRT test (1000 replicates)")
        self.alrt.setChecked(s.value("iqtree/alrt", True) in (True, "true"))
        form.addRow("", self.alrt)

        self.threads = QSpinBox(); self.threads.setRange(0, 256)
        self.threads.setSpecialValueText(f"auto (up to {PH.auto_threads_max()})")
        self.threads.setToolTip("auto: IQ-TREE times a few thread counts first and picks the fastest, "
                                "using at most all cores but one. On small alignments that timing can take "
                                "longer than the search itself; set a number to skip it.")
        self.threads.setValue(int(s.value("iqtree/threads", 0)))
        self.seed = QSpinBox(); self.seed.setRange(0, 2**31 - 1); self.seed.setSpecialValueText("random")
        trow = QHBoxLayout(); trow.addWidget(self.threads); trow.addWidget(QLabel("Random seed")); trow.addWidget(self.seed)
        trow.addStretch(1)
        form.addRow("Threads", _wrap(trow))
        self.extra = QLineEdit()
        form.addRow("Extra IQ-TREE arguments", self.extra)

    def _model_changed(self, *_):
        self.model_edit.setEnabled(self.model_box.currentData() == "")

    def _boot_changed(self, *_, keep=False):
        _, flag, default, minimum = PH.BOOTSTRAP_METHODS[self.boot.currentIndex()]
        self.reps.setEnabled(bool(flag))
        if flag:
            self.reps.setMinimum(minimum)
            if not keep:
                self.reps.setValue(default)

    def _model_string(self) -> str:
        return self.model_box.currentData() or self.model_edit.text().strip()

    def _save_settings(self):
        s = self.settings
        s.setValue(f"iqtree/model_{self.seq_type}", self._model_string() or "MFP")
        s.setValue("iqtree/bootstrap", self.boot.currentIndex())
        s.setValue("iqtree/replicates", self.reps.value())
        s.setValue("iqtree/alrt", self.alrt.isChecked())
        s.setValue("iqtree/threads", self.threads.value())
        s.setValue("iqtree/out_dir", self.out_dir.text())

    def running(self) -> bool:
        return self.proc is not None and self.proc.state() != QProcess.NotRunning

    def start(self):
        if self.running():
            return
        if not self._model_string():
            QMessageBox.warning(self, "IQ-TREE", "Type a model name, or pick one from the list."); return
        boot = self.boot.currentIndex()
        scope = self.resolve_scope(needs_support=bool(boot))
        if scope is None:
            return
        idx, cols, og = scope
        self._save_settings()
        try:
            out = self.new_output_folder()
            self.run = PH.prepare_run([self.model.rows[i].copy() for i in idx], out, self.stem,
                                      self.seq_type, cols)
        except (ValueError, OSError) as e:
            QMessageBox.warning(self, "IQ-TREE", str(e)); return
        self.out_folder = out
        args = PH.iqtree_args(self.run, model=self._model_string(), bootstrap=boot, replicates=self.reps.value(),
                              alrt=1000 if self.alrt.isChecked() else 0, threads=self.threads.value(),
                              seed=self.seed.value(), outgroup=[og] if og is not None else None,
                              extra_args=self.extra.text().split())
        self.tree = None
        self.log.clear()
        self.log.appendPlainText("$ " + " ".join([self.exe] + args) + "\n")
        self.summary.setText(f"Running on {len(idx)} sequences × {self.run.nsites} columns → {out}")
        self.proc = QProcess(self)
        self.proc.setWorkingDirectory(out)
        self.proc.setProcessChannelMode(QProcess.MergedChannels)
        env = QProcessEnvironment.systemEnvironment()
        for k, v in PH.IQTREE_ENV.items():
            env.insert(k, v)
        self.proc.setProcessEnvironment(env)
        self.proc.readyReadStandardOutput.connect(self._read)
        self.proc.finished.connect(self._finished)
        self.proc.errorOccurred.connect(self._error)
        self.set_running(True)
        self.proc.start(self.exe, args)

    def _read(self):
        data = bytes(self.proc.readAllStandardOutput()).decode("utf-8", "replace")
        if data:
            self.log.moveCursor(QTextCursor.End)
            self.log.insertPlainText(data)
            self.log.ensureCursorVisible()

    def _error(self, err):
        if err == QProcess.FailedToStart:
            self.set_running(False)
            self.summary.setText(f"Could not start {self.exe}: {self.proc.errorString()}")

    def _finished(self, code, status):
        self._read()
        if self._stopping:
            self._stopping = False
            self.summary.setText(f"Stopped. Partial output is in {self.run.out_dir}")
        elif status != QProcess.NormalExit or code != 0:
            self.summary.setText(f"IQ-TREE failed (exit {code}). See the log above and {self.run.prefix}.log")
        else:
            try:
                self.tree = PH.finish_run(self.run)
            except Exception as e:       # noqa: BLE001
                self.summary.setText(str(e))
            else:
                r = self.tree
                bits = [f"Tree: {r.tree_path}"]
                if r.best_model:
                    bits.append(f"Model: {r.best_model}")
                if r.log_likelihood:
                    bits.append(f"log-likelihood {r.log_likelihood}")
                self.summary.setText("Finished.  " + "   ".join(bits)
                                     + f"\nFull report: {r.report_path}")
                self._status(f"IQ-TREE finished: {r.tree_path}")
        self.set_running(False)

    def stop(self):
        if self.running():
            self._stopping = True
            self.proc.kill()        # IQ-TREE is a console program: terminate() is ignored on Windows

    def wait_stopped(self):
        if self.proc is not None:
            self.proc.waitForFinished(5000)


# ------------------------------------------------------------------ Neighbor joining
class _NJWorker(QThread):
    progress = Signal(int, int)
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, kwargs, parent=None):
        super().__init__(parent)
        self.kwargs = kwargs
        self.cancel = False

    def run(self):
        try:
            res = NJ.run_nj(**self.kwargs, progress=lambda d, t: self.progress.emit(d, t),
                            cancelled=lambda: self.cancel)
        except NJ.Cancelled:
            self.failed.emit("")
        except (ValueError, OSError) as e:
            self.failed.emit(str(e))
        except Exception as e:          # noqa: BLE001 - show anything rather than dying silently
            self.failed.emit(f"{type(e).__name__}: {e}")
        else:
            self.done.emit(res)


class NJDialog(TreeDialogBase):
    TITLE = "Neighbor-joining tree"
    PROGRAM = "Neighbor joining"
    SETTINGS = "nj"
    FOLDER_SUFFIX = "_nj"

    def __init__(self, parent, settings, model, selected_rows, col_range):
        self.worker: _NJWorker | None = None
        super().__init__(parent, settings, model, selected_rows, col_range)
        self.status_label.setText("Neighbor joining (Saitou & Nei 1987), computed by NeoEdit. "
                                  "Gaps and ambiguous bases are treated as missing data.")
        self.log.setPlaceholderText("The run report appears here.")

    def build_options(self, form):
        s = self.settings
        self.dist_box = QComboBox()
        for label, code in (NJ.PROTEIN_MODELS if self.seq_type == "protein" else NJ.DNA_MODELS):
            self.dist_box.addItem(label, code)
        i = self.dist_box.findData(s.value(f"nj/model_{self.seq_type}",
                                           "poisson" if self.seq_type == "protein" else "k2p"))
        self.dist_box.setCurrentIndex(max(0, i))
        form.addRow("Distance", self.dist_box)
        self.gap_box = QComboBox()
        for label, code in NJ.GAP_MODES:
            self.gap_box.addItem(label, code)
        self.gap_box.setCurrentIndex(max(0, self.gap_box.findData(s.value("nj/gaps", "pairwise"))))
        self.gap_box.setToolTip("Pairwise: each pair of sequences uses every site both of them have.\n"
                                "Complete: only sites with no gap or ambiguity in any sequence are used.")
        form.addRow("Gaps / missing data", self.gap_box)
        self.boot = QCheckBox("Bootstrap")
        self.boot.setChecked(s.value("nj/bootstrap", True) in (True, "true"))
        self.reps = QSpinBox(); self.reps.setRange(10, 100000)
        self.reps.setValue(int(s.value("nj/replicates", 1000)))
        self.boot.toggled.connect(self.reps.setEnabled)
        self.reps.setEnabled(self.boot.isChecked())
        self.seed = QSpinBox(); self.seed.setRange(0, 2**31 - 1); self.seed.setSpecialValueText("random")
        brow = QHBoxLayout()
        for w in (self.boot, self.reps, QLabel("replicates"), QLabel("   Random seed"), self.seed):
            brow.addWidget(w)
        brow.addStretch(1)
        form.addRow("Branch support", _wrap(brow))

    def _save_settings(self):
        s = self.settings
        s.setValue(f"nj/model_{self.seq_type}", self.dist_box.currentData())
        s.setValue("nj/gaps", self.gap_box.currentData())
        s.setValue("nj/bootstrap", self.boot.isChecked())
        s.setValue("nj/replicates", self.reps.value())
        s.setValue("nj/out_dir", self.out_dir.text())

    def running(self) -> bool:
        return self.worker is not None and self.worker.isRunning()

    def start(self):
        if self.running():
            return
        reps = self.reps.value() if self.boot.isChecked() else 0
        scope = self.resolve_scope(needs_support=bool(reps))
        if scope is None:
            return
        idx, cols, og = scope
        self._save_settings()
        self.out_folder = None
        self.tree = None
        self.log.clear()
        out = self.new_output_folder()
        kwargs = dict(rows=[self.model.rows[i].copy() for i in idx], out_dir=out, stem=self.stem,
                      seq_type=self.seq_type, model=self.dist_box.currentData(),
                      gaps=self.gap_box.currentData(), replicates=reps, seed=self.seed.value(),
                      outgroup=og, columns=cols)
        self.summary.setText(f"Computing distances for {len(idx)} sequences…")
        self.busy.setRange(0, 0)
        self.worker = _NJWorker(kwargs, self)
        self.worker.progress.connect(self._progress)
        self.worker.done.connect(lambda res, out=out: self._done(res, out))
        self.worker.failed.connect(self._failed)
        self.worker.finished.connect(lambda: self.set_running(False))
        self.set_running(True)
        self.worker.start()

    def _progress(self, done, total):
        self.busy.setRange(0, total)
        self.busy.setValue(done)
        self.summary.setText(f"Bootstrap replicate {done} of {total}…")

    def _done(self, res, out):
        self.tree = res
        self.out_folder = out
        with open(res.report_path, encoding="utf-8") as fh:
            self.log.setPlainText(fh.read())
        self.summary.setText(f"Finished.  Tree: {res.tree_path}\nDistances: {res.distance_path}")
        self._status(f"Neighbor-joining tree written: {res.tree_path}")

    def _failed(self, msg):
        if not msg:
            self.summary.setText("Stopped.")
            return
        self.summary.setText("Could not build the tree.")
        QMessageBox.warning(self, "Neighbor joining", msg)

    def stop(self):
        if self.running():
            self.worker.cancel = True

    def wait_stopped(self):
        if self.worker is not None:
            self.worker.wait(10000)
