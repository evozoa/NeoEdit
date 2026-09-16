"""Build a maximum-likelihood tree with IQ-TREE 3.

Modeless: IQ-TREE runs in a QProcess (minutes to hours) with its log streamed into the
dialog, and the alignment stays editable meanwhile (the input was written when the run
started). NeoEdit doesn't draw trees; the result is a Newick file to open in FigTree/iTOL."""
from __future__ import annotations

import os

from PySide6.QtCore import Qt, QProcess, QProcessEnvironment, QUrl
from PySide6.QtGui import QDesktopServices, QFont, QGuiApplication, QTextCursor
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QComboBox, QLineEdit, QLabel,
                               QPlainTextEdit, QPushButton, QSpinBox, QCheckBox, QFileDialog, QMessageBox,
                               QWidget, QProgressBar)

from ...analysis import phylo as PH
from .import_dialog import default_download_dir

ITOL_UPLOAD = "https://itol.embl.de/upload.cgi"


class IQTreeDialog(QDialog):
    def __init__(self, parent, settings, model, selected_rows: list[int], col_range: tuple[int, int] | None):
        super().__init__(parent)
        self.setWindowTitle("Build tree with IQ-TREE")
        self.setModal(False)
        self.resize(760, 720)
        self.settings = settings
        self.model = model
        self.selected_rows = selected_rows
        self.col_range = col_range
        self.proc: QProcess | None = None
        self.run: PH.IQTreeRun | None = None
        self.tree: PH.TreeResult | None = None
        self.seq_type = "protein" if model.seq_type == "protein" else "dna"

        lay = QVBoxLayout(self)
        self.exe = PH.find_iqtree(settings.value("exe/IQ-TREE") or None)
        problem = None
        if self.exe:
            ver = PH.iqtree_version(self.exe)
            problem = PH.version_problem(ver)
            status = f"IQ-TREE {ver}: {self.exe}" + (f"\n{problem}." if problem else "")
        else:
            status = ("IQ-TREE not found.  " + PH.iqtree_install_hint()
                      + "\nSet its location in Edit > Preferences.")
        lab = QLabel(status); lab.setWordWrap(True)
        lay.addWidget(lab)

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

        self.model_box = QComboBox()
        for label, m in PH.MODEL_PRESETS[self.seq_type]:
            self.model_box.addItem(label, m)
        self.model_box.addItem("Other (type an IQ-TREE model below)", "")
        self.model_edit = QLineEdit()
        self.model_edit.setPlaceholderText("e.g. TIM2+F+R3  or  MFP+MERGE")
        saved = settings.value(f"iqtree/model_{self.seq_type}", "MFP")
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
        self.boot.setCurrentIndex(int(settings.value("iqtree/bootstrap", 1)))
        self.reps.setValue(int(settings.value("iqtree/replicates", 1000)))
        self.boot.currentIndexChanged.connect(self._boot_changed)
        self._boot_changed(keep=True)
        brow = QHBoxLayout(); brow.addWidget(self.boot, 1); brow.addWidget(QLabel("replicates")); brow.addWidget(self.reps)
        form.addRow("Branch support", self._wrap(brow))
        self.alrt = QCheckBox("Also SH-aLRT test (1000 replicates)")
        self.alrt.setChecked(settings.value("iqtree/alrt", True) in (True, "true"))
        form.addRow("", self.alrt)

        self.outgroup = QComboBox()
        self.outgroup.addItem("(none — unrooted tree)", -1)
        for i, r in enumerate(model.rows):
            self.outgroup.addItem(r.name, i)
        form.addRow("Outgroup", self.outgroup)

        self.threads = QSpinBox(); self.threads.setRange(0, 256); self.threads.setSpecialValueText(f"auto (up to {PH.auto_threads_max()})")
        self.threads.setToolTip("auto: IQ-TREE times a few thread counts first and picks the fastest, "
                                "using at most all cores but one. On small alignments that timing can take "
                                "longer than the search itself; set a number to skip it.")
        self.threads.setValue(int(settings.value("iqtree/threads", 0)))
        self.seed = QSpinBox(); self.seed.setRange(0, 2**31 - 1); self.seed.setSpecialValueText("random")
        trow = QHBoxLayout(); trow.addWidget(self.threads); trow.addWidget(QLabel("Random seed")); trow.addWidget(self.seed)
        trow.addStretch(1)
        form.addRow("Threads", self._wrap(trow))
        self.extra = QLineEdit()
        form.addRow("Extra IQ-TREE arguments", self.extra)

        base = settings.value("iqtree/out_dir", "") or (os.path.dirname(model.path) if model.path
                                                         else default_download_dir())
        self.out_dir = QLineEdit(base)
        b = QPushButton("…"); b.clicked.connect(self._browse)
        orow = QHBoxLayout(); orow.addWidget(self.out_dir, 1); orow.addWidget(b)
        form.addRow("Save results in", self._wrap(orow))
        stem = os.path.splitext(os.path.basename(model.path))[0] if model.path else "alignment"
        self.stem = PH.safe_stem(stem)
        form.addRow("", QLabel(f"A new folder {self.stem}_iqtree is created there for each run."))
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
        self.b_run.setEnabled(bool(self.exe) and not problem)
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

    @staticmethod
    def _wrap(layout) -> QWidget:
        layout.setContentsMargins(0, 0, 0, 0)
        w = QWidget(); w.setLayout(layout)
        return w

    def _model_changed(self, *_):
        self.model_edit.setEnabled(self.model_box.currentData() == "")

    def _boot_changed(self, *_, keep=False):
        _, flag, default, minimum = PH.BOOTSTRAP_METHODS[self.boot.currentIndex()]
        self.reps.setEnabled(bool(flag))
        if flag:
            self.reps.setMinimum(minimum)
            if not keep:
                self.reps.setValue(default)

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "Save IQ-TREE results in", self.out_dir.text())
        if d:
            self.out_dir.setText(d)

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

    # ------------------------------------------------------------ run
    def running(self) -> bool:
        return self.proc is not None and self.proc.state() != QProcess.NotRunning

    def start(self):
        if self.running():
            return
        if not self._model_string():
            QMessageBox.warning(self, "IQ-TREE", "Type a model name, or pick one from the list."); return
        idx = (self.selected_rows if self.rows_scope.currentIndex() == 1 else list(range(self.model.nrows)))
        boot = self.boot.currentIndex()
        if boot and len(idx) < 4:
            QMessageBox.warning(self, "IQ-TREE", "Branch support needs at least four sequences."); return
        og = self.outgroup.currentData()
        if og is not None and og >= 0 and og not in idx:
            QMessageBox.warning(self, "IQ-TREE", "The outgroup is not among the sequences being analysed."); return
        cols = self.col_range if self.cols_scope.currentIndex() == 1 else None
        self._save_settings()
        try:
            out = PH.unique_dir(os.path.join(self.out_dir.text() or default_download_dir(), self.stem + "_iqtree"))
            self.run = PH.prepare_run([self.model.rows[i].copy() for i in idx], out, self.stem,
                                      self.seq_type, cols)
        except (ValueError, OSError) as e:
            QMessageBox.warning(self, "IQ-TREE", str(e)); return
        args = PH.iqtree_args(self.run, model=self._model_string(), bootstrap=boot, replicates=self.reps.value(),
                              alrt=1000 if self.alrt.isChecked() else 0, threads=self.threads.value(),
                              seed=self.seed.value(), outgroup=[idx.index(og)] if og is not None and og >= 0 else None,
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
        self._set_running(True)
        self.proc.start(self.exe, args)

    def _set_running(self, on: bool):
        self.b_run.setEnabled(not on)
        self.b_stop.setEnabled(on)
        self.busy.setVisible(on)
        self.b_folder.setEnabled(self.run is not None and not on)
        for w in (self.b_copy, self.b_itol):
            w.setEnabled(self.tree is not None and not on)

    def _read(self):
        data = bytes(self.proc.readAllStandardOutput()).decode("utf-8", "replace")
        if data:
            self.log.moveCursor(QTextCursor.End)
            self.log.insertPlainText(data)
            self.log.ensureCursorVisible()

    def _error(self, err):
        if err == QProcess.FailedToStart:
            self._set_running(False)
            self.summary.setText(f"Could not start {self.exe}: {self.proc.errorString()}")

    def _finished(self, code, status):
        self._read()
        if getattr(self, "_stopping", False):
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
        self._set_running(False)

    def stop(self):
        if self.running():
            self._stopping = True
            self.proc.kill()        # IQ-TREE is a console program: terminate() is ignored on Windows

    def closeEvent(self, ev):
        if self.running():
            if QMessageBox.question(self, "IQ-TREE", "IQ-TREE is still running. Stop it?") != QMessageBox.Yes:
                ev.ignore(); return
            self.stop()
            self.proc.waitForFinished(5000)
        super().closeEvent(ev)

    # ------------------------------------------------------------ results
    def _status(self, msg: str):
        parent = self.parent()
        if parent is not None and hasattr(parent, "statusBar"):
            parent.statusBar().showMessage(msg, 8000)

    def open_folder(self):
        if self.run:
            QDesktopServices.openUrl(QUrl.fromLocalFile(self.run.out_dir))

    def copy_tree(self):
        if self.tree:
            QGuiApplication.clipboard().setText(self.tree.newick)
            self._status("Tree copied to the clipboard (Newick)")

    def open_itol(self):
        if self.tree:
            self.copy_tree()
            QDesktopServices.openUrl(QUrl(ITOL_UPLOAD))
