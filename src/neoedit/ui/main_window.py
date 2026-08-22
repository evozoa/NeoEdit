from __future__ import annotations

import os
import sys
import traceback

from PySide6.QtCore import Qt, QSettings, QSize
from PySide6.QtGui import QAction, QActionGroup, QKeySequence, QIcon
from PySide6.QtWidgets import (QMainWindow, QFileDialog, QMessageBox, QLabel, QToolBar, QDockWidget,
                               QApplication, QInputDialog, QMenu, QTabWidget, QWidget, QVBoxLayout, QComboBox,
                               QSpinBox, QProgressDialog, QSplitter)

from ..model import AlignmentModel, SequenceRow, Feature
from ..model import io as mio
from ..model import colors as C
from ..analysis import translate as T
from ..analysis import external as EXT
from .alignment_view import AlignmentView
from .feature_track import FeaturePanel
from .genome_panel import GenomePanel
from .circular_view import CircularView
from ..genome.fasta_index import IndexedFasta
from ..genome import annotations as GA
from ..genome.projection import RefProjection
from .icons import icon, app_icon as icon_app
from .dialogs.translate_dialog import TranslateDialog
from .dialogs.orf_dialog import ORFFinderDialog
from .dialogs.primer_dialog import PrimerDialog
from .dialogs.design_dialog import DesignDialog
from .dialogs.restriction_dialog import RestrictionDialog
from .dialogs.misc_dialogs import (FindDialog, StatsDialog, IdentityDialog, PlotDialog, AlignDialog,
                                   PreferencesDialog, NewSequenceDialog, ConsensusDialog)
from .dialogs.common import TextDialog
from .dialogs.import_dialog import ImportDialog
from .. import __version__

FILE_FILTER = ";;".join(
    ["All sequence files (*.fasta *.fas *.fa *.fna *.faa *.aln *.phy *.nex *.nexus *.sto *.gb *.gbk *.embl *.msf *.bio *.txt)"]
    + [f"{lbl} (" + " ".join("*" + e for e in exts) + ")" for lbl, _, exts, _ in mio.FORMATS]
    + ["All files (*)"])


class MainWindow(QMainWindow):
    def __init__(self, paths: list[str] | None = None):
        super().__init__()
        self.setWindowTitle("NeoEdit")
        self.setWindowIcon(icon_app())
        self.resize(1200, 750)
        self.settings = QSettings("neoedit", os.environ.get("NEOEDIT_SETTINGS", "neoedit"))
        self.model = AlignmentModel()
        self.view = AlignmentView(self.model)
        # genome context panel (tiers 1-2) above the grid (tier 3); hidden until a genome is opened
        self.genome_panel = GenomePanel()
        self.genome_panel.hide()
        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.addWidget(self.genome_panel)
        self.splitter.addWidget(self.view)
        self.splitter.setStretchFactor(0, 0); self.splitter.setStretchFactor(1, 1)
        self.setCentralWidget(self.splitter)
        self.genome: IndexedFasta | None = None
        self.genome_contig: str | None = None
        self._proj: RefProjection | None = None
        self._proj_row: int = -1
        self._ref_ungapped: str | None = None
        self.annotation: GA.Annotation | None = None
        self.synteny_blocks: list = []
        self.genome_panel.contigSelected.connect(self._genome_load_contig)
        self.genome_panel.focusRequested.connect(self._genome_focus)
        self.genome_panel.geneActivated.connect(self._genome_open_gene)
        self.genome_panel.openRegionRequested.connect(self._genome_open_region)
        self.genome_panel.insertionClicked.connect(self._goto_columns)
        self.view.horizontalScrollBar().valueChanged.connect(self._grid_scrolled)
        self.view.horizontalScrollBar().rangeChanged.connect(lambda *_: self._grid_scrolled())
        self.find_dlg = None
        self._children = []  # keep non-modal dialogs alive

        self.features_panel = FeaturePanel(self.model)
        self.dock = QDockWidget("Features", self)
        self.dock.setWidget(self.features_panel)
        self.addDockWidget(Qt.BottomDockWidgetArea, self.dock)
        self.dock.hide()
        self.features_panel.featureSelected.connect(self._goto_feature)
        self.features_panel.featuresChanged.connect(self.view.viewport().update)

        self._build_actions()
        self._build_menus()
        self._build_toolbar()
        self._build_status()
        self.view.cursorChanged.connect(self._update_status)
        self.view.selectionChanged.connect(self._update_status)
        self.view.modeChanged.connect(self._update_status)
        self.view.featureActivated.connect(lambda f: None)
        self.model.add_listener(self._on_model)
        self.view.translation_provider = self._translation_regions
        self._update_status()
        self._restore()
        self.view.trans_table = int(self.settings.value("default_table", 1))
        if paths:
            for p in paths:
                self.open_path(p)

    # ------------------------------------------------------------ actions
    def _act(self, text, slot, shortcut=None, checkable=False, tip=None):
        a = QAction(text, self)
        if isinstance(shortcut, (list, tuple)):      # several bindings; the first is shown in menus
            a.setShortcuts([QKeySequence(k) for k in shortcut])
        elif shortcut:
            a.setShortcut(QKeySequence(shortcut))
        a.setCheckable(checkable)
        if tip:
            a.setStatusTip(tip)
        if slot:
            a.triggered.connect(slot)
        return a

    def _build_actions(self):
        A = self._act
        self.a_new = A("&New alignment", self.new_alignment, "Ctrl+N")
        self.a_open = A("&Open…", self.open_file, "Ctrl+O")
        self.a_import = A("&Import sequences into current…", self.import_file)
        self.a_import_remote = A("Import from &NCBI / Ensembl / UCSC…", self.import_remote, "Ctrl+Shift+I",
                                 tip="Fetch records from NCBI Entrez, Ensembl (pinned to release 116) or the UCSC Genome Browser by accession, gene, region or search")
        self.a_save = A("&Save", self.save_file, "Ctrl+S")
        self.a_saveas = A("Save &As…", self.save_file_as, "Ctrl+Shift+S")
        self.a_export_sel = A("Export selected sequences…", self.export_selected)
        self.a_quit = A("&Quit", self.close, "Ctrl+Q")

        self.a_undo = A("&Undo", self.undo, "Ctrl+Z")
        self.a_redo = A("&Redo", self.redo, "Ctrl+Shift+Z")
        self.a_copy = A("&Copy (FASTA)", self.copy_sel, "Ctrl+C")
        self.a_copy_raw = A("Copy sequence text only", self.copy_raw, "Ctrl+Shift+C")
        self.a_paste = A("&Paste sequences", self.paste_seqs, "Ctrl+V")
        self.a_selall = A("Select &all", self.view.select_all, "Ctrl+A")
        self.a_find = A("&Find…", self.find, "Ctrl+F")
        self.a_findnext = A("Find &next", self.find_next, "F3")
        self.a_goto = A("&Go to position…", self.goto, "Ctrl+G")
        self.a_prefs = A("&Preferences…", self.preferences)

        self.mode_group = QActionGroup(self)
        self.a_mode_slide = A("&Select / Slide mode", lambda: self._set_mode("slide"), "F5", True,
                              "Box-select residues; drag a selection to slide it over gaps (Shift: move whole downstream sequence)")
        self.a_mode_edit = A("&Edit mode", lambda: self._set_mode("edit"), "F6", True,
                             "Place the cursor and type residues (Insert / Overwrite)")
        self.a_mode_grab = A("&Grab && Drag mode", lambda: self._set_mode("grab"), "F7", True,
                             "Grab a single residue and drag it (Shift: drag everything downstream of it)")
        for a in (self.a_mode_slide, self.a_mode_edit, self.a_mode_grab):
            self.mode_group.addAction(a)
        self.a_mode_slide.setChecked(True)
        self.sub_group = QActionGroup(self)
        self.a_mode_insert = A("&Insert", lambda: self.view.set_edit_submode("insert"), None, True, "Typing inserts residues (Edit mode)")
        self.a_mode_over = A("&Overwrite", lambda: self.view.set_edit_submode("overwrite"), None, True, "Typing overwrites residues (Edit mode)")
        self.sub_group.addAction(self.a_mode_insert); self.sub_group.addAction(self.a_mode_over)
        self.a_mode_insert.setChecked(True)
        self.a_downstream = A("Sliding moves &downstream sequence by default", self._toggle_downstream, None, True,
                              "When on, sliding/grabbing moves the entire sequence downstream (Shift reverses); when off, gaps are crunched/opened (Shift reverses)")

        self.a_zoom_in = A("Zoom &in", lambda: self.view.zoom(1), "Ctrl+=")
        self.a_zoom_out = A("Zoom &out", lambda: self.view.zoom(-1), "Ctrl+-")
        self.a_font = A("&Font…", self.choose_font)
        self.a_crisp = A("&Crisp text (no anti-aliasing)", lambda on: self.view.set_text_style(crisp=on), None, True); self.a_crisp.setChecked(True)
        self.weight_group = QActionGroup(self)
        self.a_w_reg = A("Regular", lambda: self.view.set_text_style(weight="regular"), None, True)
        self.a_w_semi = A("Semi-bold", lambda: self.view.set_text_style(weight="semibold"), None, True)
        self.a_w_bold = A("Bold", lambda: self.view.set_text_style(weight="bold"), None, True)
        for a in (self.a_w_reg, self.a_w_semi, self.a_w_bold):
            self.weight_group.addAction(a)
        self.a_w_reg.setChecked(True)
        self.a_row_more = A("Increase &line spacing", lambda: self.view.set_spacing(row_pad=self.view.row_pad + 1), "Ctrl+Shift+Up")
        self.a_row_less = A("Decrease line spacing", lambda: self.view.set_spacing(row_pad=self.view.row_pad - 1), "Ctrl+Shift+Down")
        self.a_col_more = A("Increase &character spacing", lambda: self.view.set_spacing(col_pad=self.view.col_pad + 1), "Ctrl+Shift+Right")
        self.a_col_less = A("Decrease character spacing", lambda: self.view.set_spacing(col_pad=self.view.col_pad - 1), "Ctrl+Shift+Left")
        self.a_pinned_ref = A("Show pinned &reference row", self._toggle_pinned_ref, None, True,
                              tip="Keep the reference sequence (the row the chromosome map and gene models follow) visible above the grid")
        self.a_pinned_ref.setChecked(True)
        self.a_translation = A("Show &translation under DNA", self._toggle_translation, "Ctrl+T", True)
        self.a_features = A("Show &features", self._toggle_features, None, True); self.a_features.setChecked(True)
        self.a_shade_feat = A("Shade feature &regions", self._toggle_shade_features, None, True,
                              "Tint residues covered by a feature; off by default so gene models "
                              "do not wash out the sequence")
        self.a_dock = self.dock.toggleViewAction(); self.a_dock.setText("Features &panel")
        self.a_dots = A("Identities as &dots", self._toggle_dots, None, True)
        self.view_group = QActionGroup(self)
        self.a_normal = A("&Normal view (colored letters)", lambda: self._set_inverse(False), None, True)
        self.a_inverse = A("In&verse view (colored backgrounds)", lambda: self._set_inverse(True), "Ctrl+I", True)
        self.a_normal.setIcon(icon("normal_view")); self.a_inverse.setIcon(icon("inverse_view"))
        self.view_group.addAction(self.a_normal); self.view_group.addAction(self.a_inverse)
        self.a_normal.setChecked(True)
        self.a_toggle_inverse = A("Toggle inverse view", lambda: self._set_inverse(not self.a_inverse.isChecked()), "Ctrl+I")
        self.a_inverse.setShortcut(QKeySequence())
        self.color_group = QActionGroup(self)
        self.a_col_scheme = A("Color by &scheme", lambda: self._set_color_mode("scheme"), None, True)
        self.a_col_ident = A("Shade &identities", lambda: self._set_color_mode("identity"), None, True)
        self.a_col_none = A("&No color", lambda: self._set_color_mode("none"), None, True)
        for a in (self.a_col_scheme, self.a_col_ident, self.a_col_none):
            self.color_group.addAction(a)
        self.a_col_scheme.setChecked(True)
        self.ref_group = QActionGroup(self)
        self.a_ref_cons = A("…relative to consensus", lambda: self._set_ref("consensus"), None, True)
        self.a_ref_first = A("…relative to first sequence", lambda: self._set_ref("first"), None, True)
        self.ref_group.addAction(self.a_ref_cons); self.ref_group.addAction(self.a_ref_first)
        self.a_ref_cons.setChecked(True)
        self.a_threshold = A("Consensus threshold (identity coloring)…", self._set_threshold)

        # sequence ops
        self.a_revcomp = A("&Reverse complement", lambda: self.model.reverse_complement(self.view.target_rows()),
                           ("Ctrl+Shift+R", "Ctrl+R"))      # Ctrl+Shift+R as in BioEdit; Ctrl+R kept as an alias
        self.a_rev = A("Re&verse", lambda: self.model.reverse(self.view.target_rows()))
        self.a_comp = A("&Complement", lambda: self.model.complement(self.view.target_rows()))
        self.a_upper = A("&Uppercase", lambda: self.model.to_upper(self.view.target_rows()))
        self.a_lower = A("&Lowercase", lambda: self.model.to_lower(self.view.target_rows()))
        self.a_rmgaps = A("Remove &gaps from sequence(s)", lambda: self.model.remove_gaps(self.view.target_rows()))
        self.a_rna = A("DNA → RNA", lambda: self.model.dna_to_rna(self.view.target_rows()))
        self.a_dna = A("RNA → DNA", lambda: self.model.rna_to_dna(self.view.target_rows()))
        self.a_translate = A("&Translate…", self.translate, "Ctrl+Shift+T")
        self.a_sixframe = A("Six-frame translation report", self.six_frame)
        self.a_rename = A("Re&name sequence…", self.rename)
        self.a_newseq = A("&New sequence…", self.new_sequence, "Ctrl+Shift+N")
        self.a_delseq = A("&Delete sequence(s)", self.delete_seqs)
        self.a_dupseq = A("Du&plicate sequence(s)", lambda: self.model.duplicate_rows(self.view.target_rows()))
        self.a_up = A("Move sequence(s) &up", lambda: self._move(-1), "Ctrl+Up")
        self.a_down = A("Move sequence(s) &down", lambda: self._move(1), "Ctrl+Down")
        self.a_settype = A("Set sequence type…", self.set_type)
        self.a_pin = A("&Pin as reference", self.pin_reference, "Ctrl+Shift+E",
                       tip="Anchor the gene models, chromosome map and coordinates to this sequence")
        self.a_circular = A("&Circular molecule", self._toggle_circular, None, True,
                            tip="Treat the sequence as circular (mitogenome, plasmid)")
        self.a_grp_set = A("&Assign selected to group…", self.group_assign, "Ctrl+Shift+A")
        self.a_grp_clear = A("&Remove selected from group", lambda: self.model.set_group(self.view.target_rows(), ""))
        self.a_grp_color = A("Set group &color…", self.group_color)
        self.a_blast = A("BLAST selected sequence (NCBI web)", self.blast)

        # alignment ops
        self.a_align = A("&Align with MAFFT…", self.align_external, "Ctrl+M")
        self.a_rm_gapcols = A("Remove gap-only &columns", self.model_call("remove_gap_only_columns"))
        self.a_pad = A("&Pad sequences to equal length", self.model_call("pad_to_equal_length"))
        self.a_insgapcol = A("Insert gap column at cursor", lambda: self.model.insert_gap_columns(self.view.cur_col, 1), "Ctrl+Space")
        self.a_delgapcol = A("Delete gap column at cursor", lambda: self.model.delete_gap_columns(self.view.cur_col, 1), "Ctrl+Delete")
        self.a_extract = A("Extract selected columns to new window", self.extract_cols)

        # analysis
        self.a_orf = A("&ORF finder…", self.orf_finder, "Ctrl+Shift+O")
        self.a_primer = A("&Primer design (single template)…", self.primer_design, "Ctrl+Shift+P")
        self.a_design = A("&Conserved / discriminating primers…", self.design_primers, "Ctrl+Shift+D",
                          tip="Design primers across the alignment: universal for a group, or discriminating against another group (eDNA)")
        self.a_restrict = A("&Restriction sites…", self.restriction_sites, "Ctrl+Shift+X",
                            tip="Restriction map of one sequence, or which enzymes distinguish sequences in the alignment")
        self.a_stats = A("Sequence &statistics", self.stats)
        self.a_ident = A("&Identity matrix", self.identity)
        self.a_plot = A("&Conservation / entropy plot", self.plot)
        self.a_consensus = A("&Consensus sequence…", self.consensus_tool,
                             tip="Majority or IUPAC-degenerate consensus of all or selected sequences; add it as a row, copy or save it")

        self.a_about = A("&About", self.about)

        # genome
        self.a_g_ref = A("Open Gen&Bank reference (mitogenome/plasmid)…", self.genome_open_reference,
                         tip="Open an annotated GenBank record: sequence becomes the reference row, its features populate the gene view")
        self.a_g_add = A("Add se&quences anchored to reference…", self.genome_add_anchored,
                         tip="MAFFT --add --keeplength: align new sequences to the current alignment without changing its columns")
        self.a_g_open = A("Open &genome FASTA…", self.genome_open, "Ctrl+Shift+G",
                          tip="Open a (large, indexed) genome FASTA and pick a contig/chromosome")
        self.a_g_ann = A("Load &annotation (GFF3/GTF/BED/GenBank)…", self.genome_load_annotation)
        self.a_g_syn = A("Load &synteny blocks (PAF)…", self.genome_load_synteny)
        self.a_g_panel = A("Show genome &panel", self._toggle_genome_panel, "Ctrl+Shift+B", True)
        self.a_g_goto = A("Go to &region / gene…", self.genome_goto, "Ctrl+J")
        self.a_g_openreg = A("Open current region in &new editor window", lambda: self._genome_open_region(*self.genome_panel.window()))
        self.a_g_clear = A("&Close genome (keep sequence)", self.genome_close)
        self.a_g_orfclear = A("Clear ORF &tracks", lambda: self.genome_panel.clear_orf_tracks())
        self.a_g_circ = A("Circular &map (mitogenome / plasmid)…", self.circular_map, "Ctrl+Shift+M",
                          tip="Circular gene map with GC content and GC skew; export as SVG or PNG")

        for name, act in (("open", self.a_open), ("save", self.a_save), ("undo", self.a_undo), ("redo", self.a_redo),
                          ("mode_slide", self.a_mode_slide), ("mode_edit", self.a_mode_edit), ("mode_grab", self.a_mode_grab),
                          ("mode_insert", self.a_mode_insert), ("mode_overwrite", self.a_mode_over), ("downstream", self.a_downstream), ("zoom_in", self.a_zoom_in), ("zoom_out", self.a_zoom_out),
                          ("translation", self.a_translation), ("orf", self.a_orf), ("primer", self.a_primer),
                          ("align", self.a_align), ("features", self.a_dock),
                          ("row_more", self.a_row_more), ("row_less", self.a_row_less),
                          ("col_more", self.a_col_more), ("col_less", self.a_col_less)):
            act.setIcon(icon(name))
        self.a_toggle_inverse.setShortcutContext(Qt.WindowShortcut)
        self.addAction(self.a_toggle_inverse)

    def model_call(self, name):
        return lambda: getattr(self.model, name)()

    def _build_menus(self):
        mb = self.menuBar()
        f = mb.addMenu("&File")
        for a in (self.a_new, self.a_open, self.a_import, self.a_import_remote, None, self.a_save, self.a_saveas, self.a_export_sel, None):
            f.addAction(a) if a else f.addSeparator()
        self.recent_menu = f.addMenu("Open &recent")
        f.addSeparator(); f.addAction(self.a_quit)

        e = mb.addMenu("&Edit")
        for a in (self.a_undo, self.a_redo, None, self.a_copy, self.a_copy_raw, self.a_paste, self.a_selall, None,
                  self.a_find, self.a_findnext, self.a_goto, None,
                  self.a_mode_slide, self.a_mode_edit, self.a_mode_grab, self.a_mode_insert, self.a_mode_over, self.a_downstream, None, self.a_prefs):
            e.addAction(a) if a else e.addSeparator()

        v = mb.addMenu("&View")
        for a in (self.a_zoom_in, self.a_zoom_out, self.a_font, self.a_crisp, None, self.a_row_more, self.a_row_less, self.a_col_more, self.a_col_less, None, self.a_pinned_ref, self.a_translation, self.a_features, self.a_shade_feat, self.a_dock, None):
            v.addAction(a) if a else v.addSeparator()
        wm = v.addMenu("Text &weight")
        for a in (self.a_w_reg, self.a_w_semi, self.a_w_bold):
            wm.addAction(a)
        self.scheme_menu = v.addMenu("Color &scheme")
        self.scheme_group = QActionGroup(self)
        self._rebuild_scheme_menu()
        for a in (self.a_col_scheme, self.a_col_ident, self.a_col_none, self.a_normal, self.a_inverse, self.a_dots, None, self.a_ref_cons, self.a_ref_first, self.a_threshold):
            v.addAction(a) if a else v.addSeparator()
        tm = v.addMenu("Amino acid line")
        tm.addAction(self.a_translation)
        self.a_trans_feat = self._act("Use annotated ORFs / CDS (own frame and code)",
                                      self._toggle_trans_features, None, True)
        self.a_trans_feat.setChecked(True)
        tm.addAction(self.a_trans_feat)
        tm.addSeparator()
        self.frame_actions = QActionGroup(self)
        for i, lab in enumerate(["+1", "+2", "+3", "-1", "-2", "-3"]):
            a = self._act(f"Elsewhere: frame {lab}", lambda _, k=i: self._set_frame(k), None, True)
            self.frame_actions.addAction(a); tm.addAction(a)
            if i == 0:
                a.setChecked(True)
        self.a_trans_blank = self._act("Elsewhere: leave blank", self._toggle_trans_fill, None, True)
        self.a_trans_blank.setChecked(True)
        tm.addAction(self.a_trans_blank)
        tm.addSeparator()
        self.a_trans_table = self._act("Genetic code outside ORFs…", self._set_overlay_table)
        tm.addAction(self.a_trans_table)

        s = mb.addMenu("&Sequence")
        self.group_menu = QMenu("&Groups", self)
        self.group_menu.aboutToShow.connect(self._fill_group_menu)
        for a in (self.a_newseq, self.a_rename, self.a_delseq, self.a_dupseq, self.a_up, self.a_down, None,
                  self.a_revcomp, self.a_rev, self.a_comp, self.a_upper, self.a_lower, self.a_rmgaps, self.a_rna, self.a_dna, None,
                  self.a_translate, self.a_sixframe, None, self.a_pin, self.a_settype, self.a_circular, self.a_blast):
            s.addAction(a) if a else s.addSeparator()
        s.addSeparator(); s.addMenu(self.group_menu)

        al = mb.addMenu("&Alignment")
        for a in (self.a_align, None, self.a_insgapcol, self.a_delgapcol, self.a_rm_gapcols, self.a_pad, None, self.a_consensus, self.a_extract):
            al.addAction(a) if a else al.addSeparator()

        an = mb.addMenu("A&nalysis")
        for a in (self.a_orf, self.a_primer, self.a_design, self.a_restrict, None, self.a_stats, self.a_ident, self.a_plot):
            an.addAction(a) if a else an.addSeparator()

        gm = mb.addMenu("&Genome")
        for a in (self.a_g_ref, self.a_g_open, self.a_g_ann, self.a_g_syn, self.a_g_add, None, self.a_g_panel, self.a_g_circ, self.a_g_goto, self.a_g_openreg, self.a_g_orfclear, None, self.a_g_clear):
            gm.addAction(a) if a else gm.addSeparator()

        h = mb.addMenu("&Help")
        h.addAction(self._act("Keyboard shortcuts", self.shortcuts))
        h.addAction(self.a_about)

        # right-click behaviour (BioEdit-style) + context menu
        rc = e.addMenu("&Right-click action")
        self.rc_group = QActionGroup(self)
        self.rc_actions = {}
        for key, label in self.view.RIGHT_CLICK_ACTIONS:
            a = self._act(label, lambda _, k=key: self._set_right_click(k), None, True)
            self.rc_group.addAction(a); rc.addAction(a); self.rc_actions[key] = a
        self.rc_actions[self.view.right_click_action].setChecked(True)
        rc.addSeparator()
        rc.addAction(self._act("(Shift+right-click always shows the context menu)", None))
        self.view.contextMenuWanted.connect(self._context_menu)

    def _set_right_click(self, key):
        self.view.right_click_action = key
        self.settings.setValue("right_click_action", key)
        self.rc_actions[key].setChecked(True)
        if hasattr(self, "rc_combo"):
            self.rc_combo.blockSignals(True)
            self.rc_combo.setCurrentIndex([k for k, _ in self.view.RIGHT_CLICK_ACTIONS].index(key))
            self.rc_combo.blockSignals(False)
        self._update_status()

    def _set_mode(self, mode):
        self.view.set_mode(mode)
        {"slide": self.a_mode_slide, "edit": self.a_mode_edit, "grab": self.a_mode_grab}[mode].setChecked(True)
        edit = mode == "edit"
        self.a_mode_insert.setVisible(edit); self.a_mode_over.setVisible(edit)
        self.settings.setValue("mode", mode)
        self._update_status()

    def _toggle_downstream(self, on):
        self.view.slide_downstream_default = on
        self.settings.setValue("slide_downstream", on)
        self._update_status()

    def _context_menu(self, pos):
        m = QMenu(self)
        m.addAction(self.a_pin)
        m.addMenu(self.group_menu)
        m.addSeparator()
        for a in (self.a_copy, self.a_paste, None, self.a_revcomp, self.a_translate, self.a_rmgaps, None,
                  self.a_insgapcol, self.a_delgapcol, None, self.a_orf, self.a_primer, self.a_blast, None, self.a_rename, self.a_delseq):
            m.addAction(a) if a else m.addSeparator()
        m.exec(self.view.viewport().mapToGlobal(pos))

    def _rebuild_scheme_menu(self):
        self.scheme_menu.clear()
        for a in self.scheme_group.actions():
            self.scheme_group.removeAction(a)
        for name in C.schemes_for(self.model.seq_type):
            a = self._act(name, lambda _, n=name: self.view.set_scheme(n), None, True)
            a.setChecked(name == self.view.scheme_name)
            self.scheme_group.addAction(a)
            self.scheme_menu.addAction(a)

    def _build_toolbar(self):
        tb = QToolBar("Main")
        tb.setIconSize(QSize(22, 22))
        tb.setToolButtonStyle(Qt.ToolButtonIconOnly)
        tb.setMovable(False)
        self.addToolBar(tb)
        for a in (self.a_open, self.a_save, None, self.a_undo, self.a_redo, None,
                  self.a_mode_slide, self.a_mode_edit, self.a_mode_grab, self.a_mode_insert, self.a_mode_over, self.a_downstream, None,
                  self.a_normal, self.a_inverse, None,
                  self.a_zoom_in, self.a_zoom_out, self.a_row_more, self.a_row_less, self.a_col_more, self.a_col_less, None,
                  self.a_translation, self.a_dock, None, self.a_orf, self.a_primer, self.a_align):
            tb.addAction(a) if a else tb.addSeparator()
        for a in tb.actions():
            if a.text() and not a.toolTip():
                a.setToolTip(a.text().replace("&", ""))
        self.a_normal.setToolTip("Normal view: colored letters on neutral background")
        self.a_inverse.setToolTip("Inverse view: white letters on colored background (Ctrl+I)")
        self.addToolBarBreak()
        tb2 = QToolBar("Options")
        tb2.setMovable(False)
        self.addToolBar(tb2)
        tb = tb2
        tb.addWidget(QLabel(" Right-click: "))
        self.rc_combo = QComboBox()
        for key, label in self.view.RIGHT_CLICK_ACTIONS:
            self.rc_combo.addItem(label, key)
        self.rc_combo.currentIndexChanged.connect(lambda i: self._set_right_click(self.rc_combo.itemData(i)))
        tb.addWidget(self.rc_combo)
        tb.addSeparator()
        tb.addWidget(QLabel(" Color: "))
        self.scheme_combo = QComboBox()
        self._fill_scheme_combo()
        self.scheme_combo.currentTextChanged.connect(self._scheme_combo_changed)
        tb.addWidget(self.scheme_combo)

    def _fill_scheme_combo(self):
        self.scheme_combo.blockSignals(True)
        self.scheme_combo.clear()
        for name in C.schemes_for(self.model.seq_type):
            self.scheme_combo.addItem(name)
        self.scheme_combo.setCurrentText(self.view.scheme_name)
        self.scheme_combo.blockSignals(False)

    def _scheme_combo_changed(self, name):
        if name:
            self.view.set_scheme(name)
            for a in self.scheme_group.actions():
                a.setChecked(a.text() == name)

    def _build_status(self):
        sb = self.statusBar()
        self.lbl_pos = QLabel(); self.lbl_sel = QLabel(); self.lbl_mode = QLabel(); self.lbl_info = QLabel()
        for w in (self.lbl_pos, self.lbl_sel, self.lbl_mode):
            sb.addPermanentWidget(w)
        sb.addWidget(self.lbl_info, 1)

    # ------------------------------------------------------------ status
    def _update_status(self, *_):
        m, v = self.model, self.view
        if hasattr(self, "a_circular"):
            self.a_circular.blockSignals(True); self.a_circular.setChecked(bool(m.circular)); self.a_circular.blockSignals(False)
        if m.nrows:
            r, c = v.cur_row, v.cur_col
            seq = m.rows[r].seq
            if c < len(seq):
                upto = c + 1
                ug = upto - seq.count("-", 0, upto) - seq.count(".", 0, upto) - seq.count("~", 0, upto)
            else:
                ug = 0
            ch = seq[c] if c < len(seq) else ""
            self.lbl_pos.setText(f"  {m.rows[r].name}  col {c + 1}  (residue {ug})  [{ch}]  ")
        else:
            self.lbl_pos.setText("")
        s = v.selection()
        if s:
            self.lbl_sel.setText(f"  sel: rows {s[0] + 1}-{s[1] + 1}, cols {s[2] + 1}-{s[3] + 1} ({s[3] - s[2] + 1})  ")
        else:
            self.lbl_sel.setText("")
        self.lbl_mode.setText(f"  {v.mode_text()}  ")
        if v.show_translation:
            frames = ["+1", "+2", "+3", "-1", "-2", "-3"]
            bits = []
            if v.trans_from_features:
                bits.append("ORFs")
            if v.trans_fill == "all" or not v.trans_from_features:
                bits.append(f"frame {frames[v.trans_frame]} code {v.trans_table}")
            self.lbl_mode.setText(f"  {v.mode_text()}  |  aa: {', '.join(bits) or 'off'}  ")
        self.lbl_info.setText(f"{m.nrows} sequences, {m.width} columns, {m.seq_type.upper()}"
                              + (f"  —  {os.path.basename(m.path)}" if m.path else ""))
        self.a_undo.setEnabled(m.can_undo()); self.a_undo.setText(f"&Undo {m.undo_text()}".strip())
        self.a_redo.setEnabled(m.can_redo()); self.a_redo.setText(f"&Redo {m.redo_text()}".strip())
        title = "NeoEdit" + (f" — {os.path.basename(m.path)}" if m.path else " — untitled") + (" *" if m.dirty else "")
        self.setWindowTitle(title)

    def _on_model(self, what):
        self._proj = None
        self._proj_row = -1
        self._ins_dirty = True
        self._ref_ungapped = None
        self._update_status()
        if what == "type":
            self._rebuild_scheme_menu(); self._fill_scheme_combo()

    # ------------------------------------------------------------ file ops
    def _maybe_save(self) -> bool:
        if not self.model.dirty:
            return True
        r = QMessageBox.question(self, "Unsaved changes", "Save changes to the current alignment?",
                                 QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
        if r == QMessageBox.Cancel:
            return False
        if r == QMessageBox.Save:
            return self.save_file()
        return True

    def _set_model(self, model: AlignmentModel):
        self.model.remove_listener(self._on_model)
        self.model = model
        self.model.add_listener(self._on_model)
        self.view.set_model(model)
        self.features_panel.set_model(model)
        self._rebuild_scheme_menu(); self._fill_scheme_combo()
        self._update_status()
        if model.features:
            self.dock.show()

    def new_alignment(self):
        if not self._maybe_save():
            return
        self._set_model(AlignmentModel())

    def open_file(self):
        if not self._maybe_save():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open alignment", self.settings.value("last_dir", ""), FILE_FILTER)
        if path:
            self.open_path(path)

    def open_path(self, path: str):
        try:
            if os.path.getsize(path) > 50_000_000 and path.lower().endswith((".fa", ".fasta", ".fna")):
                r = QMessageBox.question(self, "Large FASTA", "This FASTA is large. Open it as a genome (indexed, one contig at a time)?",
                                         QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
                if r == QMessageBox.Cancel:
                    return
                if r == QMessageBox.Yes:
                    self.open_genome_path(path); return
        except OSError:
            pass
        if self.genome:
            self.genome_close()
        try:
            model = mio.load(path)
        except Exception as e:
            QMessageBox.critical(self, "Open failed", f"{path}\n\n{e}\n\n{traceback.format_exc()[-800:]}")
            return
        self._set_model(model)
        self._reset_reference_mode()
        self.settings.setValue("last_dir", os.path.dirname(path))
        self._add_recent(path)
        msg = f"Opened {path} ({model.nrows} sequences, format {model.format})"
        # An annotated GenBank/EMBL record is a reference: populate the genome view from
        # its own features so gene models (and any ORF tracks) have somewhere to live.
        if model.format in ("genbank", "embl"):
            if self._adopt_annotated_records(path, range(model.nrows)):
                msg += f" — genome view: {self.annotation.count()} annotated features"
        self.statusBar().showMessage(msg, 6000)

    def import_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import sequences", self.settings.value("last_dir", ""), FILE_FILTER)
        if path:
            self.import_path(path)

    def import_path(self, path: str) -> int:
        """Append every sequence in `path` to the current alignment; returns the number added."""
        try:
            m = mio.load(path)
        except Exception as e:
            QMessageBox.critical(self, "Import failed", str(e)); return 0
        offset = self.model.nrows
        self.model.begin_batch("Import")
        for r in m.rows:
            self.model.add_row(r)
        self.model.end_batch()
        # keep the GenBank/EMBL features of the imported records, re-homed onto their new rows
        for f in m.features:
            f.row += offset
            self.model.features.append(f)
        if m.features:
            self.model._emit("features")
        msg = f"Imported {m.nrows} sequence(s) from {path}"
        # An annotated record arrives with a chromosome map. The map and the gene models follow
        # the reference row, so the record is pinned as the reference to put its features on screen.
        if m.format in ("genbank", "embl"):
            if self._adopt_annotated_records(path, range(offset, offset + m.nrows)):
                msg += (f" — pinned {self.model.rows[self.model.ref_row].name} as the reference; "
                        f"genome view: {self.annotation.count()} annotated features")
        self.statusBar().showMessage(msg, 5000)
        return m.nrows

    def import_remote(self):
        """Modeless NCBI / Ensembl importer; every record is added to the current alignment via import_path."""
        dlg = getattr(self, "_import_dialog", None)
        if dlg is None:
            dlg = ImportDialog(self, self.settings, self.import_path)
            self._import_dialog = dlg
        dlg.show(); dlg.raise_(); dlg.activateWindow()

    def save_file(self) -> bool:
        if not self.model.path or self.model.format == "genbank" and self.model.nrows > 1:
            return self.save_file_as()
        try:
            mio.save(self.model, self.model.path, self.model.format)
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e)); return False
        self._update_status()
        self.statusBar().showMessage(f"Saved {self.model.path}", 3000)
        return True

    def save_file_as(self) -> bool:
        filters = [f"{lbl} (" + " ".join("*" + e for e in exts) + ")" for lbl, _, exts, _ in mio.FORMATS]
        start = self.model.path or os.path.join(self.settings.value("last_dir", ""), "alignment.fasta")
        path, chosen = QFileDialog.getSaveFileName(self, "Save alignment as", start, ";;".join(filters))
        if not path:
            return False
        idx = filters.index(chosen) if chosen in filters else 0
        fmt = mio.FORMATS[idx][1]
        exts = mio.FORMATS[idx][2]
        if not os.path.splitext(path)[1]:
            path += exts[0]
        try:
            mio.save(self.model, path, fmt)
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e)); return False
        self.settings.setValue("last_dir", os.path.dirname(path))
        self._add_recent(path)
        self._update_status()
        return True

    def export_selected(self):
        rows = self.view.target_rows()
        if not rows:
            return
        sub = AlignmentModel([self.model.rows[i].copy() for i in rows], self.model._seq_type)
        old, self.model = self.model, sub
        try:
            self.save_file_as()
        finally:
            self.model = old
            self._update_status()

    def _recent(self) -> list[str]:
        v = self.settings.value("recent", [])
        if isinstance(v, str):
            v = [v]
        return [p for p in (v or []) if isinstance(p, str) and os.path.sep in p]

    def _add_recent(self, path):
        path = os.path.abspath(path)
        rec = [p for p in self._recent() if p != path]
        rec.insert(0, path)
        self.settings.setValue("recent", rec[:10])
        self._fill_recent()

    def _fill_recent(self):
        self.recent_menu.clear()
        for p in self._recent():
            self.recent_menu.addAction(self._act(p, lambda _, q=p: (self._maybe_save() and self.open_path(q))))

    def _restore(self):
        self._fill_recent()
        g = self.settings.value("geometry")
        if g:
            self.restoreGeometry(g)
        fs = self.settings.value("font_size")
        if fs:
            self.view.font_size = int(fs); self.view._apply_font()
        fam = self.settings.value("font_family")
        if fam:
            self.view.set_font_family(fam)
        md = self.settings.value("mode", "slide")
        self._set_mode(md if md in self.view.MODES else "slide")
        sd = self.settings.value("slide_downstream", False)
        sd = sd in (True, "true", "True", 1, "1")
        self.a_downstream.setChecked(sd); self.view.slide_downstream_default = sd
        rck = self.settings.value("right_click_action")
        if rck in dict(self.view.RIGHT_CLICK_ACTIONS):
            self._set_right_click(rck)
        tf = self.settings.value("trans_frame")
        if tf is not None:
            self.view.trans_frame = int(tf)
        tt = self.settings.value("trans_table")
        if tt is not None:
            self.view.trans_table = int(tt)
        tff = self.settings.value("trans_from_features")
        if tff is not None:
            on = tff in (True, "true", "True", 1, "1")
            self.a_trans_feat.setChecked(on); self.view.trans_from_features = on
        tfl = self.settings.value("trans_fill")
        if tfl:
            self.view.trans_fill = tfl
            self.a_trans_blank.setChecked(tfl == "regions")
        tw = self.settings.value("text_weight")
        if tw in ("regular", "semibold", "bold"):
            {"regular": self.a_w_reg, "semibold": self.a_w_semi, "bold": self.a_w_bold}[tw].setChecked(True)
            self.view.set_text_style(weight=tw)
        cr = self.settings.value("crisp_text")
        if cr is not None:
            cr = cr in (True, "true", "True", 1, "1")
            self.a_crisp.setChecked(cr); self.view.set_text_style(crisp=cr)
        rp, cp = self.settings.value("row_pad"), self.settings.value("col_pad")
        if rp is not None or cp is not None:
            self.view.set_spacing(int(rp) if rp is not None else None, int(cp) if cp is not None else None)
        sh = self.settings.value("shade_features", False)
        sh = sh in (True, "true", "True", 1, "1")
        self.a_shade_feat.setChecked(sh); self.view.shade_features = sh
        inv = self.settings.value("inverse_view", False)
        inv = inv in (True, "true", "True", 1, "1")
        self._set_inverse(inv)
        sch = self.settings.value("scheme")
        if sch and sch in C.schemes_for(self.model.seq_type):
            self.view.set_scheme(sch); self._fill_scheme_combo()

    def closeEvent(self, e):
        if not self._maybe_save():
            e.ignore(); return
        self.settings.setValue("geometry", self.saveGeometry())
        self.settings.setValue("font_size", self.view.font_size)
        self.settings.setValue("font_family", self.view.font_family)
        self.settings.setValue("scheme", self.view.scheme_name)
        self.settings.setValue("text_weight", self.view.text_weight)
        self.settings.setValue("crisp_text", self.view.crisp_text)
        self.settings.setValue("row_pad", self.view.row_pad)
        self.settings.setValue("col_pad", self.view.col_pad)
        # close any non-modal dialogs / secondary windows we opened
        for w in list(self._children):
            try:
                w.close()
            except RuntimeError:
                pass
        self._children.clear()
        if self.find_dlg is not None:
            self.find_dlg.close()
        dlg = getattr(self, "_import_dialog", None)
        if dlg is not None:
            dlg.close()
        e.accept()

    # ------------------------------------------------------------ edit ops
    def undo(self):
        self.model.undo()

    def redo(self):
        self.model.redo()

    def copy_sel(self):
        txt = self.view.selected_text()
        if not txt and self.model.nrows:
            r = self.model.rows[self.view.cur_row]
            txt = f">{r.name}\n{r.seq}\n"
        QApplication.clipboard().setText(txt)

    def copy_raw(self):
        s = self.view.selection()
        if s:
            r0, r1, c0, c1 = s
            QApplication.clipboard().setText("\n".join(self.model.rows[i].seq[c0:c1 + 1] for i in range(r0, r1 + 1)))
        elif self.model.nrows:
            QApplication.clipboard().setText(self.model.rows[self.view.cur_row].seq)

    def paste_seqs(self):
        txt = QApplication.clipboard().text().strip()
        if not txt:
            return
        import re
        if txt.startswith(">"):
            try:
                m = mio.loads(txt, "fasta")
            except Exception as e:
                QMessageBox.warning(self, "Paste", f"Could not parse FASTA: {e}"); return
            rows = m.rows
        else:
            seq = re.sub(r"[\s\d]", "", txt)
            if self.view.typing is not None and self.model.nrows:
                # paste into current row at cursor
                if self.view.typing == "insert":
                    self.model.insert_text(self.view.cur_row, self.view.cur_col, seq)
                else:
                    self.model.overwrite(self.view.cur_row, self.view.cur_col, seq)
                self.view.set_cursor(self.view.cur_row, self.view.cur_col + len(seq))
                return
            rows = [SequenceRow(f"pasted_{self.model.nrows + 1}", seq)]
        self.model.begin_batch("Paste")
        for r in rows:
            self.model.add_row(r)
        self.model.end_batch()

    def find(self):
        if self.find_dlg is None:
            self.find_dlg = FindDialog(self)
            self.find_dlg.findNext.connect(self._do_find)
            self.find_dlg.findInNames.connect(self._find_name)
        self.find_dlg.show(); self.find_dlg.raise_(); self.find_dlg.pattern.setFocus()

    def find_next(self):
        if self.find_dlg and self.find_dlg.pattern.text():
            d = self.find_dlg
            self._do_find(d.pattern.text(), d.ignore_gaps.isChecked(), d.regex.isChecked(), d.case.isChecked())
        else:
            self.find()

    def _do_find(self, pattern, ignore_gaps, regex, case):
        if not pattern or not self.model.nrows:
            return
        try:
            hit = self.model.find(pattern, self.view.cur_row, self.view.cur_col + 1, ignore_gaps, regex, case)
        except Exception as e:
            self.statusBar().showMessage(f"Bad pattern: {e}", 4000); return
        if hit:
            r, s, e = hit
            self.view.select_region(r, r, s, e - 1)
            self.statusBar().showMessage(f"Found at {self.model.rows[r].name} col {s + 1}", 3000)
        else:
            self.statusBar().showMessage("Not found", 3000)

    def _find_name(self, text):
        t = text.lower()
        start = self.view.cur_row + 1
        order = list(range(start, self.model.nrows)) + list(range(0, start))
        for r in order:
            if t in self.model.rows[r].name.lower():
                self.view.sel_rows = {r}; self.view.anchor = None
                self.view.set_cursor(r, self.view.cur_col)
                self.view.sel_rows = {r}; self.view.viewport().update()
                return
        self.statusBar().showMessage("Name not found", 3000)

    def goto(self):
        if not self.model.nrows:
            return
        n, ok = QInputDialog.getInt(self, "Go to", "Column (alignment position):", self.view.cur_col + 1, 1, max(1, self.model.width))
        if ok:
            self.view.set_cursor(self.view.cur_row, n - 1)

    def preferences(self):
        PreferencesDialog(self.settings, self).exec()

    # ------------------------------------------------------------ view ops
    def _toggle_pinned_ref(self, on):
        self.view.show_reference = on; self.view._refresh()

    def _toggle_translation(self, on):
        self.view.show_translation = on; self.view._refresh()

    def _toggle_features(self, on):
        self.view.show_features = on; self.view.viewport().update()

    def _set_inverse(self, on: bool):
        (self.a_inverse if on else self.a_normal).setChecked(True)
        self.view.color_target = "background" if on else "text"
        self.settings.setValue("inverse_view", on)
        self.view.viewport().update()

    def _toggle_shade_features(self, on):
        self.view.shade_features = on
        self.settings.setValue("shade_features", on)
        self.view.viewport().update()

    def _toggle_dots(self, on):
        self.view.dots_for_identity = on; self.view.viewport().update()

    def _set_color_mode(self, mode):
        self.view.color_mode = mode; self.view.viewport().update()

    def _set_ref(self, ref):
        self.view.identity_ref = ref; self.view.viewport().update()

    def _set_threshold(self):
        v, ok = QInputDialog.getDouble(self, "Consensus threshold", "Fraction of residues that must agree:",
                                       self.view.shade_threshold, 0.0, 1.0, 2)
        if ok:
            self.view.shade_threshold = v; self.view._consensus_cache = None; self.view.viewport().update()

    def choose_font(self):
        from PySide6.QtWidgets import QFontDialog
        from PySide6.QtGui import QFont
        cur = QFont(self.view.font_family, self.view.font_size)
        font, ok = QFontDialog.getFont(cur, self, "Alignment font", QFontDialog.MonospacedFonts)
        if ok:
            self.view.font_size = font.pointSize()
            self.view.set_font_family(font.family())

    def _set_frame(self, k):
        self.view.trans_frame = k
        self.view.trans_fill = "all"
        if hasattr(self, "a_trans_blank"):
            self.a_trans_blank.setChecked(False)
        self.settings.setValue("trans_frame", k)
        self.view.viewport().update(); self._update_status()

    def _toggle_trans_features(self, on):
        self.view.trans_from_features = on
        self.settings.setValue("trans_from_features", on)
        self.view.viewport().update(); self._update_status()

    def _toggle_trans_fill(self, on):
        self.view.trans_fill = "regions" if on else "all"
        self.settings.setValue("trans_fill", self.view.trans_fill)
        self.view.viewport().update(); self._update_status()

    def _set_overlay_table(self):
        tables = T.codon_tables()
        items = [f"{i}: {n}" for i, n in tables]
        cur = next((k for k, (i, _) in enumerate(tables) if i == self.view.trans_table), 0)
        s, ok = QInputDialog.getItem(self, "Genetic code", "Table:", items, cur, False)
        if ok:
            self.view.trans_table = int(s.split(":")[0])
            self.settings.setValue("trans_table", self.view.trans_table)
            self.view.viewport().update(); self._update_status()

    # ------------------------------------------------------------ sequence ops
    def rename(self):
        if not self.model.nrows:
            return
        r = self.view.cur_row
        name, ok = QInputDialog.getText(self, "Rename", "Name:", text=self.model.rows[r].name)
        if ok and name:
            self.model.rename(r, name)

    def new_sequence(self):
        d = NewSequenceDialog(self)
        if d.exec():
            self.model.begin_batch("New sequence")
            for name, seq in d.records():
                if seq:
                    self.model.add_row(SequenceRow(name, seq))
            self.model.end_batch()

    def delete_seqs(self):
        rows = self.view.target_rows()
        if not rows:
            return
        if len(rows) > 1 and QMessageBox.question(self, "Delete", f"Delete {len(rows)} sequences?") != QMessageBox.Yes:
            return
        self.model.remove_rows(rows)
        self.view.clear_selection()

    def _move(self, d):
        rows = self.view.target_rows()
        self.model.move_rows(rows, d)
        if self.view.sel_rows:
            self.view.sel_rows = {r + d for r in rows}
        self.view.set_cursor(self.view.cur_row + d, self.view.cur_col)

    def group_assign(self):
        rows = self.view.target_rows()
        if not rows:
            return
        existing = self.model.groups()
        name, ok = QInputDialog.getItem(self, "Assign to group", "Group name (pick or type a new one):",
                                        existing + ["include", "exclude"] if not existing else existing, 0, True)
        if ok and name.strip():
            self.model.set_group(rows, name.strip())

    def group_color(self):
        groups = self.model.groups()
        if not groups:
            return
        name, ok = QInputDialog.getItem(self, "Group color", "Group:", groups, 0, False)
        if not ok:
            return
        from PySide6.QtWidgets import QColorDialog
        from PySide6.QtGui import QColor
        c = QColorDialog.getColor(QColor(self.model.group_color(name)), self, f"Color for '{name}'")
        if c.isValid():
            self.model.set_group_color(name, c.name())

    def _fill_group_menu(self):
        m = self.group_menu
        m.clear()
        m.addAction(self.a_grp_set)
        m.addAction(self.a_grp_clear)
        m.addAction(self.a_grp_color)
        groups = self.model.groups()
        if groups:
            m.addSeparator()
            for g in groups:
                n = len(self.model.group_rows(g))
                sub = m.addMenu(f"{g}  ({n})")
                sub.addAction(self._act("Select rows", lambda _, gg=g: self._group_select(gg)))
                sub.addAction(self._act("Assign selected rows here", lambda _, gg=g: self.model.set_group(self.view.target_rows(), gg)))
                sub.addAction(self._act("Dissolve group", lambda _, gg=g: self.model.set_group(self.model.group_rows(gg), "")))

    def _group_select(self, name):
        rows = self.model.group_rows(name)
        if rows:
            self.view.sel_rows = set(rows)
            self.view.anchor = None
            self.view.cur_row = rows[0]
            self.view.selectionChanged.emit()
            self.view.viewport().update()

    def pin_reference(self):
        if not self.model.nrows:
            return
        row = self.view.cur_row
        self.model.set_reference(row)
        self._proj = None; self._proj_row = -1; self._ref_ungapped = None
        name = self.model.rows[row].name
        msg = f"Pinned {name} as the reference sequence"
        # Features belong to the reference: re-point the map and the gene models at this row.
        if self._annotated(name) or self.genome_panel.isVisible() or self.genome_contig:
            self._enter_reference_mode(self.annotation)
            n = len(self.annotation.genes_by_seq.get(name, [])) if self.annotation else 0
            msg += f" — {n} annotated features" if n else " — no annotation for this sequence"
        self.view.viewport().update()
        self.statusBar().showMessage(msg, 5000)

    def _toggle_circular(self, on):
        self.model.circular = on
        self.view.viewport().update()

    def set_type(self):
        s, ok = QInputDialog.getItem(self, "Sequence type", "Type:", ["auto", "dna", "rna", "protein"], 0, False)
        if ok:
            self.model.seq_type = None if s == "auto" else s

    def blast(self):
        if not self.model.nrows:
            return
        seq = self.model.rows[self.view.cur_row].seq
        s = self.view.selection()
        if s and s[0] == s[1]:
            seq = seq[s[2]:s[3] + 1]
        prog = "blastn" if self.model.is_nucleotide() else "blastp"
        EXT.open_blast(seq, prog)

    def translate(self):
        rows = self.view.target_rows()
        if not rows:
            return
        d = TranslateDialog(self, int(self.settings.value("default_table", 1)))
        if not d.exec():
            return
        v = d.values()
        out_rows = []
        from Bio.Seq import Seq
        for r in rows:
            row = self.model.rows[r]
            frames = [v["frame"]] if v["frame"] < 6 else list(range(6))
            for fr in frames:
                seq = row.seq
                if fr >= 3:
                    seq = str(Seq(seq).reverse_complement())
                f = fr % 3
                if v["keep_align"] and not v["to_stop"] and v["frame"] < 6:
                    aa = T.translate_aligned(seq, v["table"], f)
                else:
                    aa = T.translate_gapped(seq, v["table"], f, v["to_stop"])
                suffix = "" if len(frames) == 1 else f"_{'+' if fr < 3 else '-'}{f + 1}"
                out_rows.append(SequenceRow(row.name + suffix, aa, row.description))
        if v["output"] == 0:
            w = MainWindow()
            w._set_model(AlignmentModel(out_rows, "protein"))
            w.show(); self._children.append(w)
        elif v["output"] == 1:
            txt = "\n".join(f">{r.name}\n{r.seq}" for r in out_rows)
            TextDialog(self, "Translation", txt).exec()
        else:
            self.model.begin_batch("Translate")
            for r, nr in zip(rows, out_rows):
                self.model.set_sequence(r, nr.seq, "Translate")
            self.model.end_batch()
            self.model.seq_type = "protein"

    def six_frame(self):
        if not self.model.nrows:
            return
        row = self.model.rows[self.view.cur_row]
        table = int(self.settings.value("default_table", 1))
        sf = T.six_frame(row.seq, table)
        txt = f"Six-frame translation of {row.name} (table {table})\n\n" + "\n\n".join(f"Frame {k}:\n{v}" for k, v in sf.items())
        TextDialog(self, "Six-frame translation", txt).exec()

    # ------------------------------------------------------------ alignment ops
    def align_external(self):
        if self.model.nrows < 2:
            QMessageBox.information(self, "Align", "Need at least two sequences."); return
        d = AlignDialog(self, self.settings)
        if not d.exec():
            return
        v = d.values()
        rows = self.view.target_rows() if v.pop("sel_only") else list(range(self.model.nrows))
        if len(rows) < 2:
            rows = list(range(self.model.nrows))
        exe = self.settings.value("exe/MAFFT") or None
        prog = QProgressDialog("Running MAFFT…", None, 0, 0, self)
        prog.setWindowModality(Qt.WindowModal); prog.show(); QApplication.processEvents()
        try:
            aligned = EXT.run_mafft([self.model.rows[i] for i in rows], exe, seq_type=self.model.seq_type, **v)
        except Exception as e:
            prog.close()
            QMessageBox.critical(self, "MAFFT", str(e)); return
        prog.close()
        self.model.begin_batch("Align (MAFFT)")
        for i, r in zip(rows, aligned):
            self.model.set_sequence(i, r.seq)
            if r.description != self.model.rows[i].description:
                self.model.rows[i].description = r.description
        self.model.end_batch()
        self.statusBar().showMessage(f"MAFFT finished: {len(rows)} sequences aligned", 5000)

    def extract_cols(self):
        s = self.view.selection()
        if not s:
            return
        sub = self.model.column_slice(s[2], s[3] + 1)
        sub.rows = [sub.rows[i] for i in range(s[0], s[1] + 1)]
        w = MainWindow(); w._set_model(sub); w.show(); self._children.append(w)

    # ------------------------------------------------------------ analysis
    def orf_finder(self):
        if not self.model.nrows:
            return
        rows = self.view.target_rows()
        d = ORFFinderDialog(self.model, rows, self, int(self.settings.value("default_table", 1)))
        d.orfSelected.connect(lambda r, s, e: self.view.select_region(r, r, s, e - 1))
        d.featuresReady.connect(self._add_features)
        d.orfsFound.connect(self._orfs_to_track)
        d.show(); self._children.append(d)

    def primer_design(self):
        if not self.model.nrows:
            return
        if not self.model.is_nucleotide():
            QMessageBox.information(self, "Primer design", "Primer design needs nucleotide sequences."); return
        rows = self.view.target_rows()
        s = self.view.selection()
        cols = (s[2], s[3]) if s and (s[2] != s[3]) else None
        d = PrimerDialog(self.model, rows, cols, self)
        d.pairSelected.connect(lambda r, s_, e: self.view.select_region(r, r, s_, e - 1))
        d.featuresReady.connect(self._add_features)
        d.show(); self._children.append(d)

    def restriction_sites(self):
        if not self.model.nrows:
            return
        if not self.model.is_nucleotide():
            QMessageBox.information(self, "Restriction sites", "Needs nucleotide sequences."); return
        d = RestrictionDialog(self.model, self.view, self)
        d.featuresReady.connect(self._add_features)
        d.siteSelected.connect(lambda r, a, b: self.view.select_region(r, r, a, b - 1))
        d.show(); self._children.append(d)

    def design_primers(self):
        if self.model.nrows < 2:
            QMessageBox.information(self, "Design", "Load an alignment of at least two sequences first."); return
        if not self.model.is_nucleotide():
            QMessageBox.information(self, "Design", "Primer design needs nucleotide sequences."); return
        d = DesignDialog(self.model, self.view, self)
        d.featuresReady.connect(self._add_features)
        d.regionSelected.connect(lambda a, b: self.view.select_region(0, self.model.nrows - 1, a, b - 1))
        d.show(); self._children.append(d)

    ORF_TRACK_COLORS = {1: "#7d3c98", 2: "#0e7490", 3: "#b45309", 4: "#166534", 5: "#9d174d"}
    # distinct colors for ORFs named after a known peptide (stable per name)
    NAMED_ORF_COLORS = ["#dc2626", "#ea580c", "#ca8a04", "#16a34a", "#0891b2", "#2563eb", "#7c3aed", "#db2777"]

    def _named_color(self, name: str) -> str:
        h = 0
        for ch in name.lower():
            h = (h * 31 + ord(ch)) & 0xFFFFFFFF
        return self.NAMED_ORF_COLORS[h % len(self.NAMED_ORF_COLORS)]

    def _orfs_to_track(self, orfs, table: int):
        """Show ORF-finder results as their own track under the gene models."""
        if not self.genome_panel.isVisible() and not orfs:
            return
        ref_row = self.ref_index()
        proj = self.proj()
        items = []
        for o in orfs:
            if o.row != ref_row:
                continue
            s0 = proj.col_to_ref(o.start) if not proj.identity else o.start
            e0 = proj.col_to_ref(o.end) if not proj.identity else o.end
            partial = ("5'" if o.partial5 else "") + ("3'" if o.partial3 else "")
            tip = (f"<b>{o.name or 'ORF'}</b> {s0 + 1:,}-{e0:,} ({'+' if o.strand > 0 else '-'}{o.frame})<br>"
                   f"{o.length_aa} aa, start {o.start_codon}, stop {o.stop_codon or 'none'}"
                   f"{', partial ' + partial if partial else ''}<br>genetic code {o.table}"
                   + (f"<br>{o.extra['note']}" if o.extra.get('note') else ""))
            matched = "similar to" in (o.extra.get("note") or "")
            color = self._named_color(o.name) if (matched and o.name) else None
            items.append((s0, e0, o.strand, o.name or f"{o.length_aa}aa", tip, color))
        if not items:
            return
        label = f"code {table}"
        self.genome_panel.add_orf_track(label, items, self.ORF_TRACK_COLORS.get(table, "#7d3c98"))
        if not self.genome_panel.isVisible() and self.model.nrows:
            if self.genome_contig is None:
                self._enter_reference_mode(self.annotation)
            else:
                self.a_g_panel.setChecked(True); self._toggle_genome_panel(True)
            self.genome_panel.add_orf_track(label, items, self.ORF_TRACK_COLORS.get(table, "#7d3c98"))
        self.statusBar().showMessage(f"ORF track '{label}': {len(items)} ORFs shown under the gene models", 5000)

    def _add_features(self, feats):
        self.model.features.extend(feats)
        self.features_panel.refresh()
        if not self.dock.isVisible():
            self.dock.show()
            self.resizeDocks([self.dock], [170], Qt.Vertical)
        self.view.viewport().update()

    def _goto_feature(self, f: Feature):
        self.view.select_region(f.row, f.row, f.start, f.end - 1)

    def stats(self):
        if self.model.nrows:
            StatsDialog(self.model, self.view.target_rows() if len(self.view.target_rows()) > 1 else range(self.model.nrows), self).exec()

    def identity(self):
        if self.model.nrows:
            rows = self.view.target_rows()
            IdentityDialog(self.model, rows if len(rows) > 1 else range(self.model.nrows), self).exec()

    def plot(self):
        if self.model.nrows:
            rows = self.view.target_rows()
            PlotDialog(self.model, rows if len(rows) > 1 else range(self.model.nrows), self).exec()

    def consensus_tool(self):
        if not self.model.nrows:
            QMessageBox.information(self, "Consensus", "Open an alignment first."); return
        ConsensusDialog(self.model, self.view.target_rows(), self.view.shade_threshold, self).exec()

    # ------------------------------------------------------------ genome
    def _toggle_genome_panel(self, on):
        self.genome_panel.setVisible(on)
        if on:
            self.splitter.setSizes([260, max(200, self.height() - 260)])

    def ref_index(self) -> int:
        r = self.model.ref_row
        return r if 0 <= r < self.model.nrows else 0

    def ref_seq(self) -> str:
        return self.model.rows[self.ref_index()].seq if self.model.nrows else ""

    def proj(self) -> RefProjection:
        seq = self.ref_seq()
        if self._proj is None or self._proj.ncols != len(seq) or self._proj_row != self.ref_index():
            self._proj = RefProjection(seq)
            self._proj_row = self.ref_index()
            self._ref_ungapped = None
        return self._proj

    def ref_ungapped(self) -> str:
        if self._ref_ungapped is None:
            self._ref_ungapped = self.model.rows[self.ref_index()].ungapped() if self.model.nrows else ""
        return self._ref_ungapped

    def insertion_carets(self) -> list[tuple[int, int, int, int]]:
        """Insertions the tracks have but the reference lacks.

        Returns (ref_position, n_columns, n_sequences, first_column) for each run of
        columns where the reference is gapped but other sequences carry residues.
        """
        m = self.model
        if m.nrows < 2:
            return []
        proj = self.proj()
        gaps = proj.gap_cols
        if not gaps:
            return []
        runs = []
        start = prev = gaps[0]
        for c in gaps[1:]:
            if c != prev + 1:
                runs.append((start, prev + 1)); start = c
            prev = c
        runs.append((start, prev + 1))
        out = []
        ri = self.ref_index()
        for a, b in runs:
            n = 0
            for i, r in enumerate(m.rows):
                if i == ri:
                    continue
                seg = r.seq[a:b]
                if any(ch not in "-.~" for ch in seg):
                    n += 1
            if n:
                out.append((proj.col_to_ref(a), b - a, n, a))
        return out

    def genome_open(self):
        if not self._maybe_save():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open genome FASTA", self.settings.value("last_genome_dir", self.settings.value("last_dir", "")),
                                              "FASTA (*.fa *.fasta *.fna *.fa.gz);;All files (*)")
        if not path:
            return
        self.open_genome_path(path)

    def open_genome_path(self, path: str, contig: str | None = None):
        if path.endswith(".gz"):
            QMessageBox.information(self, "Genome", "Please decompress the FASTA first (random access needs an uncompressed file or bgzip+index)."); return
        try:
            prog = QProgressDialog("Indexing genome… (first time only)", None, 0, 0, self)
            prog.setWindowModality(Qt.WindowModal); prog.show(); QApplication.processEvents()
            genome = IndexedFasta(path)
            prog.close()
        except Exception as e:
            QMessageBox.critical(self, "Genome", f"Could not open {path}:\n{e}"); return
        if self.genome:
            self.genome.close()
        self.genome = genome
        self.annotation = None; self.synteny_blocks = []
        self.settings.setValue("last_genome_dir", os.path.dirname(path))
        contigs = genome.contigs()
        self.genome_panel.set_contigs(contigs)
        self.genome_panel.set_annotation(None)
        self.a_g_panel.setChecked(True); self._toggle_genome_panel(True)
        # pick the contig: largest by default, let the user choose via combo; load the first (largest) now
        if contigs:
            first = contig if contig in genome.records else max(contigs, key=lambda c: c[1])[0]
            self.genome_panel.set_contigs(contigs, current=first)
            self._genome_load_contig(first)
        self._add_recent(path)

    def load_annotation_path(self, path: str, only: str | None = None):
        ann = GA.load_annotation(path, only)
        self.annotation = ann
        self.genome_panel.set_annotation(ann)
        self.view.viewport().update()
        return ann

    def _genome_load_contig(self, seqid: str):
        if not self.genome or seqid not in self.genome.records:
            return
        if self.model.dirty and not self._maybe_save():
            return
        L = self.genome.length(seqid)
        prog = QProgressDialog(f"Loading {seqid} ({L / 1e6:.1f} Mb)…", None, 0, 0, self)
        prog.setWindowModality(Qt.WindowModal); prog.show(); QApplication.processEvents()
        seq = self.genome.fetch(seqid, 0, L)
        prog.close()
        m = AlignmentModel([SequenceRow(seqid, seq, f"{os.path.basename(self.genome.path)}")], "dna")
        m.dirty = False
        self.genome_contig = seqid
        self._set_model(m)
        self.view.feature_provider = self._genome_features
        self.view.translation_provider = self._translation_regions
        self.genome_panel.set_contig(seqid, L, fetch_seq=lambda s, e, _id=seqid: self.genome.fetch(_id, s, e))
        self.genome_panel.region.fetch_var = self._variation
        self.genome_panel.set_annotation(self.annotation)
        self.genome_panel.set_synteny(self.synteny_blocks)
        self._grid_scrolled()
        self.statusBar().showMessage(f"Loaded {seqid}: {L:,} bp", 5000)

    def _enter_reference_mode(self, ann: GA.Annotation | None):
        """Use row 0 of the current model as the reference for the genome panel."""
        if not self.model.nrows:
            return
        seqid = self.model.rows[self.ref_index()].name
        self.genome_contig = seqid
        self.annotation = ann
        L = self.proj().ref_len
        self.genome_panel.set_contigs([(seqid, L)], current=seqid)
        self.genome_panel.set_contig(seqid, L, fetch_seq=lambda s, e: self.ref_ungapped()[s:e])
        self.genome_panel.region.fetch_var = self._variation
        self.genome_panel.set_annotation(ann)
        self.genome_panel.set_insertions(self.insertion_carets())
        self.view.feature_provider = self._genome_features
        self.view.translation_provider = self._translation_regions
        self.a_g_panel.setChecked(True); self._toggle_genome_panel(True)
        self.genome_panel.set_window(0, min(L, max(2000, L)))
        self._grid_scrolled()

    # -------------------------------------------------- annotated references
    def _annotated(self, name: str) -> bool:
        """Does the window's annotation carry genes for this sequence?"""
        return bool(self.annotation and self.annotation.genes_by_seq.get(name))

    def _reset_reference_mode(self):
        """Forget the annotation belonging to a model that is being replaced."""
        if self.genome:
            return                      # genome mode owns its own annotation
        self.annotation = None
        self.genome_contig = None
        self.view.feature_provider = None
        self.a_g_panel.setChecked(False)
        self.genome_panel.hide()

    def _merge_annotation(self, ann: GA.Annotation) -> GA.Annotation:
        """Fold `ann` into the window's annotation. An Annotation is keyed by seqid, so several
        annotated records coexist happily and the reference row's name selects which is shown."""
        cur = self.annotation
        if cur is None:
            ann.finalize()
            return ann
        for genes in ann.genes_by_seq.values():
            for g in genes:
                cur.add_gene(g)
        cur.finalize()
        return cur

    def _adopt_annotated_records(self, path: str, rows) -> int:
        """Pin an annotated record from `path` as the reference and show its chromosome map.

        Annotations are a property of the reference sequence: the map and the gene models follow
        whichever row is pinned, so a record has to become the reference for its features to be on
        screen. Returns how many annotated records were found among `rows`.
        """
        try:
            ann = GA.load_annotation(path)
        except Exception:
            return 0
        if ann is None or not ann.count():
            return 0
        seqids = set(ann.seqids())
        adopted = [r for r in rows if 0 <= r < self.model.nrows and self.model.rows[r].name in seqids]
        if not adopted:
            return 0
        # the annotation drives the display from here; drop the duplicate per-row copies
        drop = set(adopted)
        kept = [f for f in self.model.features if f.row not in drop]
        if len(kept) != len(self.model.features):
            self.model.features = kept
            self.model._emit("features")
        self.annotation = self._merge_annotation(ann)
        self.model.set_reference(adopted[0])
        self._proj = None; self._proj_row = -1; self._ref_ungapped = None
        self._enter_reference_mode(self.annotation)
        return len(adopted)

    def genome_open_reference(self):
        if not self._maybe_save():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open GenBank reference", self.settings.value("last_genome_dir", self.settings.value("last_dir", "")),
                                              "GenBank (*.gb *.gbk *.genbank);;All files (*)")
        if not path:
            return
        try:
            model = mio.load(path, "genbank")
            ann = GA.load_genbank(path)
        except Exception as e:
            QMessageBox.critical(self, "GenBank", str(e)); return
        model.features = []           # gene view shows these; keep grid overlay via provider
        model.dirty = False
        if self.genome:
            self.genome_close()
        self._set_model(model)
        self.settings.setValue("last_genome_dir", os.path.dirname(path))
        self._enter_reference_mode(ann)
        self._add_recent(path)
        self.statusBar().showMessage(f"Reference {model.rows[0].name}: {self.proj().ref_len:,} bp, {ann.count()} features", 6000)

    def genome_add_anchored(self):
        if not self.model.nrows:
            QMessageBox.information(self, "Anchor", "Open a reference first."); return
        path, _ = QFileDialog.getOpenFileName(self, "Add sequences anchored to reference", self.settings.value("last_dir", ""), FILE_FILTER)
        if not path:
            return
        try:
            new = mio.load(path)
        except Exception as e:
            QMessageBox.critical(self, "Open failed", str(e)); return
        if not new.rows:
            return
        prog = QProgressDialog(f"MAFFT --add: anchoring {len(new.rows)} sequence(s)…", None, 0, 0, self)
        prog.setWindowModality(Qt.WindowModal); prog.show(); QApplication.processEvents()
        try:
            rows = EXT.mafft_add(self.model.rows, new.rows, self.settings.value("exe/MAFFT") or None)
        except Exception as e:
            prog.close(); QMessageBox.critical(self, "MAFFT", str(e)); return
        prog.close()
        self.model.begin_batch("Anchor sequences")
        for r in rows[self.model.nrows:]:
            self.model.add_row(r)
        self.model.end_batch()
        self.statusBar().showMessage(f"Anchored {len(new.rows)} sequence(s) to {self.model.rows[0].name} (columns preserved)", 6000)

    def _variation(self, s: int, e: int):
        """Per-reference-position variant frequency across rows 2..N (None if not applicable)."""
        m = self.model
        if m.nrows < 2 or e - s > 200_000 or e <= s:
            return None
        proj = self.proj()
        c0, c1 = proj.span_to_cols(s, e)
        ref = self.ref_seq()
        out = [0.0] * (e - s)
        u = s
        for c in range(c0, min(c1, len(ref))):
            rc = ref[c]
            if rc in "-.~":
                continue
            rc = rc.upper()
            tot = mis = 0
            for r in (row for i, row in enumerate(m.rows) if i != self.ref_index()):
                ch = r.seq[c].upper() if c < len(r.seq) else ""
                if ch and ch not in "-.~N":
                    tot += 1
                    if ch != rc:
                        mis += 1
            if tot:
                out[u - s] = mis / tot
            u += 1
        return out

    def genome_load_annotation(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load annotation", self.settings.value("last_genome_dir", ""),
                                              "Annotations (*.gff *.gff3 *.gtf *.bed *.gb *.gbk *.gz);;All files (*)")
        if not path:
            return
        only = None
        try:
            big = os.path.getsize(path) > 50_000_000
        except OSError:
            big = False
        if big and self.genome_contig:
            r = QMessageBox.question(self, "Large annotation", f"This file is large. Load only features on {self.genome_contig}?",
                                     QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
            if r == QMessageBox.Cancel:
                return
            only = self.genome_contig if r == QMessageBox.Yes else None
        prog = QProgressDialog("Reading annotation…", None, 0, 0, self)
        prog.setWindowModality(Qt.WindowModal); prog.show(); QApplication.processEvents()
        try:
            ann = GA.load_annotation(path, only)
        except Exception as e:
            prog.close(); QMessageBox.critical(self, "Annotation", str(e)); return
        prog.close()
        self.annotation = ann
        if not self.genome and self.model.nrows:
            self._enter_reference_mode(ann)
        self.genome_panel.set_annotation(ann)
        n_here = len(ann.genes_by_seq.get(self.genome_contig or "", []))
        if ann.count() and n_here == 0:
            QMessageBox.information(self, "Annotation", f"Loaded {ann.count():,} genes, but none on '{self.genome_contig}'.\n"
                                    f"Sequence names in the file: {', '.join(list(ann.seqids())[:8])}{'…' if len(ann.seqids()) > 8 else ''}")
        self.statusBar().showMessage(f"Annotation: {ann.count():,} genes ({n_here:,} on {self.genome_contig})", 6000)
        self.view.viewport().update()

    def genome_load_synteny(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load PAF", self.settings.value("last_genome_dir", ""), "PAF (*.paf *.paf.gz);;All files (*)")
        if not path:
            return
        try:
            blocks = GA.load_paf(path, min_len=5000, min_mapq=0, query=self.genome_contig)
            if not blocks:
                blocks = GA.load_paf(path, min_len=5000)  # maybe the query names differ; show what we have
        except Exception as e:
            QMessageBox.critical(self, "PAF", str(e)); return
        self.synteny_blocks = blocks
        self.genome_panel.set_synteny(blocks)
        self.statusBar().showMessage(f"Synteny: {len(blocks):,} blocks ≥5 kb", 5000)

    def circular_map(self):
        if not self.model.nrows:
            return
        from .dialogs.circular_dialog import CircularDialog
        ann = self.annotation
        if ann is None:
            QMessageBox.information(self, "Circular map",
                                    "Load an annotation first (Genome ▸ Open GenBank reference…, or Load annotation…).\n"
                                    "The map will still draw the sequence rings without genes.")
        d = CircularDialog(self, seqid=self.genome_contig or self.model.rows[0].name,
                           length=self.proj().ref_len, ann=ann,
                           fetch_seq=lambda s, e: self.ref_ungapped()[s:e],
                           title=(self.genome_contig or self.model.rows[0].name))
        d.positionClicked.connect(lambda p: self._genome_focus(p, p))
        d.geneActivated.connect(self._genome_open_gene)
        if self.genome_panel.isVisible():
            self.genome_panel.region.windowChanged.connect(lambda a, b: d.view.set_window(a, b))
        # keep the focus wedge in sync with the grid
        self._circ = d
        d.show(); self._children.append(d)
        self._grid_scrolled()

    def genome_goto(self):
        if not self.genome_panel.isVisible():
            self.a_g_panel.setChecked(True); self._toggle_genome_panel(True)
        self.genome_panel.region_edit.setFocus(); self.genome_panel.region_edit.selectAll()

    def genome_close(self):
        if self.genome:
            self.genome.close()
        self.genome = None; self.annotation = None; self.synteny_blocks = []
        self.view.feature_provider = None
        self.a_g_panel.setChecked(False); self.genome_panel.hide()
        self.view.viewport().update()

    def _refresh_insertions(self):
        if getattr(self, "_ins_dirty", False) and self.genome_panel.isVisible():
            self.genome_panel.set_insertions(self.insertion_carets())
            self._ins_dirty = False

    def _grid_scrolled(self, *_):
        self._refresh_insertions()
        if not self.genome_panel.isVisible() and getattr(self, "_circ", None) is None:
            return
        if not self.genome_panel.isVisible():
            hs = self.view.horizontalScrollBar()
            c0 = hs.value(); c1 = c0 + self.view._visible_cols()
            proj = self.proj()
            try:
                self._circ.view.set_focus(proj.col_to_ref(c0), max(proj.col_to_ref(c1), proj.col_to_ref(c0) + 1))
            except RuntimeError:
                self._circ = None
            return
        hs = self.view.horizontalScrollBar()
        c0 = hs.value(); c1 = c0 + self.view._visible_cols()
        proj = self.proj()
        s, e = proj.col_to_ref(c0), proj.col_to_ref(c1)
        self.genome_panel.set_focus(s, max(e, s + 1))
        self.genome_panel.ensure_window_contains(s, max(e, s + 1))
        circ = getattr(self, "_circ", None)
        if circ is not None:
            try:
                circ.view.set_focus(s, max(e, s + 1))
            except RuntimeError:
                self._circ = None

    def _goto_columns(self, c0: int, c1: int):
        """Scroll the grid to alignment columns (used by insertion carets)."""
        vis = self.view._visible_cols()
        start = max(0, (c0 + c1) // 2 - vis // 2)
        self.view.horizontalScrollBar().setValue(start)
        if self.model.nrows:
            self.view.select_region(0, self.model.nrows - 1, c0, max(c0, c1 - 1))
        self.statusBar().showMessage(f"Insertion at alignment columns {c0 + 1:,}-{c1:,} "
                                     f"(absent from {self.model.rows[self.ref_index()].name})", 6000)

    def _genome_focus(self, s: int, e: int):
        """Scroll/center the grid on reference span [s,e)."""
        proj = self.proj()
        c0, c1 = proj.span_to_cols(s, max(e, s))
        vis = self.view._visible_cols()
        start = max(0, (c0 + c1) // 2 - vis // 2) if c1 - c0 < vis else c0
        self.view.horizontalScrollBar().setValue(start)
        if self.model.nrows:
            self.view.set_cursor(self.ref_index(), c0)
            self.view.horizontalScrollBar().setValue(start)

    def _translation_regions(self, row: int, c0: int, c1: int):
        """Coding regions for the amino-acid line: annotated CDS on the reference row
        plus any ORF-finder results, each with its own strand and genetic code."""
        out = []
        default_table = int(self.settings.value("default_table", 1))
        if self.annotation and self.genome_contig and row == self.ref_index():
            proj = self.proj()
            u0, u1 = proj.col_to_ref(c0), proj.col_to_ref(c1) + 1
            for g in self.annotation.overlapping(self.genome_contig, u0, u1):
                if g.biotype and g.biotype.lower() not in ("cds", "gene", "mrna", "protein_coding", "bed"):
                    continue
                t = max(g.transcripts, key=lambda t: (len(t.cds), t.end - t.start)) if g.transcripts else None
                segs = (t.cds or t.exons) if t else [(g.start, g.end)]
                table = int(g.attrs.get("transl_table") or (2 if self.model.circular and not g.cytoplasmic else default_table))
                for a, b in segs:
                    ca, cb = proj.span_to_cols(a, b) if not proj.identity else (a, b)
                    if cb > c0 and ca < c1:
                        out.append((ca, cb, g.strand, table, g.name))
        for f in self.model.features:
            if f.row == row and f.type in ("ORF", "CDS") and f.end > c0 and f.start < c1:
                out.append((f.start, f.end, f.strand, int(f.data.get("table", default_table)), f.label))
        return out

    def _genome_features(self, row: int, c0: int, c1: int):
        """Feature provider for the grid: gene models of the genome row as CDS/exon features."""
        if row != self.ref_index() or not self.annotation or not self.genome_contig:
            return []
        proj = self.proj()
        u0, u1 = proj.col_to_ref(c0), proj.col_to_ref(c1) + 1
        def P(a, b):
            return proj.span_to_cols(a, b) if not proj.identity else (a, b)
        out = []
        for g in self.annotation.overlapping(self.genome_contig, u0, u1):
            t = max(g.transcripts, key=lambda t: (len(t.cds), t.end - t.start)) if g.transcripts else None
            if t is None:
                continue
            col = "#d926a9" if getattr(g, "cytoplasmic", False) else ("#1f5fbf" if g.strand > 0 else "#c0392b")
            segs = [(a, b, "CDS") for a, b in t.cds] or [(a, b, "exon") for a, b in t.exons]
            for a, b, kind in segs:
                if b > u0 and a < u1:
                    ca, cb = P(a, b)
                    out.append(Feature(row, ca, cb, g.strand, kind, f"{g.name} {kind}", col))
            if t.cds:
                for a, b in t.utrs():
                    if b > u0 and a < u1:
                        ca, cb = P(a, b)
                        out.append(Feature(row, ca, cb, g.strand, "UTR", f"{g.name} UTR", "#9aa8c7" if g.strand > 0 else "#d9a59c"))
        return out

    def _region_features(self, seqid: str, s: int, e: int) -> list:
        """Gene features for a region re-based to 0 (for 'open region in editor')."""
        feats = []
        if not self.annotation:
            return feats
        for g in self.annotation.overlapping(seqid, s, e):
            col = "#1f5fbf" if g.strand > 0 else "#c0392b"
            for t in g.transcripts:
                for a, b in t.exons:
                    if b > s and a < e:
                        feats.append(Feature(0, max(a, s) - s, min(b, e) - s, g.strand, "exon", f"{g.name} {t.name} exon", col))
                for a, b in t.cds:
                    if b > s and a < e:
                        feats.append(Feature(0, max(a, s) - s, min(b, e) - s, g.strand, "CDS", f"{g.name} {t.name} CDS", col))
                break  # primary transcript only for overlay clarity
        return feats

    def _genome_open_region(self, s: int, e: int):
        if not self.model.nrows:
            return
        seqid = self.genome_contig or self.model.rows[self.ref_index()].name
        c0, c1 = self.proj().span_to_cols(s, e)
        seq = self.ref_seq()[c0:c1]
        m = AlignmentModel([SequenceRow(f"{seqid}:{s + 1}-{e}", seq)], "dna")
        m.features = self._region_features(seqid, s, e)
        m.dirty = False
        w = MainWindow(); w._set_model(m); w.setWindowTitle(f"NeoEdit — {seqid}:{s + 1:,}-{e:,}")
        w.show(); self._children.append(w)

    def _genome_open_gene(self, gene):
        pad = 500
        self._genome_open_region(max(0, gene.start - pad), min(self.proj().ref_len, gene.end + pad))

    # ------------------------------------------------------------ help
    def shortcuts(self):
        txt = """Navigation
  Arrows / PgUp / PgDn / Home / End   move cursor (Shift extends selection, Ctrl+Left/Right jumps 10)
  Ctrl+G                              go to column
  Ctrl+wheel                          zoom;  Shift+wheel  horizontal scroll

Gap editing (any mode)
  Space or -           insert gap at cursor (in current row; all selected rows if a block is selected)
  Ctrl+Space           insert gap column in all sequences
  Delete / Backspace   delete gap at / before cursor (only gaps in select mode)
  Ctrl+Delete          delete gap column (only if every sequence has a gap there)
  Drag a selected block left/right to slide residues over gaps (BioEdit-style)
  Right-click            insert/delete a gap in the clicked (or selected) sequence, or in all
                         unselected sequences - choose under Edit > Right-click action
  Shift+right-click      context menu

Modes (BioEdit-style; F5 / F6 / F7)
  Select / Slide   box-select; drag the selection to slide it over gaps (crunch ahead, open behind);
                   hold Shift to move the entire sequence downstream instead. The toolbar
                   "downstream" toggle swaps which of these is the default.
  Edit             place the cursor and type; Insert / Overwrite choice appears on the toolbar
                   (Insert key toggles); Delete/Backspace remove residues.
  Grab & Drag      press on a single residue and drag it; Shift drags everything downstream of it.

Selection
  Click & drag in the grid  rectangular block;  click a name  whole row (Ctrl/Shift to extend)
  Click the ruler           select a column;     Esc clears the selection

Other
  Ctrl+Shift+I import from NCBI / Ensembl / UCSC, Ctrl+C copy FASTA, Ctrl+V paste sequences, Ctrl+F find, F3 find next,
  Ctrl+Shift+R reverse complement (as in BioEdit; Ctrl+R also works),
  Ctrl+T translation overlay, Ctrl+Shift+T translate, Ctrl+Shift+O ORF finder, Ctrl+Shift+P primer design,
  Ctrl+Shift+X restriction sites, Ctrl+M align with MAFFT
"""
        TextDialog(self, "Keyboard shortcuts", txt).exec()

    def about(self):
        QMessageBox.about(self, "About NeoEdit",
                          f"<b>NeoEdit</b> {__version__}<br>A modern, open, cross-platform sequence alignment editor "
                          f"inspired by Tom Hall's BioEdit.<br><br>Python {sys.version.split()[0]}, PySide6, Biopython, primer3-py.")
