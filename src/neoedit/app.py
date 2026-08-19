"""Application entry point: `python -m neoedit [files...]` or `neoedit`."""
from __future__ import annotations

import sys


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import Qt
    app = QApplication(argv)
    app.setApplicationName("NeoEdit")
    app.setOrganizationName("neoedit")
    app.setDesktopFileName("neoedit")
    from .ui.icons import app_icon
    app.setWindowIcon(app_icon())
    from .ui.main_window import MainWindow
    import argparse
    ap = argparse.ArgumentParser(prog="neoedit", description="NeoEdit sequence alignment editor")
    ap.add_argument("files", nargs="*", help="alignment/sequence files to open")
    ap.add_argument("--genome", help="open a (large) genome FASTA in genome mode")
    ap.add_argument("--contig", help="contig/chromosome to load with --genome")
    ap.add_argument("--gff", help="annotation (GFF3/GTF/BED) to load with --genome")
    ap.add_argument("--paf", help="synteny PAF to load with --genome")
    ns = ap.parse_args(argv[1:])
    w = MainWindow(ns.files)
    w.show()
    if ns.genome:
        w.open_genome_path(ns.genome, ns.contig)
        if ns.gff:
            w.load_annotation_path(ns.gff, ns.contig)
        if ns.paf:
            from .genome import annotations as GA
            w.synteny_blocks = GA.load_paf(ns.paf, min_len=5000, query=ns.contig)
            w.genome_panel.set_synteny(w.synteny_blocks)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
