"""Application entry point: `python -m neoedit [files...]`, `neoedit`, or the frozen NeoEdit binary."""
from __future__ import annotations

import os
import sys


def _install_excepthook():
    """Windowed (frozen) builds have no console: show unhandled exceptions in a dialog."""
    import traceback

    def hook(exc_type, exc, tb):
        text = "".join(traceback.format_exception(exc_type, exc, tb))
        try:
            sys.__stderr__ and sys.__stderr__.write(text)
        except Exception:
            pass
        try:
            from PySide6.QtWidgets import QApplication, QMessageBox
            if QApplication.instance() is not None:
                QMessageBox.critical(None, "NeoEdit - unexpected error",
                                     f"{exc_type.__name__}: {exc}\n\n{text[-3000:]}")
        except Exception:
            pass

    sys.excepthook = hook


def _make_app(argv):
    from PySide6.QtCore import QEvent
    from PySide6.QtWidgets import QApplication

    class NeoEditApp(QApplication):
        """QApplication that also accepts files opened from the OS (macOS Finder double-click,
        "Open With", Dock drops), which arrive as QFileOpenEvent rather than argv."""

        def __init__(self, argv):
            super().__init__(argv)
            self.window = None
            self.pending_files: list[str] = []

        def event(self, ev):
            if ev.type() == QEvent.Type.FileOpen:
                path = ev.file()
                if path:
                    if self.window is not None:
                        self.window.open_path(path)
                    else:
                        self.pending_files.append(path)
                return True
            return super().event(ev)

    app = NeoEditApp(argv)
    app.setApplicationName("NeoEdit")
    app.setOrganizationName("neoedit")
    app.setDesktopFileName("neoedit")
    from .ui.icons import app_icon
    app.setWindowIcon(app_icon())
    return app


def _parse(argv):
    import argparse
    ap = argparse.ArgumentParser(prog="neoedit", description="NeoEdit sequence alignment editor")
    ap.add_argument("files", nargs="*", help="alignment/sequence files to open")
    ap.add_argument("--genome", help="open a (large) genome FASTA in genome mode")
    ap.add_argument("--contig", help="contig/chromosome to load with --genome")
    ap.add_argument("--gff", help="annotation (GFF3/GTF/BED) to load with --genome")
    ap.add_argument("--paf", help="synteny PAF to load with --genome")
    ap.add_argument("--self-test", metavar="REPORT.json",
                    help="check that this installation works (resources, MAFFT, network), "
                         "write a JSON report and exit 0/1")
    # Finder passes -psn_0_NNN to apps it launches; Qt may be given its own flags (-style ...).
    # Never let those kill a windowed build with an argparse SystemExit.
    args = [a for a in argv[1:] if not a.startswith("-psn_")]
    ns, unknown = ap.parse_known_args(args)
    return ns


def main(argv=None):
    import multiprocessing
    multiprocessing.freeze_support()
    argv = list(sys.argv if argv is None else argv)
    ns = _parse(argv)
    if ns.self_test:
        return self_test(ns.self_test)
    _install_excepthook()
    app = _make_app(argv)
    from .ui.main_window import MainWindow
    w = MainWindow(ns.files)
    app.window = w
    for p in app.pending_files:
        w.open_path(p)
    app.pending_files.clear()
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


def self_test(report_path: str) -> int:
    """Exercise the parts of a packaged build that fail *silently* when mis-bundled.

    Runs headless (offscreen). Writes {check: {"ok": bool, "detail": str}} to `report_path`
    and returns 0 when every check passed. Set NEOEDIT_SELFTEST_OFFLINE=1 to skip the
    network check."""
    import json
    import traceback

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    results: dict[str, dict] = {}
    app = _make_app([sys.argv[0]])       # QPixmap/QIcon need a QGuiApplication first

    def check(name, fn):
        try:
            print(f"... {name}", flush=True)
        except Exception:
            pass
        try:
            detail = fn()
            results[name] = {"ok": True, "detail": str(detail)}
        except Exception as e:                       # noqa: BLE001 - report everything
            results[name] = {"ok": False, "detail": f"{e.__class__.__name__}: {e}\n"
                                                    f"{traceback.format_exc()[-1500:]}"}

    def c_env():
        return (f"python {sys.version.split()[0]}, frozen={getattr(sys, 'frozen', False)}, "
                f"exe={sys.executable}, _MEIPASS={getattr(sys, '_MEIPASS', None)}")

    def c_colors():
        from .model import colors
        if not colors._BIOEDIT["dna"] or not colors._BIOEDIT["protein"]:
            raise RuntimeError("BioEdit colour tables not found (fell back to built-in colours): "
                               f"looked in {colors._TABLE_DIR}")
        return f"{len(colors._BIOEDIT['dna'])} DNA + {len(colors._BIOEDIT['protein'])} protein colours"

    def c_icon():
        ic = app.windowIcon()
        if ic.isNull():
            raise RuntimeError("application icon not found")
        return f"{len(ic.availableSizes())} sizes"

    def c_blosum():
        from Bio.Align import substitution_matrices
        m = substitution_matrices.load("BLOSUM62")
        return f"BLOSUM62 {m.shape}"

    def c_peptides():
        from .analysis.genbank_export import bundled_mdp_peptides
        p = bundled_mdp_peptides()
        if not p or not os.path.exists(p):
            raise RuntimeError(f"mdp_reference_peptides.fasta missing ({p})")
        return p

    def c_formats():
        from .model import io as mio
        from .model.alignment import SequenceRow
        text = ">a\nATGCATGCAT\n>b\nATGCATGCAA\n"
        m = mio.loads(text, "fasta")
        for fmt in ("fasta", "clustal", "genbank", "nexus", "phylip", "stockholm"):
            mio.dumps(m, fmt)
        return "fasta/clustal/genbank/nexus/phylip/stockholm round-trips"

    def c_primer3():
        import primer3
        from primer3 import bindings
        return f"primer3-py {primer3.__version__}: Tm(ACGTACGTACGTACGTACGT)={bindings.calc_tm('ACGTACGTACGTACGTACGT'):.1f}"

    def c_restriction():
        from Bio.Restriction import AllEnzymes, EcoRI
        return f"{len(AllEnzymes)} enzymes; EcoRI site {EcoRI.site}"

    def c_mafft():
        from .analysis import external
        from .model.alignment import SequenceRow
        exe = external.find_mafft()
        if not exe:
            raise RuntimeError("MAFFT not found (bundled dir: %s)" % external.bundled_mafft_dir())
        ver = external.mafft_version(exe)
        if ver == "unknown":
            raise RuntimeError(f"MAFFT at {exe} does not run")
        rows = [SequenceRow("s1", "ATGCATGCATGCAAGGTTCCAAGGTTAACCGG"),
                SequenceRow("s2", "ATGCATGCAAGGTTCCAAGGTTAACCGG"),
                SequenceRow("s3", "ATGCATGCATGCAAGGTTCCAAGGTTCCAACCGG")]
        out = external.run_mafft(rows, exe, strategy=1, timeout=300)
        lens = {len(r.seq) for r in out}
        if len(lens) != 1:
            raise RuntimeError(f"MAFFT output not aligned: {lens}")
        return f"{exe} ({ver}); aligned 3 seqs to {lens.pop()} columns"

    def c_window():
        import tempfile
        from .ui.main_window import MainWindow
        w = MainWindow()
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "selftest.fasta")
            with open(path, "w") as fh:
                fh.write(">seq1 self-test\nATGCATGCATGCAAGGTTCC\n>seq2\nATGCATGC--GCAAGGTTCC\n")
            w.open_path(path)
        n = w.model.nrows
        if n != 2:
            raise RuntimeError(f"opened file has {n} sequences, expected 2")
        w.close()
        return f"MainWindow opened {n} sequences"

    def c_https():
        if os.environ.get("NEOEDIT_SELFTEST_OFFLINE"):
            return "skipped (NEOEDIT_SELFTEST_OFFLINE)"
        from .remote.http import http_get, ca_bundle
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?retmode=json"
        body = http_get(url, timeout=30, retries=2)
        if b"einforesult" not in body:
            raise RuntimeError(f"unexpected reply from NCBI: {body[:200]!r}")
        return f"NCBI eutils reachable over TLS (CA bundle: {ca_bundle() or 'system'})"

    for name, fn in [("environment", c_env), ("colour_tables", c_colors), ("icon", c_icon),
                     ("blosum62", c_blosum), ("mdp_peptides", c_peptides), ("formats", c_formats),
                     ("primer3", c_primer3), ("restriction", c_restriction), ("mafft", c_mafft),
                     ("main_window", c_window), ("https", c_https)]:
        check(name, fn)

    ok = all(r["ok"] for r in results.values())
    results["_summary"] = {"ok": ok, "detail": f"{sum(r['ok'] for r in results.values())}/{len(results)} passed"}
    with open(report_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    try:
        for name, r in results.items():
            print(("PASS " if r["ok"] else "FAIL ") + name + ": " + r["detail"].splitlines()[0])
    except Exception:
        pass
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
