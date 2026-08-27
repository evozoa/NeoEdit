# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for the standalone NeoEdit build (Windows / macOS / Linux).
#   pyinstaller --noconfirm --clean packaging/neoedit.spec
# Produces dist/NeoEdit/ (one-dir) and, on macOS, dist/NeoEdit.app. Linux: packaging/make_deb.sh
# wraps dist/NeoEdit/ in a .deb (icons/menu entry come from packaging/neoedit.desktop).
# If packaging/mafft/ exists (see .github/workflows/build.yml) it is shipped inside the app
# and found by neoedit.analysis.external.bundled_mafft().
import os
import sys

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

HERE = os.path.dirname(os.path.abspath(SPEC))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)
import neoedit  # noqa: E402

VERSION = neoedit.__version__
ICON_DIR = os.path.join(SRC, "neoedit", "resources", "icons")
ICON = os.path.join(ICON_DIR, "neoedit.icns" if sys.platform == "darwin" else "neoedit.ico")

datas = [
    (os.path.join(SRC, "neoedit", "resources"), os.path.join("neoedit", "resources")),
    (os.path.join(ROOT, "LICENSE"), "."),
]
mafft_dir = os.path.join(HERE, "mafft")
if os.path.isdir(mafft_dir):
    datas.append((mafft_dir, "mafft"))
# Bio.Align.substitution_matrices.load("BLOSUM62") reads a data file by name; without this
# it silently degrades to match/mismatch scoring.
datas += collect_data_files("Bio")
# primer3-py loads its thermodynamic parameter files (primer3_config/) from the package dir;
# without them the C extension crashes on first use.
datas += collect_data_files("primer3")

hiddenimports = (
    collect_submodules("Bio.SeqIO")
    + collect_submodules("Bio.AlignIO")
    + collect_submodules("primer3")          # primer3.bindings is imported lazily
    + ["PySide6.QtSvg", "certifi"]
)

excludes = [
    "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtWebEngineQuick",
    "PySide6.QtWebChannel", "PySide6.QtWebSockets", "PySide6.QtNetwork", "PySide6.QtPrintSupport",
    "PySide6.QtCharts", "PySide6.QtDataVisualization", "PySide6.QtQml", "PySide6.QtQuick",
    "PySide6.QtQuick3D", "PySide6.QtQuickWidgets", "PySide6.Qt3DCore", "PySide6.Qt3DRender",
    "PySide6.Qt3DInput", "PySide6.Qt3DLogic", "PySide6.Qt3DAnimation", "PySide6.Qt3DExtras",
    "PySide6.QtMultimedia", "PySide6.QtMultimediaWidgets", "PySide6.QtPositioning",
    "PySide6.QtLocation", "PySide6.QtBluetooth", "PySide6.QtNfc", "PySide6.QtSensors",
    "PySide6.QtSerialPort", "PySide6.QtSerialBus", "PySide6.QtRemoteObjects", "PySide6.QtScxml",
    "PySide6.QtStateMachine", "PySide6.QtTest", "PySide6.QtSql", "PySide6.QtDBus",
    "PySide6.QtDesigner", "PySide6.QtHelp", "PySide6.QtUiTools", "PySide6.QtPdf",
    "PySide6.QtPdfWidgets", "PySide6.QtTextToSpeech", "PySide6.QtHttpServer", "PySide6.QtGraphs",
    "PySide6.QtSpatialAudio", "PySide6.QtOpenGL", "PySide6.QtOpenGLWidgets", "PySide6.QtConcurrent",
    "PySide6.QtXml", "PySide6.QtSvgWidgets", "PySide6.QtNetworkAuth", "PySide6.QtAsyncio",
    "pytest", "tkinter", "_tkinter", "IPython", "matplotlib", "scipy", "pandas",
]

a = Analysis(
    [os.path.join(HERE, "launcher.py")],
    pathex=[SRC],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    exclude_binaries=True,
    name="NeoEdit",
    debug=False,
    strip=False,
    upx=False,
    console=False,
    icon=ICON,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="NeoEdit",
)

if sys.platform == "darwin":
    EXTS = ["fasta", "fas", "fa", "fna", "faa", "aln", "bio", "gb", "gbk", "genbank", "embl",
            "phy", "phylip", "nex", "nexus", "sto", "msf", "gff", "gff3", "gtf", "bed"]
    app = BUNDLE(
        coll,
        name="NeoEdit.app",
        icon=ICON,
        bundle_identifier="org.neoedit.NeoEdit",
        info_plist={
            "CFBundleName": "NeoEdit",
            "CFBundleDisplayName": "NeoEdit",
            "CFBundleShortVersionString": VERSION,
            "CFBundleVersion": VERSION,
            "NSHighResolutionCapable": True,
            "NSHumanReadableCopyright": "MIT License. MAFFT (bundled) is GPL, (c) Kazutaka Katoh.",
            "LSMinimumSystemVersion": "12.0",
            "LSApplicationCategoryType": "public.app-category.education",
            "CFBundleDocumentTypes": [{
                "CFBundleTypeName": "Sequence file",
                "CFBundleTypeRole": "Editor",
                "LSHandlerRank": "Alternate",
                "CFBundleTypeExtensions": EXTS,
                "CFBundleTypeIconFile": os.path.basename(ICON),
            }],
        },
    )
