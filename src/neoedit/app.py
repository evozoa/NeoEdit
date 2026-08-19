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
    from .ui.main_window import MainWindow
    files = [a for a in argv[1:] if not a.startswith("-")]
    w = MainWindow(files)
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
