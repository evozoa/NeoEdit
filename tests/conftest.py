import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["NEOEDIT_SETTINGS"] = "neoedit-tests"   # keep tests out of the user's settings

import pytest


@pytest.fixture(autouse=True)
def _destroy_leftover_windows():
    """Destroy the windows a test leaves open while Qt is still running. Tests make dialogs
    with no parent; left to Python, they can be freed at interpreter exit after Qt's own
    teardown, which crashes PySide6 >= 6.12 ("shared QObject was deleted directly",
    segfaults, bus errors) after every test has passed."""
    yield
    try:
        from PySide6.QtCore import QCoreApplication, QEvent
        from PySide6.QtWidgets import QApplication
    except ImportError:
        return
    app = QApplication.instance()
    if app is None:
        return
    for w in app.topLevelWidgets():
        w.hide()             # not close(): a main window with edits would ask to save and block
        w.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.DeferredDelete)
