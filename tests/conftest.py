import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ["NEOEDIT_SETTINGS"] = "neoedit-tests"   # keep tests out of the user's settings
