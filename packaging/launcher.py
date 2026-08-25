"""PyInstaller entry point. (Not neoedit/__main__.py: that runs main() at import time.)"""
import sys

from neoedit.app import main

if __name__ == "__main__":
    sys.exit(main())
