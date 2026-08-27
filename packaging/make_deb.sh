#!/usr/bin/env bash
# Wrap the PyInstaller one-dir build in a Debian package (Ubuntu 22.04+, Debian, Mint).
#   packaging/make_deb.sh dist/NeoEdit dist/NeoEdit-Linux-x86_64.deb [VERSION]
# Installs to /opt/neoedit with a launcher-menu entry, hicolor icons and a `neoedit` command.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
APP="${1:-dist/NeoEdit}"
OUT="${2:-dist/NeoEdit-Linux-x86_64.deb}"
VERSION="${3:-$(python3 -c "import re;print(re.search(r'__version__ = \"([^\"]+)\"', open('$ROOT/src/neoedit/__init__.py').read()).group(1))")}"
[ -x "$APP/NeoEdit" ] || { echo "no PyInstaller build at $APP (expected $APP/NeoEdit)" >&2; exit 1; }
case "$(uname -m)" in x86_64) ARCH=amd64 ;; aarch64) ARCH=arm64 ;; *) ARCH="$(dpkg --print-architecture)" ;; esac

STAGE="$(mktemp -d)"
chmod 755 "$STAGE"                                   # mktemp gives 0700
trap 'rm -rf "$STAGE"' EXIT
ICONS="$ROOT/src/neoedit/resources/icons"

# /opt/neoedit: the whole one-dir bundle (NeoEdit + _internal/, incl. bundled MAFFT if any)
mkdir -p "$STAGE/opt"
cp -R "$APP" "$STAGE/opt/neoedit"
chmod -R u+rwX,go+rX,go-w "$STAGE/opt/neoedit"

# `neoedit` on the PATH, launcher-menu entry, icons
install -d "$STAGE/usr/bin"
ln -s /opt/neoedit/NeoEdit "$STAGE/usr/bin/neoedit"
install -Dm644 "$HERE/neoedit.desktop" "$STAGE/usr/share/applications/neoedit.desktop"
for S in 16 24 32 48 64 128 256; do
  install -Dm644 "$ICONS/neoedit_${S}.png" "$STAGE/usr/share/icons/hicolor/${S}x${S}/apps/neoedit.png"
done

# Copyright: MIT for NeoEdit, GPL notice for the bundled MAFFT
DOC="$STAGE/usr/share/doc/neoedit"
install -d "$DOC"
{
  echo "NeoEdit $VERSION - https://github.com/evozoa/NeoEdit"; echo
  cat "$ROOT/LICENSE"
  MAFFT="$APP/_internal/mafft"
  if [ -d "$MAFFT" ]; then
    echo; echo "----------------------------------------------------------------------"; echo
    cat "$HERE/MAFFT-NOTICE.txt"
    [ -f "$MAFFT/copyright" ] && { echo; cat "$MAFFT/copyright"; }
  fi
} > "$DOC/copyright"

# Control file. Depends = what the Qt xcb platform plugin loads from the system; Ubuntu desktops
# lack some of these (libxcb-cursor0), so apt pulls them in when the .deb is installed.
# QT_APT_LIBS is the same list the workflow installs on the runners (build.yml `env:`).
QT_LIBS="${QT_APT_LIBS:-libegl1 libgl1 libxkbcommon0 libxkbcommon-x11-0 libdbus-1-3 libfontconfig1 libglib2.0-0 libxcb-cursor0 libxcb-icccm4 libxcb-keysyms1 libxcb-shape0 libxcb-xkb1 libxcb-image0 libxcb-render-util0 libxcb-util1}"
DEPENDS="libc6 (>= 2.35), $(echo $QT_LIBS | sed 's/ /, /g')"
SIZE_KB="$(du -sk "$STAGE" | cut -f1)"
mkdir -p "$STAGE/DEBIAN"
cat > "$STAGE/DEBIAN/control" <<CONTROL
Package: neoedit
Version: $VERSION
Section: science
Priority: optional
Architecture: $ARCH
Maintainer: evozoa <evozoa@gmail.com>
Installed-Size: $SIZE_KB
Depends: $DEPENDS
Homepage: https://github.com/evozoa/NeoEdit
Description: Sequence alignment editor and genome viewer (BioEdit-style)
 NeoEdit is a modern, open-source sequence alignment editor inspired by BioEdit:
 FASTA/GenBank/EMBL import, translation, ORF and primer tools, restriction maps,
 NCBI/Ensembl/UCSC import and a genome viewer. MAFFT is bundled for alignment.
CONTROL

mkdir -p "$(dirname "$OUT")"
rm -f "$OUT"
dpkg-deb --build --root-owner-group "$STAGE" "$OUT"
dpkg-deb --info "$OUT"
dpkg-deb --contents "$OUT" | grep -v "/_internal/" | awk '{print $1, $6, $7, $8}'
if command -v desktop-file-validate >/dev/null; then desktop-file-validate "$HERE/neoedit.desktop"; fi
ls -lh "$OUT"
