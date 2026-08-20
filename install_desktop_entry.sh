#!/usr/bin/env bash
# Register NeoEdit with the desktop environment.
#
# On WSL this is what puts the NeoEdit icon on the Windows taskbar and in the
# Start Menu: WSLg only indexes /usr/share/applications, so a user-level entry in
# ~/.local/share/applications is ignored and Windows falls back to the generic
# WSL penguin. Needs sudo for the system-wide copy.
set -e
SRC="$(cd "$(dirname "$0")" && pwd)"
ICONS="$SRC/src/neoedit/resources/icons"

sudo install -Dm644 /dev/stdin /usr/share/applications/neoedit.desktop <<DESKTOP
[Desktop Entry]
Type=Application
Name=NeoEdit
GenericName=Sequence alignment editor
Comment=Sequence alignment editor and genome viewer
Exec=$SRC/launch.sh %F
Icon=neoedit
Terminal=false
StartupNotify=true
StartupWMClass=neoedit
Categories=Science;Biology;Education;
MimeType=text/x-fasta;chemical/seq-na-fasta;chemical/seq-aa-fasta;
Keywords=DNA;protein;alignment;FASTA;GenBank;BioEdit;
DESKTOP

for S in 16 24 32 48 64 128 256; do
  sudo install -Dm644 "$ICONS/neoedit_${S}.png" \
    "/usr/share/icons/hicolor/${S}x${S}/apps/neoedit.png"
done
sudo gtk-update-icon-cache -f -t /usr/share/icons/hicolor 2>/dev/null || true
sudo update-desktop-database /usr/share/applications 2>/dev/null || true

echo "Installed. On WSL, the Start Menu entry ('NeoEdit (Ubuntu)') appears within a"
echo "minute; launch from there (or pin it) to get the NeoEdit icon on the taskbar."
