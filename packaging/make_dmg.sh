#!/usr/bin/env bash
# Wrap dist/NeoEdit.app in a compressed disk image with an Applications shortcut.
#   packaging/make_dmg.sh dist/NeoEdit.app dist/NeoEdit-macOS-AppleSilicon.dmg
set -euo pipefail
APP="${1:-dist/NeoEdit.app}"
OUT="${2:-dist/NeoEdit.dmg}"
[ -d "$APP" ] || { echo "no app bundle at $APP" >&2; exit 1; }
STAGE="$(mktemp -d)"
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
rm -f "$OUT"
for attempt in 1 2 3; do            # hdiutil occasionally fails with "Resource busy" on CI
  hdiutil create -volname "NeoEdit" -srcfolder "$STAGE" -ov -format UDZO "$OUT" && break
  [ "$attempt" = 3 ] && exit 1
  sleep 5
done
rm -rf "$STAGE"
ls -lh "$OUT"
