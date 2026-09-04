# Contributing to NeoEdit

Thanks for your interest in NeoEdit. It's a small, actively-developed project — bug reports,
doc fixes, and pull requests of any size are welcome and get read promptly.

## Reporting bugs or requesting features

- Search [existing issues](https://github.com/evozoa/NeoEdit/issues) first.
- For a bug: your OS, the NeoEdit version (`neoedit --self-test`, or Help ▸ About), steps to
  reproduce, and — if you can — a small FASTA/GenBank file that reproduces it.
- For a feature request: if it's meant to match a specific BioEdit behavior, say which one
  (a menu path, a keyboard shortcut, a file-format detail); if it's new territory, describe the
  use case that's currently blocked.

## Before submitting a pull request

- For anything beyond a small, obvious fix, open an issue first to agree on the approach —
  this saves both of us from wasted work.
- Keep a PR focused on one change. Large, mixed-purpose PRs are slow to review and slow to merge.

## Development setup

```bash
git clone https://github.com/evozoa/NeoEdit
cd NeoEdit
pip install -e ".[dev]"
pytest
```

MAFFT is only needed to actually run alignments (`conda install -c bioconda mafft`,
`brew install mafft`, or the Windows installer) — not to run the test suite.

## Code layout

- `src/neoedit/model/` and `src/neoedit/analysis/` are Qt-free, plain-Python, and unit-testable
  in isolation — put logic here whenever it doesn't need a widget.
- `src/neoedit/ui/` is PySide6: the custom-painted alignment grid (`alignment_view.py`), the
  main window, dialogs.
- `src/neoedit/remote/` talks to NCBI/Ensembl/UCSC.

See the README's *Layout* section for the full picture.

## Tests

- New behavior needs a test. Model/analysis logic gets a plain pytest test; UI behavior gets an
  offscreen smoke test (`QT_QPA_PLATFORM=offscreen`) — see `tests/test_gui_smoke.py` for the
  pattern (synthetic Qt events against a real `AlignmentView`, no display needed).
- Run the full suite (`pytest`) before opening a PR. CI runs it again on every push and PR.

## Matching BioEdit behavior

Much of NeoEdit exists to reproduce a specific BioEdit convention — keyboard gap editing, a
color table, the `.bio` file format, a particular menu's wording. If a change touches something
with a BioEdit precedent, say so in the PR description; if you're deliberately diverging from
BioEdit, say why.

## Style

No enforced linter/formatter yet — match the surrounding code. (This may change; check for a
`ruff`/`black` config before assuming there still isn't one.)

## Releases

The version string is kept in sync between `pyproject.toml` and `src/neoedit/__init__.py`
(CI checks this on tag pushes). Tagging `vX.Y.Z` on `main` builds the Windows/Mac/Linux
installers and publishes a GitHub Release.

## A note on how this project is built

NeoEdit's development has involved AI pair-programming (Claude Code) alongside human design
decisions and review — every change is understood, tested, and stood behind by a human before
it ships. You're welcome to use whatever tools you like too; you're responsible for
understanding and standing behind whatever you submit.

## Support and maintenance

NeoEdit is currently maintained by MWS ([@evozoa](https://github.com/evozoa)). Issues and PRs
are reviewed on a best-effort basis — there's no SLA, but feel free to ping a thread that's
gone quiet for a couple of weeks.

Please be respectful and constructive. This is a small project, and a good experience for
contributors is worth protecting.
