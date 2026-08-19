# NeoEdit

<img src="src/neoedit/resources/icons/neoedit.png" width="96" align="right" alt="NeoEdit icon">

A modern, open-source, cross-platform (Windows / macOS / Linux) sequence alignment
editor inspired by Tom Hall's **BioEdit** — NeoEdit is a from-scratch rebuild in Python with Qt (PySide6),
Biopython and Primer3.

BioEdit is closed-source and frozen; this project recreates the parts people
actually use, on a codebase that can keep evolving:

* **Alignment editor** – BioEdit-style grid, keyboard gap editing (Space/Delete/
  Backspace), insert/overwrite typing modes, block selection, drag-to-slide
  residues over gaps, column ops, undo/redo, zoom, find (gap-tolerant/regex).
* **Coloring** – BioEdit's own color tables are the default (parsed from the
  `color.tab` files shipped with BioEdit: A green / C blue / G black / T red and the
  matching amino-acid table), plus BioEdit BLOSUM / PAM250 / Kyte-Doolittle /
  Manuel Ruiz / codon-degeneracy tables and Clustal X / Zappo / Taylor /
  hydrophobicity. BioEdit "normal" view (colored letters) by default; Ctrl+I
  toggles "inverse" view (colored backgrounds). Shade identities (black, as in
  BioEdit) vs consensus or first sequence, identities-as-dots, consensus row with
  adjustable threshold.
* **Sequence tools** – reverse-complement, complement, reverse, case, DNA↔RNA,
  translation (any NCBI genetic code, any frame, keep-aligned or six-frame report),
  live translation overlay under DNA.
* **ORF finder (MitoFinder-style)** – alternate genetic codes (default vertebrate
  mito), alternative start codons, both strands, 5′/3′ partial ORFs (incomplete
  stop codons), nested ORFs; results as features; export GFF3 / GenBank feature
  table / protein & nucleotide FASTA.
* **Primer design** – Primer3 (bundled via `primer3-py`): target region from the
  selection, size/Tm/GC constraints, salt/Mg/dNTP, hairpin/dimer metrics, primers
  drawn as features, mismatch check of each primer against every sequence in the
  alignment (for conserved-site design), CSV/FASTA export.
* **Genome context** – three synced tiers: contig/chromosome overview, gene-model
  and synteny region view (GFF3/GTF/BED/GenBank + PAF), and the alignment grid;
  indexed FASTA so whole chromosomes open instantly. GenBank references (mitogenome/
  plasmid) populate the gene view from their own features; MAFFT `--add --keeplength`
  anchors a population sample to the reference without moving its coordinates, with a
  per-site variation track.
* **Circular map** – mitogenome/plasmid view with strand-separated gene arrows, GC
  content and GC skew rings, focal wedge synced to the grid; SVG/PNG export.
* **Primer design across an alignment** – conserved (universal) or discriminating
  (eDNA) primers using sequence groups, conservation masking, 3'-weighted mismatch
  scoring, IUPAC-degenerate options and an in-silico PCR table.
* **Analysis** – identity matrix, entropy/identity conservation plots, sequence
  statistics, consensus report, MAFFT alignment (local executable; strategy/threads/--adjustdirection options),
  NCBI BLAST launcher.
* **Formats** – FASTA, Clustal, PHYLIP, NEXUS, Stockholm, GenBank, EMBL, MSF and a
  best-effort reader for legacy BioEdit `.bio` files.

## Install & run

```bash
pip install -e .            # or: pip install -e ".[dev]" for tests
neoedit examples/cox1_demo.fasta
# or
python -m neoedit
```

Requires Python ≥ 3.10. MAFFT is optional and found on PATH or set in *Edit ▸ Preferences*
(`conda install -c bioconda mafft`, `brew install mafft`, or the MAFFT Windows installer).

## Tests

```bash
pytest            # model, I/O, analysis + offscreen GUI smoke test
```

## Layout

```
src/neoedit/
  model/      alignment.py (pure-Python model + undo), io.py, colors.py
  analysis/   translate.py, orf_finder.py, primers.py, external.py   (no Qt)
  ui/         alignment_view.py (custom-painted grid), main_window.py,
              feature_track.py, dialogs/
tests/        pytest suite
examples/     demo alignment
```

Everything under `model/` and `analysis/` is Qt-free and usable as a library.

## Roadmap

Plasmid maps, restriction-site search, BLAST-based
gene naming for ORFs, PyInstaller bundles.
