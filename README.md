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
  translation (any NCBI genetic code, any frame, keep-aligned or six-frame report).
* **Amino-acid line** – an interleaved protein row under each nucleotide sequence.
  Annotated ORFs and CDS features are translated in *their own* frame with *their own*
  genetic code (mitochondrial genes with table 2, cytoplasmically-translated MDPs with
  table 1) and shaded a neutral grey, so coding regions are shown on the protein line
  rather than washing out the nucleotides; outside them you choose the frame (+1…−3)
  and code, or leave it blank.
* **ORF finder (MitoFinder-style)** – alternate genetic codes (default vertebrate
  mito), alternative start codons, both strands, 5′/3′ partial ORFs (incomplete
  stop codons), nested ORFs; results as features; export GFF3 / GenBank feature
  table / protein & nucleotide FASTA.
* **Primer design** – Primer3 (bundled via `primer3-py`): target region from the
  selection, size/Tm/GC constraints, salt/Mg/dNTP, hairpin/dimer metrics, primers
  drawn as features, mismatch check of each primer against every sequence in the
  alignment (for conserved-site design), CSV/FASTA export.
* **Pinned reference** – any sequence can be pinned (Ctrl+Shift+E) to anchor the
  chromosome map, gene models and coordinates, so all three tiers always agree; the
  pin is stored in `.bio` files as BioEdit's numbering-mask slot. Insertions the other
  sequences carry but the reference lacks appear as clickable carets that jump the
  grid to the inserted bases.
* **Genome context** – three synced tiers: contig/chromosome overview, gene-model
  and synteny region view (GFF3/GTF/BED/GenBank + PAF), and the alignment grid;
  indexed FASTA so whole chromosomes open instantly. GenBank references (mitogenome/
  plasmid) populate the gene view from their own features; MAFFT `--add --keeplength`
  anchors a population sample to the reference without moving its coordinates, with a
  per-site variation track.
* **Mitochondrial-derived peptides** – `examples/mito/NC_012920_MDP.gb` is the human
  rCRS with humanin, MOTS-c and SHLP1-6 annotated as CDS features carrying
  `/transl_table=1`; pin it as the reference and the MDPs appear by name, in their own
  color, on the map, the circle and the grid. Peptides are also shipped as
  `resources/mdp_reference_peptides.fasta` for naming ORFs in other taxa by similarity.
* **Circular map** – mitogenome/plasmid view with strand-separated gene arrows, GC
  content and GC skew rings, focal wedge synced to the grid; SVG/PNG export.
* **Primer design across an alignment** – conserved (universal) or discriminating
  (eDNA) primers using sequence groups, conservation masking, 3'-weighted mismatch
  scoring, IUPAC-degenerate options and an in-silico PCR table.
* **Import from NCBI and Ensembl** (`Ctrl+Shift+I`) – fetch GenBank/FASTA records from
  NCBI Entrez by accession/GI (with sub-range and strand) or by searching Entrez and
  ticking hits; fetch genes (symbol or stable ID), regions, cDNA/CDS/protein from the
  Ensembl REST API **pinned to release 116** (archive server first, live server as
  fallback, with the actual release shown). Ensembl genomic imports are written as GenBank
  with gene / mRNA / tRNA / rRNA / CDS features built from Ensembl's gene models, so the
  gene view and the amino-acid line work on them; minus-strand genes can be oriented to
  their own strand. Downloads are kept as files (default `~/Downloads/NeoEdit`) and appear
  in *Open recent*.
* **Restriction sites** – REBASE enzymes (Biopython) with supplier/site-length/ends
  filters; map one sequence, or compare across an alignment to find enzymes that
  *distinguish* sequences; unique cutters, non-cutters, circular molecules, digests,
  sites as features, CSV export.
* **GenBank export with two genetic codes** – ORF-finder results can be written into
  a GenBank record alongside the canonical annotation, each feature carrying its own
  `/transl_table`: mitochondrial genes keep the mitochondrial code while
  cytoplasmically-translated ORFs (mitochondrial-derived peptides such as humanin,
  MOTS-c and the SHLPs) are written with the standard code and an explanatory note,
  so any GenBank viewer shows both classes in context. ORFs can be named from
  reference peptides by sequence similarity, and appear as their own tracks
  (one per genetic code) beneath the gene models in the genome view.
* **Analysis** – identity matrix, entropy/identity conservation plots, sequence
  statistics, consensus report, MAFFT alignment (local executable; strategy/threads/--adjustdirection options),
  NCBI BLAST launcher.
* **Formats** – FASTA, Clustal, PHYLIP, NEXUS, Stockholm, GenBank, EMBL, MSF, plus
  full support for BioEdit's own formats: the binary `.bio` project file is read and
  written byte-for-byte as BioEdit 7 writes it, and GenBank can be read/written in
  BioEdit's dialect (top-level `TITLE`, residues in blocks of five).

## Install & run

```bash
pip install -e .            # or: pip install -e ".[dev]" for tests
neoedit examples/cox1_demo.fasta
# or
python -m neoedit
```

### Desktop integration (icon in the Start Menu / taskbar)

```bash
./install_desktop_entry.sh      # needs sudo
```

On WSL this is what makes Windows show the NeoEdit icon instead of the generic WSL
penguin: WSLg only indexes `/usr/share/applications`, so a user-level entry in
`~/.local/share/applications` is ignored. After running it, launch (or pin) NeoEdit
from the Start Menu entry "NeoEdit (Ubuntu)".

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
