## What's new in 0.3.0

- **Phylogenetic trees** — *Analysis ▸ Phylogeny*:
  - **Maximum likelihood (IQ-TREE)…** (Ctrl+Shift+Y). IQ-TREE 3 is included: automatic model selection or a model of your choice, bootstrap and SH-aLRT branch support, optional outgroup.
  - **Neighbor joining…** with p-distance, Jukes–Cantor, Kimura 2-parameter or Tamura–Nei distances and bootstrap support; also saves the distance matrix for Excel.
  - Each run saves its results in a new folder next to your alignment. Open the tree (`.nwk`) in [FigTree](https://github.com/rambaut/figtree/releases) or [iTOL](https://itol.embl.de/) — or use **View in iTOL…**.
- **Your choice of names for NCBI imports** — *Sequence names…* on the NCBI tab (or *File ▸ NCBI sequence names…*): drag fields such as accession, genus, species, isolate, voucher or country into the order you want.
- **Edit ▸ Select to Beginning / Select to End**. In Edit mode, **Delete** and **Backspace** now remove the selected residues.
- Clicking a gene or exon in the gene-model view selects it in the alignment.
- Dragging sequence names works as in BioEdit: click a name to select it, click it again to drag it to a new place.

## Install NeoEdit

Download the file for your computer, then follow the one-time steps below.

| Computer | Download |
|---|---|
| **Windows** (10 / 11) | `NeoEdit-Windows-Setup.exe` |
| **Mac with Apple Silicon** (M1, M2, M3, M4 — most Macs from 2021 on) | `NeoEdit-macOS-AppleSilicon.dmg` |
| **Mac with an Intel chip** (older Macs) | `NeoEdit-macOS-Intel.dmg` |
| **Linux** (Ubuntu 22.04 or newer, Mint, Debian) | `NeoEdit-Linux-x86_64.deb` |

Not sure which Mac you have? Click the Apple menu ▸ **About This Mac**: the "Chip" line says *Apple M…* (Apple Silicon) or *Intel*.

### Windows

1. Open the downloaded `NeoEdit-Windows-Setup.exe`.
2. If a blue **"Windows protected your PC"** box appears, click **More info**, then **Run anyway**. (NeoEdit is open-source but not yet code-signed, so Windows shows this once.)
3. Click **Next** through the installer. You do not need an administrator password.
4. NeoEdit is in the Start Menu. The installer can also put an icon on the desktop.

### Mac

1. Open the downloaded `.dmg` and drag **NeoEdit** onto the **Applications** folder.
2. Open **Applications** and double-click **NeoEdit**. macOS will say *"Apple could not verify NeoEdit is free of malware"* — click **Done** (not "Move to Trash").
3. Open **System Settings ▸ Privacy & Security**, scroll down to the *Security* section, and click **Open Anyway** next to the NeoEdit message. Enter your password and click **Open**.
4. That's it — from now on NeoEdit opens normally.

### Linux

1. Double-click the downloaded `NeoEdit-Linux-x86_64.deb` and click **Install** in the App Center, or run `sudo apt install ./NeoEdit-Linux-x86_64.deb` in a terminal opened in *Downloads*.
2. Press the Super key and type **NeoEdit**. Remove later with `sudo apt remove neoedit`.

MAFFT (for *Alignment ▸ Align with MAFFT*) and IQ-TREE 3 (for *Analysis ▸ Phylogeny ▸ Maximum likelihood*) are included. Everything works offline except *File ▸ Import from NCBI / Ensembl / UCSC*, which needs an internet connection.

**Problems?** Open an issue at https://github.com/evozoa/NeoEdit/issues.
