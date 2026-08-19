"""Run MAFFT locally, and launch NCBI BLAST in the browser. No Qt."""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
import tempfile
import webbrowser
from urllib.parse import quote

from ..model.alignment import SequenceRow
from ..model import io as mio

MAFFT_STRATEGIES = [
    ("Auto (MAFFT chooses by size)", ["--auto"]),
    ("L-INS-i (accurate; <~200 seqs, local homology)", ["--localpair", "--maxiterate", "1000"]),
    ("G-INS-i (accurate; global homology)", ["--globalpair", "--maxiterate", "1000"]),
    ("E-INS-i (long unalignable regions)", ["--genafpair", "--maxiterate", "1000"]),
    ("FFT-NS-i (medium)", ["--retree", "2", "--maxiterate", "1000"]),
    ("FFT-NS-2 (fast; large datasets)", ["--retree", "2", "--maxiterate", "0"]),
]

INSTALL_HINTS = {
    "Windows": "Download the Windows installer from https://mafft.cbrc.jp/alignment/software/windows.html "
               "(or `conda install -c bioconda mafft` if you use conda).",
    "Darwin": "`brew install mafft`  or  `conda install -c bioconda mafft`  "
              "(or the macOS package from https://mafft.cbrc.jp/alignment/software/macportable.html).",
    "Linux": "`sudo apt install mafft`  /  `conda install -c bioconda mafft`  "
             "(or download from https://mafft.cbrc.jp/alignment/software/linux.html).",
}


def mafft_install_hint() -> str:
    return INSTALL_HINTS.get(platform.system(), INSTALL_HINTS["Linux"])


def find_mafft(override: str | None = None) -> str | None:
    if override and os.path.exists(override):
        return override
    for cand in ("mafft", "mafft.bat", "mafft-linsi"):
        p = shutil.which(cand)
        if p:
            return p if not p.endswith("mafft-linsi") else p[: -len("-linsi")]
    # common Windows install location
    for p in (r"C:\Program Files\mafft-win\mafft.bat", r"C:\mafft-win\mafft.bat"):
        if os.path.exists(p):
            return p
    return None


def mafft_version(exe: str) -> str:
    try:
        proc = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=20)
        return (proc.stderr or proc.stdout).strip().splitlines()[0]
    except Exception:
        return "unknown"


def run_mafft(rows: list[SequenceRow], exe: str | None = None, strategy: int = 0, threads: int = 0,
              adjust_direction: bool = False, keep_order: bool = True, extra_args: list[str] | None = None,
              timeout: int = 24 * 3600, seq_type: str | None = None) -> list[SequenceRow]:
    """Align `rows` with MAFFT and return new rows (same order) with gapped sequences."""
    exe = find_mafft(exe)
    if not exe:
        raise FileNotFoundError("MAFFT executable not found.\n\n" + mafft_install_hint()
                                + "\n\nThen set its location in Edit > Preferences if it is not on PATH.")
    args = list(MAFFT_STRATEGIES[strategy][1])
    args += ["--thread", str(threads if threads > 0 else -1)]
    if adjust_direction:
        args.append("--adjustdirection")
    if keep_order:
        args.append("--inputorder")
    if seq_type in ("dna", "rna"):
        args.append("--nuc")
    elif seq_type == "protein":
        args.append("--amino")
    if extra_args:
        args += extra_args
    with tempfile.TemporaryDirectory() as td:
        inp = os.path.join(td, "in.fasta")
        with open(inp, "w") as fh:
            for i, r in enumerate(rows):
                fh.write(f">s{i}\n{r.ungapped()}\n")
        cmd = [exe, "--quiet"] + args + [inp]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0 or not proc.stdout.strip():
            raise RuntimeError(f"MAFFT failed (exit {proc.returncode}):\n{proc.stderr[-2000:]}")
        text = proc.stdout
    m = mio.loads(text, "fasta")
    by_id = {}
    flipped = set()
    for r in m.rows:
        name = r.name
        if name.startswith("_R_"):        # --adjustdirection marks reverse-complemented sequences
            name = name[3:]
            flipped.add(name)
        by_id[name] = r.seq
    result = []
    for i, r in enumerate(rows):
        new = r.copy()
        new.seq = by_id.get(f"s{i}", r.seq)
        if f"s{i}" in flipped and not new.name.endswith("_R_"):
            new.description = (new.description + " [reverse-complemented by MAFFT]").strip()
        result.append(new)
    return result


def blast_url(seq: str, program: str = "blastn") -> str:
    s = "".join(c for c in seq if c not in "-.~")
    db = "nt" if program in ("blastn", "tblastx", "tblastn") else "nr"
    return (f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={program}&DATABASE={db}"
            f"&PAGE_TYPE=BlastSearch&QUERY={quote(s)}")


def open_blast(seq: str, program: str = "blastn"):
    webbrowser.open(blast_url(seq, program))
