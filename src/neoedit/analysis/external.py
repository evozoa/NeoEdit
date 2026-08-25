"""Run MAFFT locally, and launch NCBI BLAST in the browser. No Qt."""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
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


def bundled_mafft_dir() -> str | None:
    """`mafft/` next to the frozen (PyInstaller) app's resources, else None."""
    if not getattr(sys, "frozen", False):
        return None
    base = getattr(sys, "_MEIPASS", None) or os.path.dirname(sys.executable)
    d = os.path.join(base, "mafft")
    return d if os.path.isdir(d) else None


def bundled_mafft() -> str | None:
    """The MAFFT launcher shipped inside a packaged NeoEdit (Windows: mafft.bat from the
    mafft-win package; macOS: the `mafft` script + `libexec/` from Homebrew), else None."""
    d = bundled_mafft_dir()
    if not d:
        return None
    for name in ("mafft.bat", "mafft"):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None


def is_bundled_mafft(exe: str | None) -> bool:
    d = bundled_mafft_dir()
    return bool(exe and d and os.path.abspath(exe).startswith(os.path.abspath(d)))


def _run(cmd: list[str], timeout: float) -> subprocess.CompletedProcess:
    """subprocess.run for MAFFT: no console window in a windowed Windows build, no stdin,
    and MAFFT_BINARIES pointed at the bundled `libexec/` so the `mafft` script finds its
    binaries wherever the app is installed."""
    env = None
    libexec = os.path.join(os.path.dirname(cmd[0]), "libexec")
    if os.path.isdir(libexec):
        env = {**os.environ, "MAFFT_BINARIES": libexec}
    kw = dict(capture_output=True, text=True, timeout=timeout, stdin=subprocess.DEVNULL, env=env)
    if os.name == "nt":
        kw["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    return subprocess.run(cmd, **kw)


def find_mafft(override: str | None = None) -> str | None:
    if override and os.path.exists(override):
        return override
    b = bundled_mafft()
    if b:
        return b
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
        proc = _run([exe, "--version"], timeout=20)
        lines = [ln.strip() for ln in (proc.stderr + proc.stdout).splitlines() if ln.strip()]
        # mafft.bat on Windows prints "Active code page: 65001" before the version line
        return next((ln for ln in lines if ln.lower().startswith("v") or "mafft" in ln.lower()), lines[0])
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
        proc = _run(cmd, timeout=timeout)
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


def mafft_add(existing: list[SequenceRow], new: list[SequenceRow], exe: str | None = None,
              threads: int = 0, adjust_direction: bool = True, timeout: int = 3600) -> list[SequenceRow]:
    """Anchor `new` sequences to an existing alignment with MAFFT --add --keeplength.

    The existing alignment's columns are preserved exactly (insertions unique to the
    added sequences are trimmed), so reference coordinates stay valid. Returns
    existing rows unchanged + the added rows, gapped to the alignment length.
    """
    exe = find_mafft(exe)
    if not exe:
        raise FileNotFoundError("MAFFT executable not found.\n\n" + mafft_install_hint())
    with tempfile.TemporaryDirectory() as td:
        ref = os.path.join(td, "ref.fasta")
        add = os.path.join(td, "add.fasta")
        with open(ref, "w") as fh:
            for i, r in enumerate(existing):
                fh.write(f">e{i}\n{r.seq if any(c in '-.~' for c in r.seq) else r.seq}\n")
        with open(add, "w") as fh:
            for i, r in enumerate(new):
                fh.write(f">a{i}\n{r.ungapped()}\n")
        cmd = [exe, "--quiet", "--thread", str(threads if threads > 0 else -1), "--keeplength"]
        if adjust_direction:
            cmd.append("--adjustdirection")
        cmd += ["--add", add, ref]
        proc = _run(cmd, timeout=timeout)
        if proc.returncode != 0 or not proc.stdout.strip():
            raise RuntimeError(f"MAFFT --add failed (exit {proc.returncode}):\n{proc.stderr[-2000:]}")
    m = mio.loads(proc.stdout, "fasta")
    by_id = {}
    flipped = set()
    for r in m.rows:
        name = r.name
        if name.startswith("_R_"):
            name = name[3:]; flipped.add(name)
        by_id[name] = r.seq
    out = [r.copy() for r in existing]          # columns preserved -> keep originals verbatim
    for i, r in enumerate(new):
        nr = r.copy()
        nr.seq = by_id.get(f"a{i}", r.seq)
        if f"a{i}" in flipped:
            nr.description = (nr.description + " [reverse-complemented by MAFFT]").strip()
        out.append(nr)
    return out


def blast_url(seq: str, program: str = "blastn") -> str:
    s = "".join(c for c in seq if c not in "-.~")
    db = "nt" if program in ("blastn", "tblastx", "tblastn") else "nr"
    return (f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={program}&DATABASE={db}"
            f"&PAGE_TYPE=BlastSearch&QUERY={quote(s)}")


def open_blast(seq: str, program: str = "blastn"):
    webbrowser.open(blast_url(seq, program))
