"""Run external aligners (MAFFT / MUSCLE / Clustal Omega) and BLAST. No Qt."""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import webbrowser
from urllib.parse import quote

from ..model.alignment import AlignmentModel, SequenceRow
from ..model import io as mio

ALIGNERS = {
    "MAFFT": {"exe": ["mafft"], "cmd": lambda exe, inp, out: [exe, "--auto", "--quiet", inp], "stdout": True},
    "MUSCLE (v5)": {"exe": ["muscle", "muscle5"], "cmd": lambda exe, inp, out: [exe, "-align", inp, "-output", out], "stdout": False},
    "MUSCLE (v3)": {"exe": ["muscle3", "muscle"], "cmd": lambda exe, inp, out: [exe, "-in", inp, "-out", out, "-quiet"], "stdout": False},
    "Clustal Omega": {"exe": ["clustalo"], "cmd": lambda exe, inp, out: [exe, "-i", inp, "-o", out, "--force", "--outfmt=fasta"], "stdout": False},
}


def find_executable(name: str, override: str | None = None) -> str | None:
    if override and os.path.exists(override):
        return override
    for cand in ALIGNERS[name]["exe"]:
        p = shutil.which(cand)
        if p:
            return p
    return None


def run_aligner(name: str, rows: list[SequenceRow], exe: str | None = None, extra_args: list[str] | None = None,
                timeout: int = 3600) -> list[SequenceRow]:
    spec = ALIGNERS[name]
    exe = find_executable(name, exe)
    if not exe:
        raise FileNotFoundError(f"{name} executable not found on PATH. Set it in Preferences.")
    with tempfile.TemporaryDirectory() as td:
        inp = os.path.join(td, "in.fasta")
        out = os.path.join(td, "out.fasta")
        # write ungapped, with safe numeric ids to survive aligners' name mangling
        with open(inp, "w") as fh:
            for i, r in enumerate(rows):
                fh.write(f">s{i}\n{r.ungapped()}\n")
        cmd = spec["cmd"](exe, inp, out)
        if extra_args:
            cmd[1:1] = extra_args
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0:
            raise RuntimeError(f"{name} failed (exit {proc.returncode}):\n{proc.stderr[-2000:]}")
        text = proc.stdout if spec["stdout"] else open(out).read()
    m = mio.loads(text, "fasta")
    by_id = {r.name: r.seq for r in m.rows}
    result = []
    for i, r in enumerate(rows):
        new = r.copy()
        new.seq = by_id.get(f"s{i}", r.seq)
        result.append(new)
    return result


def blast_url(seq: str, program: str = "blastn") -> str:
    s = "".join(c for c in seq if c not in "-.~")
    db = "nt" if program in ("blastn", "tblastx", "tblastn") else "nr"
    return (f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={program}&DATABASE={db}"
            f"&PAGE_TYPE=BlastSearch&QUERY={quote(s)}")


def open_blast(seq: str, program: str = "blastn"):
    webbrowser.open(blast_url(seq, program))
