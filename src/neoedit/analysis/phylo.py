"""Maximum-likelihood trees with IQ-TREE 3 (an external program). No Qt.

A run lives in its own output folder: the alignment is written there under safe ids
(s0, s1, ...) because IQ-TREE rejects names that collide once it has sanitized them, IQ-TREE
writes its usual files next to it, and `finish_run` writes `<stem>.nwk`, the ML tree with the
real sequence names back in, for FigTree / iTOL / R.
"""
from __future__ import annotations

import os
import platform
import re
import shutil
import subprocess
from dataclasses import dataclass, field

from ..model.alignment import SequenceRow
from .external import bundled_tool_dir

IQTREE_SITE = "https://iqtree.github.io/"

# (label, -m value); the first entry is the default for that data type
MODEL_PRESETS = {
    "dna": [
        ("ModelFinder: choose the best-fit model (BIC)", "MFP"),
        ("GTR+F+I+G4", "GTR+F+I+G4"),
        ("GTR+F+G4", "GTR+F+G4"),
        ("HKY+F+G4", "HKY+F+G4"),
        ("K2P+G4", "K2P+G4"),
        ("JC (Jukes-Cantor)", "JC"),
    ],
    "protein": [
        ("ModelFinder: choose the best-fit model (BIC)", "MFP"),
        ("LG+G4", "LG+G4"),
        ("WAG+G4", "WAG+G4"),
        ("JTT+G4", "JTT+G4"),
        ("mtVer+F+G4 (vertebrate mitochondrial)", "mtVer+F+G4"),
        ("mtREV+F+G4 (mitochondrial)", "mtREV+F+G4"),
    ],
}

# (label, flag, default replicates, minimum replicates)
BOOTSTRAP_METHODS = [
    ("None", None, 0, 0),
    ("Ultrafast bootstrap (UFBoot)", "-B", 1000, 1000),
    ("Standard nonparametric bootstrap (slow)", "-b", 100, 1),
]

INSTALL_HINTS = {
    "Windows": f"Download the Windows zip from {IQTREE_SITE} and point Edit > Preferences at bin\\iqtree3.exe.",
    "Darwin": f"`conda install -c bioconda iqtree`, or download the macOS zip from {IQTREE_SITE}.",
    "Linux": f"`conda install -c bioconda iqtree`, or download the Linux tar.gz from {IQTREE_SITE}.",
}

_NEWICK_SPECIAL = re.compile(r"[\s()\[\]':;,]")
_SAFE_LABEL = re.compile(r"(?<=[(,])s(\d+)(?=[:,)])")


def iqtree_install_hint() -> str:
    return INSTALL_HINTS.get(platform.system(), INSTALL_HINTS["Linux"])


def bundled_iqtree() -> str | None:
    """The IQ-TREE executable shipped inside a packaged NeoEdit, else None."""
    d = bundled_tool_dir("iqtree")
    if not d:
        return None
    for name in ("iqtree3.exe", "iqtree3"):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None


def find_iqtree(override: str | None = None) -> str | None:
    if override and os.path.exists(override):
        return override
    b = bundled_iqtree()
    if b:
        return b
    for cand in ("iqtree3", "iqtree2", "iqtree"):
        p = shutil.which(cand)
        if p:
            return p
    return None


# The official static builds carry LLVM's OpenMP runtime, which prints "libarcher.so: cannot open
# shared object file" on every run unless its tool interface is switched off.
IQTREE_ENV = {"OMP_TOOL": "disabled"}


def auto_threads_max() -> int:
    """Upper bound for `-T AUTO`: all cores but one, so NeoEdit and the rest of the machine stay responsive."""
    return max(1, (os.cpu_count() or 1) - 1)


def _popen_kw() -> dict:
    kw = dict(stdin=subprocess.DEVNULL, env={**os.environ, **IQTREE_ENV})
    if os.name == "nt":
        kw["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    return kw


def iqtree_version(exe: str) -> str:
    """e.g. "3.1.4", or "unknown" if `exe` does not run."""
    try:
        proc = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=30, **_popen_kw())
    except Exception:
        return "unknown"
    m = re.search(r"IQ-TREE(?: multicore)? version (\S+)", proc.stdout + proc.stderr)
    return m.group(1) if m else "unknown"


def version_problem(version: str) -> str | None:
    """Why this IQ-TREE can't be used, or None. IQ-TREE 1.x has a different command line."""
    m = re.match(r"(\d+)", version)
    if not m:
        return "could not determine the IQ-TREE version"
    if int(m.group(1)) < 2:
        return f"IQ-TREE {version} is too old (NeoEdit needs IQ-TREE 2 or 3)"
    return None


def newick_label(name: str) -> str:
    """A sequence name as a Newick taxon label, single-quoted when it needs to be."""
    if name and not _NEWICK_SPECIAL.search(name):
        return name
    return "'" + name.replace("'", "''") + "'"


def rename_newick(text: str, names: list[str]) -> str:
    """Replace the safe ids (s0, s1, ...) that IQ-TREE saw with the real names."""
    def sub(m):
        i = int(m.group(1))
        return newick_label(names[i]) if i < len(names) else m.group(0)
    return _SAFE_LABEL.sub(sub, text)


def unique_dir(path: str) -> str:
    """`path`, or `path_2`, `path_3`, ... whichever does not exist yet."""
    if not os.path.exists(path):
        return path
    n = 2
    while os.path.exists(f"{path}_{n}"):
        n += 1
    return f"{path}_{n}"


def safe_stem(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._") or "alignment"


@dataclass
class IQTreeRun:
    out_dir: str
    stem: str
    names: list[str]
    seq_type: str
    nsites: int

    @property
    def prefix(self) -> str:
        return os.path.join(self.out_dir, self.stem)

    @property
    def input_path(self) -> str:
        return self.prefix + ".input.fasta"

    @property
    def tree_path(self) -> str:
        return self.prefix + ".nwk"


@dataclass
class TreeResult:
    tree_path: str
    newick: str
    report_path: str
    best_model: str = ""
    log_likelihood: str = ""
    files: list[str] = field(default_factory=list)


def prepare_run(rows: list[SequenceRow], out_dir: str, stem: str, seq_type: str,
                columns: tuple[int, int] | None = None) -> IQTreeRun:
    """Write the (aligned) `rows` to a new folder as IQ-TREE input.

    `columns` = (start, end), end exclusive, restricts the run to those alignment columns."""
    if len(rows) < 3:
        raise ValueError("IQ-TREE needs at least three sequences.")
    seqs = [r.seq if columns is None else r.seq[columns[0]:columns[1]] for r in rows]
    if len({len(s) for s in seqs}) != 1:
        raise ValueError("The sequences are not aligned (they have different lengths).\n\n"
                         "Align them first (Alignment > Align with MAFFT) or pad them to equal length.")
    if not seqs[0]:
        raise ValueError("No alignment columns to analyse.")
    os.makedirs(out_dir, exist_ok=True)
    run = IQTreeRun(out_dir, stem, [r.name for r in rows],
                    "protein" if seq_type == "protein" else "dna", len(seqs[0]))
    gap = str.maketrans(".~", "--")
    with open(run.input_path, "w", encoding="utf-8", newline="\n") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">s{i}\n{s.translate(gap)}\n")
    with open(run.prefix + ".names.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("id\tname\n")
        for i, name in enumerate(run.names):
            fh.write(f"s{i}\t{name}\n")
    return run


def iqtree_args(run: IQTreeRun, model: str = "MFP", bootstrap: int = 0, replicates: int = 0,
                alrt: int = 0, threads: int = 0, seed: int = 0, outgroup: list[int] | None = None,
                extra_args: list[str] | None = None) -> list[str]:
    """Command-line arguments (without the executable). `bootstrap` indexes BOOTSTRAP_METHODS;
    `threads` 0 = let IQ-TREE decide, up to auto_threads_max(); `seed` 0 = random; `outgroup` = row indexes into run.names."""
    args = ["-s", run.input_path, "--prefix", run.prefix,
            "-st", "AA" if run.seq_type == "protein" else "DNA",
            "-m", model.strip() or "MFP",
            "-T", str(threads) if threads > 0 else "AUTO"]
    if threads <= 0:
        args += ["--threads-max", str(auto_threads_max())]
    _, flag, _, minimum = BOOTSTRAP_METHODS[bootstrap]
    if flag:
        args += [flag, str(max(replicates, minimum))]
    if alrt > 0:
        args += ["--alrt", str(alrt)]
    if seed > 0:
        args += ["--seed", str(seed)]
    if outgroup:
        args += ["-o", ",".join(f"s{i}" for i in outgroup)]
    if extra_args:
        args += extra_args
    return args


def _grep(text: str, pattern: str) -> str:
    m = re.search(pattern, text, re.M)
    return m.group(1).strip() if m else ""


def finish_run(run: IQTreeRun) -> TreeResult:
    """After IQ-TREE exits successfully: write `<stem>.nwk` with the real names and summarize."""
    treefile = run.prefix + ".treefile"
    if not os.path.exists(treefile):
        raise RuntimeError(f"IQ-TREE did not write {treefile}; see {run.prefix}.log")
    with open(treefile, encoding="utf-8") as fh:
        newick = rename_newick(fh.read().strip(), run.names)
    with open(run.tree_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(newick + "\n")
    contree = run.prefix + ".contree"
    if os.path.exists(contree):
        with open(contree, encoding="utf-8") as fh:
            con = rename_newick(fh.read().strip(), run.names)
        with open(run.prefix + ".consensus.nwk", "w", encoding="utf-8", newline="\n") as fh:
            fh.write(con + "\n")
    report = run.prefix + ".iqtree"
    text = ""
    if os.path.exists(report):
        with open(report, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    return TreeResult(
        tree_path=run.tree_path, newick=newick, report_path=report,
        best_model=_grep(text, r"^Best-fit model according to \w+:\s*(\S+)")
        or _grep(text, r"^Model of substitution:\s*(\S+)"),
        log_likelihood=_grep(text, r"^Log-likelihood of the tree:\s*(\S+)"),
        files=sorted(os.listdir(run.out_dir)))


def run_iqtree(run: IQTreeRun, exe: str | None = None, timeout: float = 24 * 3600, **opts) -> TreeResult:
    """Run IQ-TREE to completion (blocking). The GUI runs the same command asynchronously."""
    exe = find_iqtree(exe)
    if not exe:
        raise FileNotFoundError("IQ-TREE executable not found.\n\n" + iqtree_install_hint())
    proc = subprocess.run([exe] + iqtree_args(run, **opts), capture_output=True, text=True,
                          timeout=timeout, cwd=run.out_dir, **_popen_kw())
    if proc.returncode != 0:
        raise RuntimeError(f"IQ-TREE failed (exit {proc.returncode}):\n{(proc.stdout + proc.stderr)[-2000:]}")
    return finish_run(run)
