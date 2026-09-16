"""Neighbor-joining trees (Saitou & Nei 1987) from evolutionary distances, with bootstrap. No Qt.

Distances follow the usual textbook formulas (as in MEGA and R's ape::dist.dna). Gaps and
ambiguity codes count as missing data: *pairwise deletion* skips a site only for the pairs
that lack it, *complete deletion* drops every site that is missing in any sequence. The
alignment is reduced to unique site patterns first, so a bootstrap replicate is just a new
weight vector over the patterns.
"""
from __future__ import annotations

import csv
import os
import secrets
from collections import Counter
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from ..model.alignment import SequenceRow
from .phylo import tree_labels, write_names

DNA_MODELS = [
    ("p-distance (proportion of differences)", "p"),
    ("Jukes–Cantor (JC69)", "jc69"),
    ("Kimura 2-parameter (K2P)", "k2p"),
    ("Tamura–Nei (TN93)", "tn93"),
]
PROTEIN_MODELS = [
    ("p-distance (proportion of differences)", "p"),
    ("Poisson correction", "poisson"),
]
MODEL_NAMES = {code: label for label, code in DNA_MODELS + PROTEIN_MODELS}
GAP_MODES = [
    ("Pairwise deletion", "pairwise"),
    ("Complete deletion", "complete"),
]
CITATION = ("Saitou N, Nei M (1987) The neighbor-joining method: a new method for reconstructing "
            "phylogenetic trees. Mol Biol Evol 4:406-425.")

_AMINO = "ACDEFGHIKLMNPQRSTVWY"
_DNA_CODES = np.full(256, -1, dtype=np.int8)
for _i, _c in enumerate("ACGT"):
    _DNA_CODES[ord(_c)] = _DNA_CODES[ord(_c.lower())] = _i
_DNA_CODES[ord("U")] = _DNA_CODES[ord("u")] = 3
_AA_CODES = np.full(256, -1, dtype=np.int8)
for _i, _c in enumerate(_AMINO):
    _AA_CODES[ord(_c)] = _AA_CODES[ord(_c.lower())] = _i


class Cancelled(Exception):
    pass


def encode(seqs: list[str], seq_type: str) -> np.ndarray:
    """(n, L) int8 matrix: A,C,G,T/U = 0-3 (or the 20 amino acids = 0-19); anything else -1."""
    table = _AA_CODES if seq_type == "protein" else _DNA_CODES
    raw = np.frombuffer("".join(seqs).encode("ascii", "replace"), dtype=np.uint8)
    return table[raw].reshape(len(seqs), -1)


def site_patterns(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unique columns of X -> (patterns (n, P), counts (P,))."""
    if X.shape[1] == 0:
        return X, np.zeros(0, dtype=np.int64)
    cols, counts = np.unique(X.T, axis=0, return_counts=True)
    return np.ascontiguousarray(cols.T), counts


def pair_counts(X: np.ndarray, W: np.ndarray, dna: bool) -> dict[str, np.ndarray]:
    """Weighted per-pair site counts for the R weight vectors in W (R, P).

    Returns arrays of shape (R, n, n): "valid" (both sites present), "diff", and for DNA
    "ag" (A<->G transitions), "ct" (C<->T transitions) and "tv" (transversions)."""
    n = X.shape[0]
    R = W.shape[0]
    keys = ("valid", "diff", "ag", "ct", "tv") if dna else ("valid", "diff")
    out = {k: np.zeros((R, n, n)) for k in keys}
    Wt = W.T.astype(np.float32)          # counts stay exact below 2**24 sites
    present = X >= 0
    for i in range(n - 1):
        xi, rest = X[i], X[i + 1:]
        valid = present[i] & present[i + 1:]
        diff = valid & (rest != xi)
        masks = {"valid": valid, "diff": diff}
        if dna:
            ag = diff & (((xi == 0) & (rest == 2)) | ((xi == 2) & (rest == 0)))
            ct = diff & (((xi == 1) & (rest == 3)) | ((xi == 3) & (rest == 1)))
            masks.update(ag=ag, ct=ct, tv=diff & ~ag & ~ct)
        for k in keys:
            c = masks[k].astype(np.float32) @ Wt           # (n-i-1, R)
            out[k][:, i, i + 1:] = c.T
            out[k][:, i + 1:, i] = c.T
    return out


def base_frequencies(X: np.ndarray, W: np.ndarray) -> np.ndarray:
    """(R, 4) A/C/G/T frequencies over all sequences, for the R weight vectors."""
    comp = np.stack([(X == b).sum(axis=0) for b in range(4)], axis=1).astype(np.float64)  # (P, 4)
    f = W.astype(np.float64) @ comp
    return f / f.sum(axis=1, keepdims=True)


def _patterns_as_weights(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    P, c = site_patterns(X)
    return P, c[None, :]


def distances(counts: dict[str, np.ndarray], model: str, freqs: np.ndarray | None = None) -> np.ndarray:
    """(R, n, n) distances. NaN where a pair shares no sites; inf where the correction is
    undefined (the sequences are too different for the model)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        valid = counts["valid"]
        p = counts["diff"] / valid
        if model == "p":
            d = p
        elif model == "jc69":
            d = -0.75 * np.log(1 - p * 4 / 3)
        elif model == "poisson":
            d = -np.log(1 - p)
        elif model == "k2p":
            P = (counts["ag"] + counts["ct"]) / valid
            Q = counts["tv"] / valid
            d = -0.5 * np.log(1 - 2 * P - Q) - 0.25 * np.log(1 - 2 * Q)
        elif model == "tn93":
            gA, gC, gG, gT = (freqs[:, k][:, None, None] for k in range(4))
            if np.any(freqs <= 0):
                raise ValueError("Tamura–Nei needs all four bases present in the alignment.")
            gR, gY = gA + gG, gC + gT
            P1 = counts["ag"] / valid
            P2 = counts["ct"] / valid
            Q = counts["tv"] / valid
            d = (-2 * gA * gG / gR * np.log(1 - gR / (2 * gA * gG) * P1 - Q / (2 * gR))
                 - 2 * gC * gT / gY * np.log(1 - gY / (2 * gC * gT) * P2 - Q / (2 * gY))
                 - 2 * (gR * gY - gA * gG * gY / gR - gC * gT * gR / gY) * np.log(1 - Q / (2 * gR * gY)))
        else:
            raise ValueError(f"unknown distance model {model!r}")
    d = np.where(valid > 0, d, np.nan)
    d = np.where(np.isnan(d) & (valid > 0), np.inf, d)        # log of a negative number: saturated
    d = np.where(d < 0, 0.0, d)                              # -0.0 and rounding
    idx = np.arange(d.shape[1])
    d[:, idx, idx] = 0.0
    return d


class Tree:
    """An unrooted tree: leaves are 0..n-1, internal nodes n, n+1, ...; `adj[node]` lists
    (neighbour, branch length)."""

    def __init__(self, n: int):
        self.n = n
        self.adj: dict[int, list[tuple[int, float]]] = {}
        self.center = None

    def add_edge(self, a: int, b: int, length: float):
        self.adj.setdefault(a, []).append((b, length))
        self.adj.setdefault(b, []).append((a, length))

    def _walk(self, root: int):
        """Nodes in DFS preorder from `root`, with each node's parent."""
        order, parent, stack = [], {root: None}, [root]
        while stack:
            u = stack.pop()
            order.append(u)
            for v, _ in reversed(self.adj[u]):
                if v != parent[u]:
                    parent[v] = u
                    stack.append(v)
        return order, parent

    def splits(self) -> set[int]:
        """Non-trivial bipartitions as leaf bitmasks, normalized to the side without leaf 0."""
        full = (1 << self.n) - 1
        order, parent = self._walk(self.center)
        mask: dict[int, int] = {}
        out = set()
        for u in reversed(order):
            m = 1 << u if u < self.n else 0
            for v, _ in self.adj[u]:
                if v != parent[u]:
                    m |= mask[v]
            mask[u] = m
            if parent[u] is not None:
                s = m ^ full if m & 1 else m
                if 2 <= bin(s).count("1") <= self.n - 2:
                    out.add(s)
        return out

    def path_lengths(self) -> np.ndarray:
        """(n, n) leaf-to-leaf path lengths (for tests)."""
        D = np.zeros((self.n, self.n))
        for a in range(self.n):
            dist, stack = {a: 0.0}, [a]
            while stack:
                u = stack.pop()
                for v, w in self.adj[u]:
                    if v not in dist:
                        dist[v] = dist[u] + w
                        stack.append(v)
            D[a] = [dist[b] for b in range(self.n)]
        return D

    def newick(self, labels: list[str], support: dict[int, float] | None = None,
               outgroup: int | None = None) -> str:
        """Newick text, drawn from the last joining node, or from the outgroup's neighbour with
        the outgroup listed first. `labels` must be Newick-safe (phylo.tree_labels);
        `support` maps split bitmasks to percentages."""
        full = (1 << self.n) - 1
        root = self.adj[outgroup][0][0] if outgroup is not None else self.center
        order, parent = self._walk(root)
        length = {}
        for u in order:
            for v, w in self.adj[u]:
                if parent.get(v) == u:
                    length[v] = w
        text: dict[int, str] = {}
        mask: dict[int, int] = {}
        for u in reversed(order):
            kids = [v for v, _ in self.adj[u] if v != parent[u]]
            if outgroup is not None and u == root:
                kids.sort(key=lambda v: v != outgroup)
            if u < self.n:
                s, m = labels[u], 1 << u
            else:
                s = "(" + ",".join(f"{text[v]}:{_fmt(length[v])}" for v in kids) + ")"
                m = 0
                for v in kids:
                    m |= mask[v]
                if support is not None and u != root:
                    key = m ^ full if m & 1 else m
                    if key in support:
                        s += f"{support[key]:.0f}"
            text[u], mask[u] = s, m
        return text[root] + ";"


def _fmt(x: float) -> str:
    return f"{x:.8f}".rstrip("0").rstrip(".") if x else "0"


def neighbor_joining(D: np.ndarray) -> Tree:
    """Saitou & Nei neighbor joining on a full symmetric (n, n) matrix, n >= 3.
    Negative branch-length estimates are set to zero."""
    D = np.array(D, dtype=np.float64)
    n = D.shape[0]
    if n < 3:
        raise ValueError("Neighbor joining needs at least three sequences.")
    tree = Tree(n)
    nodes = list(range(n))
    nxt = n
    while len(nodes) > 3:
        k = len(nodes)
        r = D.sum(axis=1)
        Q = (k - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = divmod(int(np.argmin(Q)), k)
        li = 0.5 * D[i, j] + (r[i] - r[j]) / (2 * (k - 2))
        lj = D[i, j] - li
        u = nxt; nxt += 1
        tree.add_edge(u, nodes[i], max(li, 0.0))
        tree.add_edge(u, nodes[j], max(lj, 0.0))
        du = 0.5 * (D[i] + D[j] - D[i, j])
        D[i, :] = du
        D[:, i] = du
        D[i, i] = 0.0
        D = np.delete(np.delete(D, j, axis=0), j, axis=1)
        nodes[i] = u
        del nodes[j]
    a, b, c = nodes
    center = nxt
    la = 0.5 * (D[0, 1] + D[0, 2] - D[1, 2])
    lb = 0.5 * (D[0, 1] + D[1, 2] - D[0, 2])
    lc = 0.5 * (D[0, 2] + D[1, 2] - D[0, 1])
    for node, w in ((a, la), (b, lb), (c, lc)):
        tree.add_edge(center, node, max(w, 0.0))
    tree.center = center
    return tree


@dataclass
class NJResult:
    tree_path: str
    newick: str
    report_path: str
    distance_path: str
    nsites: int
    seed: int = 0
    replicates: int = 0
    warnings: list[str] = field(default_factory=list)


def _undefined_pairs(d: np.ndarray, names: list[str], what) -> list[str]:
    ii, jj = np.where(np.triu(what(d), 1))
    return [f"{names[i]} – {names[j]}" for i, j in zip(ii[:5], jj[:5])] + (["…"] if len(ii) > 5 else [])


def run_nj(rows: list[SequenceRow], out_dir: str, stem: str, seq_type: str, model: str = "k2p",
           gaps: str = "pairwise", replicates: int = 0, seed: int = 0, outgroup: int | None = None,
           columns: tuple[int, int] | None = None,
           progress: Callable[[int, int], None] | None = None,
           cancelled: Callable[[], bool] | None = None) -> NJResult:
    """Distances + NJ (+ bootstrap) for aligned `rows`; writes the tree, the distance matrix
    and a short report into `out_dir`. `outgroup` is an index into `rows`."""
    protein = seq_type == "protein"
    allowed = dict(PROTEIN_MODELS if protein else DNA_MODELS)
    if model not in allowed.values():
        raise ValueError(f"The {model} distance is not available for {'protein' if protein else 'nucleotide'} data.")
    if len(rows) < 3:
        raise ValueError("Neighbor joining needs at least three sequences.")
    if replicates and len(rows) < 4:
        raise ValueError("Bootstrap support needs at least four sequences.")
    seqs = [r.seq if columns is None else r.seq[columns[0]:columns[1]] for r in rows]
    if len({len(s) for s in seqs}) != 1:
        raise ValueError("The sequences are not aligned (they have different lengths).\n\n"
                         "Align them first (Alignment > Align with MAFFT) or pad them to equal length.")
    names = [r.name for r in rows]
    X = encode(seqs, "protein" if protein else "dna")
    total_sites = X.shape[1]
    # TN93 base frequencies come from every selected column, before any deletion (as ape::dist.dna)
    freqs = base_frequencies(*_patterns_as_weights(X)) if model == "tn93" and not protein else None
    if gaps == "complete":
        X = X[:, (X >= 0).all(axis=0)]
    X, counts = site_patterns(X)
    nsites = int(counts.sum())
    if nsites == 0:
        raise ValueError("No sites left to compare" + (" after complete deletion (every site has a gap or "
                         "ambiguity in some sequence); try pairwise deletion." if gaps == "complete" else "."))

    def dist(W, f=None):
        c = pair_counts(X, W, dna=not protein)
        if model == "tn93" and f is None:
            f = base_frequencies(X, W)          # bootstrap: re-estimated from the resampled sites
        return distances(c, model, f)

    D = dist(counts[None, :], freqs)[0]
    if np.isnan(D).any():
        raise ValueError("Some pairs of sequences share no comparable sites, so their distance is undefined:\n"
                         + "\n".join(_undefined_pairs(D, names, np.isnan)))
    if np.isinf(D).any():
        raise ValueError(f"The {MODEL_NAMES[model]} distance is undefined for some pairs (the sequences are too "
                         "different for this correction):\n" + "\n".join(_undefined_pairs(D, names, np.isinf))
                         + "\n\nUse the p-distance, or remove the most divergent sequences.")
    tree = neighbor_joining(D)

    warnings: list[str] = []
    support = None
    if replicates:
        seed = seed or secrets.randbelow(2**31 - 1) + 1
        rng = np.random.default_rng(seed)
        probs = counts / nsites
        cap = 2 * float(D.max()) if D.max() > 0 else 1.0
        tally: Counter[int] = Counter()
        capped = 0
        n = len(rows)
        batch = max(1, min(25, int(2e8 // (8 * 5 * n * n))))
        done = 0
        while done < replicates:
            if cancelled and cancelled():
                raise Cancelled()
            W = rng.multinomial(nsites, probs, size=min(batch, replicates - done))
            for Db in dist(W):
                if cancelled and cancelled():
                    raise Cancelled()
                bad = ~np.isfinite(Db)
                if bad.any():
                    capped += 1
                    Db = np.where(bad, cap, Db)
                tally.update(neighbor_joining(Db).splits())
                done += 1
                if progress:
                    progress(done, replicates)
        support = {s: 100.0 * tally[s] / replicates for s in tree.splits()}
        if capped:
            warnings.append(f"{capped} bootstrap replicate(s) had undefined distances; those were set to "
                            f"{cap:.4g} (twice the largest observed distance).")

    labels = tree_labels(names)
    newick = tree.newick(labels, support, outgroup)
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, stem)
    tree_path = prefix + ".nj.nwk"
    with open(tree_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(newick + "\n")
    write_names(prefix + ".names.tsv", labels, names)
    distance_path = prefix + ".distances.csv"
    with open(distance_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([""] + labels)
        for label, row in zip(labels, D):
            w.writerow([label] + [f"{x:.6f}" for x in row])
    report_path = prefix + ".nj.txt"
    lines = [
        "Neighbor-joining tree (NeoEdit)",
        "",
        f"Sequences:        {len(rows)} ({'protein' if protein else 'nucleotide'})",
        f"Alignment sites:  {total_sites}" + (f" (columns {columns[0] + 1}-{columns[1]})" if columns else ""),
        f"Sites used:       {nsites}" + (" (complete deletion)" if gaps == "complete" else
                                        " (pairwise deletion: each pair uses the sites both sequences have)"),
        f"Distance:         {MODEL_NAMES[model]}",
        "Gaps/ambiguities: treated as missing data",
        *(["Base frequencies: from all selected columns (bootstrap replicates: from the resampled sites)"]
          if model == "tn93" else []),
        f"Bootstrap:        {replicates} replicates, random seed {seed}" if replicates else "Bootstrap:        none",
        f"Outgroup:         {names[outgroup]}" if outgroup is not None else "Outgroup:         none (unrooted tree)",
        "Negative branch-length estimates are set to zero.",
        "",
        f"Tree:             {os.path.basename(tree_path)}",
        f"Distance matrix:  {os.path.basename(distance_path)}",
        f"Label -> name:    {stem}.names.tsv",
    ]
    lines += [""] + [f"Warning: {w}" for w in warnings] if warnings else []
    lines += ["", "Please cite:", "  " + CITATION]
    if replicates:
        lines.append("  Felsenstein J (1985) Confidence limits on phylogenies: an approach using the "
                     "bootstrap. Evolution 39:783-791.")
    with open(report_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return NJResult(tree_path, newick, report_path, distance_path, nsites,
                    seed if replicates else 0, replicates, warnings)
