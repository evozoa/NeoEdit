"""Phylogenetic trees for the tree viewer: Newick/NEXUS reading and writing, rerooting,
ladderizing and layout. No Qt.

A support value (bootstrap, "53.5/59", ...) written as an internal node label belongs to the
branch above that node. Rerooting keeps it on that branch: the tree is turned into an
unrooted graph whose edges carry the lengths and labels, and rebuilt from the new root.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

_SPECIAL = set("()[]':;,")


@dataclass(eq=False)
class Node:
    name: str = ""                 # tip name (internal nodes: usually empty)
    length: float | None = None    # length of the branch above this node
    label: str = ""                # internal node label: support for the branch above
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = field(default=None, repr=False)

    @property
    def is_tip(self) -> bool:
        return not self.children

    def add(self, child: "Node") -> "Node":
        child.parent = self
        self.children.append(child)
        return child

    def walk(self):
        """Nodes in preorder (iterative, safe for deep trees)."""
        stack = [self]
        while stack:
            n = stack.pop()
            yield n
            stack.extend(reversed(n.children))

    def postorder(self):
        return list(reversed(list(self._rev_pre())))

    def _rev_pre(self):
        stack = [self]
        while stack:
            n = stack.pop()
            yield n
            stack.extend(n.children)

    def tips(self) -> list["Node"]:
        return [n for n in self.walk() if n.is_tip]


# ---------------------------------------------------------------- reading
class TreeError(ValueError):
    pass


def _tokens(text: str):
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif c == "[":                                   # comment, e.g. [&R] or [&label=...]
            j = text.find("]", i)
            if j < 0:
                raise TreeError("unterminated [comment] in tree")
            i = j + 1
        elif c in "(),:;":
            yield c
            i += 1
        elif c == "'":
            buf, i = [], i + 1
            while True:
                if i >= n:
                    raise TreeError("unterminated quoted label in tree")
                if text[i] == "'":
                    if i + 1 < n and text[i + 1] == "'":
                        buf.append("'"); i += 2; continue
                    i += 1
                    break
                buf.append(text[i]); i += 1
            yield ("Q", "".join(buf))
        else:
            j = i
            while j < n and not text[j].isspace() and text[j] not in "(),:;[":
                j += 1
            yield ("W", text[i:j])
            i = j


def parse_newick(text: str) -> Node:
    """The first tree in `text`."""
    toks = list(_tokens(text))
    if not toks:
        raise TreeError("no tree found")
    root = cur = Node()
    stack: list[Node] = []
    expect_len = False
    after_close = False
    for t in toks:
        if t == "(":
            stack.append(cur)
            cur = cur.add(Node())
            after_close = False
        elif t == ",":
            if not stack:
                raise TreeError("unexpected ',' in tree")
            cur = stack[-1].add(Node())
            after_close = False
        elif t == ")":
            if not stack:
                raise TreeError("unbalanced ')' in tree")
            cur = stack.pop()
            after_close = True
        elif t == ":":
            expect_len = True
        elif t == ";":
            break
        else:
            kind, word = t
            if expect_len:
                try:
                    cur.length = float(word)
                except ValueError:
                    raise TreeError(f"bad branch length {word!r}") from None
                expect_len = False
            elif after_close:
                cur.label = word
            else:
                cur.name = word
    if stack:
        raise TreeError("unbalanced '(' in tree")
    if not root.children:
        raise TreeError("not a Newick tree (no parentheses)")
    return root


def _nexus_trees(text: str) -> tuple[str, dict[str, str]]:
    """First tree string and the translate table of a NEXUS TREES block."""
    m = re.search(r"begin\s+trees\s*;(.*?)end\s*;", text, re.I | re.S)
    if not m:
        raise TreeError("NEXUS file has no TREES block")
    block = m.group(1)
    trans: dict[str, str] = {}
    tm = re.search(r"translate\s+(.*?);", block, re.I | re.S)
    if tm:
        for entry in re.split(r",\s*(?=(?:[^']*'[^']*')*[^']*$)", tm.group(1)):
            parts = entry.strip().split(None, 1)
            if len(parts) == 2:
                trans[parts[0]] = parts[1].strip().strip("'").replace("''", "'")
    t = re.search(r"tree\s+[^=]+=\s*(.*?;)", block, re.I | re.S)
    if not t:
        raise TreeError("NEXUS TREES block has no tree")
    return t.group(1), trans


def read_tree(path: str) -> Node:
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    if text.lstrip()[:6].upper() == "#NEXUS":
        newick, trans = _nexus_trees(text)
        root = parse_newick(newick)
        for tip in root.tips():
            tip.name = trans.get(tip.name, tip.name)
        return root
    return parse_newick(text)


def names_file_for(tree_path: str) -> str | None:
    """The `<stem>.names.tsv` (label -> full name) NeoEdit wrote next to a tree, if any."""
    d, base = os.path.split(os.path.abspath(tree_path))
    try:
        cands = [f for f in os.listdir(d) if f.endswith(".names.tsv")]
    except OSError:
        return None
    best = [f for f in cands if base.startswith(f[:-len(".names.tsv")] + ".")]
    pick = max(best, key=len) if best else (cands[0] if len(cands) == 1 else None)
    return os.path.join(d, pick) if pick else None


def read_names(path: str) -> dict[str, str]:
    out = {}
    with open(path, encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t", 1)
            if i == 0 and parts[0] in ("label", "id"):
                continue
            if len(parts) == 2:
                out[parts[0]] = parts[1]
    return out


# ---------------------------------------------------------------- writing
def _quote(s: str) -> str:
    if s and not any(c in _SPECIAL or c.isspace() for c in s):
        return s
    return "'" + s.replace("'", "''") + "'"


def _fmt_len(x: float) -> str:
    return f"{x:.10f}".rstrip("0").rstrip(".") if x else "0"


def to_newick(root: Node) -> str:
    out: dict[int, str] = {}
    for n in root.postorder():
        if n.children:
            s = "(" + ",".join(out[id(c)] for c in n.children) + ")" + (_quote(n.label) if n.label else "")
        else:
            s = _quote(n.name)
        if n.length is not None and n is not root:
            s += ":" + _fmt_len(n.length)
        out[id(n)] = s
    return out[id(root)] + ";"


# ---------------------------------------------------------------- rerooting
def _unrooted(root: Node):
    """Undirected graph {node: {neighbour: (length, label)}}, with degree-2 nodes (a bifurcating
    root) merged away."""
    adj: dict[Node, dict[Node, tuple[float, str]]] = {}
    for n in root.walk():
        adj.setdefault(n, {})
        for c in n.children:
            e = (c.length or 0.0, c.label if c.children else "")
            adj[n][c] = e
            adj.setdefault(c, {})[n] = e
    for n in list(adj):
        if len(adj[n]) == 2 and not n.name:
            (a, (la, sa)), (b, (lb, sb)) = adj[n].items()
            e = (la + lb, sa or sb)
            del adj[a][n], adj[b][n], adj[n]
            adj[a][b] = e
            adj[b][a] = e
    return adj


def _build(adj, node: Node, parent: Node | None) -> Node:
    """A rooted copy of the subtree at `node`, directed away from `parent`."""
    copies: dict[Node, Node] = {}
    stack = [(node, parent)]
    order = []
    while stack:
        u, p = stack.pop()
        copies[u] = Node(name=u.name)
        order.append((u, p))
        for v in reversed(list(adj[u])):
            if v is not p:
                stack.append((v, u))
    for u, p in order:
        if p is not None and u is not node:
            length, label = adj[p][u]
            cu = copies[u]
            cu.length = length
            cu.label = label if len(adj[u]) > 1 else ""
            copies[p].add(cu)
    return copies[node]


def reroot(root: Node, node: Node, position: float | None = None) -> Node:
    """A new tree rooted on the branch above `node`, `position` from `node` (default: middle)."""
    if node.parent is None:
        return root
    adj = _unrooted(root)
    a = node
    # the branch above `node` may have been merged (node's parent was a degree-2 root)
    if node.parent in adj[a]:
        b = node.parent
    else:
        b = next(v for v in adj[a] if v not in node.children)
    length, label = adj[a][b]
    x = length / 2 if position is None else max(0.0, min(length, position))
    new = Node()
    left = _build(adj, a, b)
    right = _build(adj, b, a)
    left.length, right.length = x, length - x
    if len(adj[a]) > 1:
        left.label = label
    if len(adj[b]) > 1:
        right.label = label
    new.add(left)
    new.add(right)
    return new


def _tip_distances(adj, start: Node) -> dict[Node, tuple[float, Node | None]]:
    """Distance from `start` to every node, with the previous node on the path."""
    dist = {start: (0.0, None)}
    stack = [start]
    while stack:
        u = stack.pop()
        for v, (length, _) in adj[u].items():
            if v not in dist:
                dist[v] = (dist[u][0] + length, u)
                stack.append(v)
    return dist


def midpoint_root(root: Node) -> Node:
    """Root at the middle of the longest tip-to-tip path."""
    adj = _unrooted(root)
    tips = [n for n in adj if len(adj[n]) <= 1]
    if len(tips) < 3:
        return root
    d0 = _tip_distances(adj, tips[0])
    far1 = max(tips, key=lambda t: d0[t][0])
    d1 = _tip_distances(adj, far1)
    far2 = max(tips, key=lambda t: d1[t][0])
    half = d1[far2][0] / 2
    # walk back from far2 towards far1 until the midpoint is passed
    u = far2
    while True:
        prev = d1[u][1]
        if d1[prev][0] <= half:
            break
        u = prev
    edge_len = adj[u][prev][0]
    new = Node()
    left = _build(adj, u, prev)
    right = _build(adj, prev, u)
    left.length = d1[u][0] - half
    right.length = edge_len - left.length
    label = adj[u][prev][1]
    if len(adj[u]) > 1:
        left.label = label
    if len(adj[prev]) > 1:
        right.label = label
    new.add(left)
    new.add(right)
    return new


def ladderize(root: Node, big_last: bool = True):
    """Order children by clade size, in place."""
    size = {}
    for n in root.postorder():
        size[id(n)] = 1 if n.is_tip else sum(size[id(c)] for c in n.children)
    for n in root.walk():
        n.children.sort(key=lambda c: size[id(c)], reverse=not big_last)


def rotate(node: Node):
    node.children.reverse()


def mrca(nodes: list[Node]) -> Node | None:
    if not nodes:
        return None
    paths = []
    for n in nodes:
        p, x = [], n
        while x is not None:
            p.append(x); x = x.parent
        paths.append(p[::-1])
    common = None
    for group in zip(*paths):
        if all(g is group[0] for g in group):
            common = group[0]
        else:
            break
    return common


# ---------------------------------------------------------------- layout
@dataclass
class Layout:
    x: dict[int, float]         # id(node) -> distance from the root (or depth for a cladogram)
    y: dict[int, float]         # id(node) -> row (tips 0, 1, 2, ...)
    width: float                # largest x
    ntips: int
    has_lengths: bool


def layout(root: Node, cladogram: bool = False) -> Layout:
    has_lengths = any(n.length for n in root.walk() if n is not root)
    x: dict[int, float] = {}
    y: dict[int, float] = {}
    if cladogram or not has_lengths:
        height = {}
        for n in root.postorder():
            height[id(n)] = 0 if n.is_tip else 1 + max(height[id(c)] for c in n.children)
        top = height[id(root)]
        for n in root.walk():
            x[id(n)] = top - height[id(n)]
    else:
        for n in root.walk():
            x[id(n)] = 0.0 if n.parent is None else x[id(n.parent)] + max(0.0, n.length or 0.0)
    row = 0
    for n in root.postorder():
        if n.is_tip:
            y[id(n)] = row
            row += 1
        else:
            y[id(n)] = (y[id(n.children[0])] + y[id(n.children[-1])]) / 2
    return Layout(x, y, max(x.values()) if x else 0.0, row, has_lengths and not cladogram)
