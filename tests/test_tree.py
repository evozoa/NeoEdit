"""Tree viewer model: Newick/NEXUS parsing, rerooting that keeps supports on their branches."""
import io
import os
import random

import pytest

from neoedit.analysis import tree as T

TREE = "((a:1,b:2)90:0.5,(c:1,'d e':1)[&note]80:0.5,f:3);"


def splits(root):
    """{frozenset(tips on the side without the first tip, sorted): support label}"""
    names = sorted(t.name for t in root.tips())
    ref = names[0]
    out = {}
    for n in root.walk():
        if n is root or n.is_tip:
            continue
        side = {t.name for t in n.tips()}
        if ref in side:
            side = set(names) - side
        if 2 <= len(side) <= len(names) - 2:
            key = frozenset(side)
            out[key] = out.get(key) or n.label
    return out


def path_lengths(root):
    """Tip-to-tip distances (unchanged by rerooting)."""
    tips = root.tips()
    dist = {}
    for t in tips:
        d, x = 0.0, t
        up = {}
        while x is not None:
            up[id(x)] = d
            d += x.length or 0.0
            x = x.parent
        for u in tips:
            d2, y = 0.0, u
            while id(y) not in up:
                d2 += y.length or 0.0
                y = y.parent
            dist[(t.name, u.name)] = up[id(y)] + d2
    return dist


def same_distances(a, b):
    return a.keys() == b.keys() and all(abs(a[k] - b[k]) < 1e-9 for k in a)


def test_parse_and_write():
    r = T.parse_newick(TREE)
    assert [t.name for t in r.tips()] == ["a", "b", "c", "d e", "f"]
    assert [n.label for n in r.children] == ["90", "80", ""]
    assert r.children[2].length == 3
    assert T.to_newick(r) == "((a:1,b:2)90:0.5,(c:1,'d e':1)80:0.5,f:3);"
    assert T.to_newick(T.parse_newick(T.to_newick(r))) == T.to_newick(r)
    iq = T.parse_newick("(A:0.1,(B:0.2,C:0.3)53.5/59:0.05,D:0);")
    assert iq.children[1].label == "53.5/59"
    assert T.parse_newick("((x,y),z);").tips()[0].length is None
    assert T.parse_newick("(it''s:1,'q''t':2,w:3);").tips()[1].name == "q't"
    for bad in ("(a,b", "a,b);", "(a:x,b);", "abc;", ""):
        with pytest.raises(T.TreeError):
            T.parse_newick(bad)


def test_matches_biopython_on_real_output():
    from Bio import Phylo
    text = ("(NC_005129.2_Elephas:0.0266,((NC_073531.1_Rhynchocyon_s:0.0392,NC_082504.1_Rhynchocyon_p:0.0305)"
            "53.5/59:0.0105,NC_073534.1_Rhynchocyon_c:0.0256)100/100:0.2954,NC_1_x:0.1);")
    ours = T.parse_newick(text)
    bio = Phylo.read(io.StringIO(text), "newick")
    assert [t.name for t in ours.tips()] == [t.name for t in bio.get_terminals()]
    assert [round(t.length, 6) for t in ours.tips()] == [round(t.branch_length, 6) for t in bio.get_terminals()]


def test_reroot_keeps_supports_and_distances():
    r = T.parse_newick(TREE)
    before_s, before_d = splits(r), path_lengths(r)
    for i, node in enumerate(n for n in list(r.walk()) if n is not r):
        new = T.reroot(r, node)
        assert splits(new) == before_s, T.to_newick(new)
        assert same_distances(path_lengths(new), before_d)
        assert len(new.children) == 2
        assert sorted(t.name for t in new.children[0].tips()) == sorted(t.name for t in node.tips())
    a = r.tips()[0]
    assert T.to_newick(T.reroot(r, a)) == "(a:0.5,(((c:1,'d e':1)80:0.5,f:3)90:0.5,b:2):0.5);"
    assert T.to_newick(T.reroot(r, a, position=0.2)).startswith("(a:0.2,")
    # rerooting twice, and on a bifurcating root whose two branches are one edge
    twice = T.reroot(r, a)
    f = next(t for t in twice.tips() if t.name == "f")
    assert splits(T.reroot(twice, f)) == before_s
    rb = T.parse_newick("((a:1,b:1)70:0.5,(c:1,d:1)70:0.5);")
    assert T.to_newick(T.reroot(rb, rb.tips()[0])) == "(a:0.5,(b:1,(c:1,d:1)70:1):0.5);"


def test_random_trees_reroot():
    rng = random.Random(3)
    for trial in range(30):
        n = rng.randint(4, 25)
        nodes = [T.Node(name=f"t{i}", length=rng.uniform(0.01, 1)) for i in range(n)]
        while len(nodes) > 3:
            a, b = rng.sample(nodes, 2)
            nodes.remove(a); nodes.remove(b)
            p = T.Node(length=rng.uniform(0.01, 1), label=str(rng.randint(0, 100)))
            p.add(a); p.add(b)
            nodes.append(p)
        root = T.Node()
        for x in nodes:
            root.add(x)
        s0, d0 = splits(root), path_lengths(root)
        cur = root
        for _ in range(4):
            target = rng.choice([x for x in cur.walk() if x is not cur])
            cur = T.reroot(cur, target)
            assert splits(cur) == s0
            assert same_distances(path_lengths(cur), d0)
        mid = T.midpoint_root(cur)
        assert splits(mid) == s0
        assert same_distances(path_lengths(mid), d0)
        depths = [sum_to_root(t) for t in mid.tips()]
        left = max(sum_to_root(t) for t in mid.children[0].tips())
        right = max(sum_to_root(t) for t in mid.children[1].tips())
        assert left == pytest.approx(right)                        # the midpoint balances the longest path
        assert max(depths) == pytest.approx(max(d0.values()) / 2)


def sum_to_root(t):
    d = 0.0
    while t.parent is not None:
        d += t.length or 0.0
        t = t.parent
    return d


def test_ladderize_rotate_mrca_layout():
    r = T.parse_newick("(((a:1,b:1):1,c:1):1,d:1,(e:1,(f:1,(g:1,h:1):1):1):1);")
    T.ladderize(r)
    assert [t.name for t in r.tips()] == ["d", "c", "a", "b", "e", "f", "g", "h"]
    T.rotate(r)
    assert [t.name for t in r.tips()][0] == "e"
    tips = {t.name: t for t in r.tips()}
    assert T.mrca([tips["g"], tips["f"]]) is tips["f"].parent
    assert T.mrca([tips["a"], tips["d"]]) is r
    L = T.layout(r)
    assert L.ntips == 8 and L.has_lengths and L.width == 4
    assert [L.y[id(r.tips()[i])] for i in range(8)] == list(range(8))
    assert L.x[id(tips["h"])] == 4 and L.x[id(tips["d"])] == 1
    C = T.layout(r, cladogram=True)
    assert not C.has_lengths and {C.x[id(t)] for t in r.tips()} == {C.width}
    assert not T.layout(T.parse_newick("((a,b),c);")).has_lengths      # no lengths: drawn as a cladogram


def test_nexus_and_names_file(tmp_path):
    nex = tmp_path / "t.nex"
    nex.write_text("#NEXUS\nbegin taxa;\n dimensions ntax=3;\nend;\nbegin trees;\n translate\n"
                   "  1 'Homo sapiens',\n  2 Pan_troglodytes,\n  3 'Bob''s ape'\n ;\n"
                   " tree tree_1 = [&R] ((1:0.1,2:0.2)95:0.05,3:0.3);\nend;\n")
    r = T.read_tree(str(nex))
    assert [t.name for t in r.tips()] == ["Homo sapiens", "Pan_troglodytes", "Bob's ape"]
    run = tmp_path / "run"
    run.mkdir()
    (run / "demo.nj.nwk").write_text("(a:1,b:1,c:1);")
    (run / "demo.treefile").write_text("(a:1,b:1,c:1);")
    (run / "demo.names.tsv").write_text("label\tname\na\tA full name\nb\tB, with comma\n")
    for f in ("demo.nj.nwk", "demo.treefile"):
        p = T.names_file_for(str(run / f))
        assert p == str(run / "demo.names.tsv")
    assert T.read_names(p) == {"a": "A full name", "b": "B, with comma"}
    assert T.names_file_for(str(tmp_path / "t.nex")) is None


# ---------------------------------------------------------------- viewer window
def _app():
    import sys
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def test_tree_window(tmp_path):
    app = _app()
    from PySide6.QtCore import Qt, QPointF
    from PySide6.QtTest import QTest
    from neoedit.ui.tree_view import TreeWindow, nice_length
    from neoedit.ui.tree_view import spaced
    assert [nice_length(x) for x in (0.013, 0.04, 0.3, 7)] == [0.01, 0.05, 0.2, 5]
    assert spaced("NC_012920.1_Homo_sapiens_mitochondrion") == "NC_012920.1 Homo sapiens mitochondrion"
    assert spaced("Pimephales_promelas_MT455673.1") == "Pimephales promelas MT455673.1"
    assert spaced("WP_000001.1_x") == "WP_000001.1 x"
    assert spaced("Pan_troglodytes") == "Pan troglodytes" and spaced("a__b") == "a  b"
    run = tmp_path / "run"
    run.mkdir()
    (run / "demo.nj.nwk").write_text(
        "(Homo_sapiens:0.1,(Pan_troglodytes:0.05,Pan_paniscus:0.06)97:0.02,(Mus_musculus:0.3,Rattus_rattus:0.28)100:0.4);\n")
    (run / "demo.names.tsv").write_text("label\tname\nHomo_sapiens\tNC_012920.1 Homo sapiens mitochondrion\n")
    picked = []
    w = TreeWindow(select_names=lambda labels, names: picked.append((labels, names)) or len(labels))
    w.load(str(run / "demo.nj.nwk"))
    w.resize(800, 500); w.show(); QTest.qWait(30)
    c = w.canvas
    assert "demo.nj.nwk" in w.windowTitle() and w.a_full.isEnabled() and w.a_full.isChecked()
    assert c.tip_text(c.root.tips()[0]) == "NC_012920.1 Homo sapiens mitochondrion"
    assert c.tip_text(c.root.tips()[1]) == "Pan troglodytes"        # no full name: label with spaces
    w.a_full.trigger()
    assert c.tip_text(c.root.tips()[0]) == "Homo sapiens"
    img = c.grab()
    assert img.width() > 300 and c.height() > 100
    # click the branch above the (Mus, Rattus) clade: it is selected
    mouse = next(n for n in c.root.walk() if n.label == "100")
    rect = next(r for r, n in c._hits if n is mouse)
    QTest.mouseClick(c, Qt.LeftButton, Qt.NoModifier, rect.center().toPoint())
    assert c.selected is mouse and w.a_root.isEnabled() and w.a_select.isEnabled()
    assert "Clade of 2 sequences, support 100" in w.statusBar().currentMessage()
    w.a_select.trigger()
    assert picked[-1][0] == ["Mus_musculus", "Rattus_rattus"] and "Homo_sapiens" in picked[-1][1]
    w.a_root.trigger()
    assert c.selected is None and len(c.root.children) == 2
    assert sorted(t.name for t in c.root.children[0].tips()) == ["Mus_musculus", "Rattus_rattus"]
    before = T.to_newick(c.root)
    w.a_ladder.trigger(); w.a_mid.trigger()
    w.a_reset.trigger()
    assert T.to_newick(c.root) == "(Homo_sapiens:0.1,(Pan_troglodytes:0.05,Pan_paniscus:0.06)97:0.02," \
                                  "(Mus_musculus:0.3,Rattus_rattus:0.28)100:0.4);"
    assert before != T.to_newick(c.root)
    # rotate a clade, cladogram, zoom, exports
    pan = next(n for n in c.root.walk() if n.label == "97")
    c.select(pan); w.a_rotate.trigger()
    assert [t.name for t in pan.tips()] == ["Pan_paniscus", "Pan_troglodytes"]
    w.a_clado.trigger(); assert c._layout and not c._layout.has_lengths
    h0, w0 = c.sizeHint().height(), c.sizeHint().width()
    w.a_zin.trigger(); app.processEvents()
    assert c.sizeHint().height() > h0 and c.sizeHint().width() > w0
    w.a_fit.trigger()
    c.export_image(str(tmp_path / "t.svg"))
    c.export_image(str(tmp_path / "t.png"))
    assert "Pan paniscus" in (tmp_path / "t.svg").read_text()
    assert (tmp_path / "t.png").stat().st_size > 1000
    # clicking empty space clears the selection
    QTest.mouseClick(c, Qt.LeftButton, Qt.NoModifier, QPointF(c.width() - 3, 3).toPoint())
    assert c.selected is None and not w.a_root.isEnabled()
    w.close()


def test_main_window_trees(tmp_path):
    _app()
    from PySide6.QtTest import QTest
    from neoedit.ui.main_window import MainWindow
    from neoedit.model.alignment import AlignmentModel, SequenceRow
    from neoedit.ui.dialogs.tree_dialogs import NJDialog
    from neoedit.ui.tree_view import TreeWindow
    w = MainWindow()
    seqs = ["ACGTACGTACGTAAGGTTCC", "ACGTACGTACGTAAGGTTCA", "ACGAACGTTCGTAAGCTTCC", "TCGAACCTTCGTTAGCTACC"]
    names = ["NC_1.1 Homo sapiens, mito", "NC_2.1 Pan troglodytes", "AB3.1 Mus musculus", "OQ4.1 Rattus"]
    w._set_model(AlignmentModel([SequenceRow(n, s, id=n.split()[0]) for n, s in zip(names, seqs)]))
    # tree labels, underscores-for-spaces, accessions and names.tsv all find their rows
    names_map = {"x1": "AB3.1 Mus musculus"}
    assert w.select_tree_names(["NC_1.1_Homo_sapiens_mito", "NC_2.1 Pan troglodytes", "x1", "OQ4.1", "nope"],
                               names_map) == 4
    assert w.view.sel_rows == {0, 1, 2, 3}
    assert w.select_tree_names(["NC_2.1_Pan_troglodytes"]) == 1 and w.view.sel_rows == {1}
    # a finished NJ run opens the viewer
    d = NJDialog(w, w.settings, w.model, [], None)
    d.out_dir.setText(str(tmp_path))
    d.reps.setValue(20)
    d.start()
    for _ in range(200):
        if not d.running():
            break
        QTest.qWait(25)
    QTest.qWait(50)
    viewers = [c for c in w._children if isinstance(c, TreeWindow)]
    assert d.tree is not None and len(viewers) == 1 and viewers[0].isVisible()
    tv = viewers[0]
    assert sorted(t.name for t in tv.canvas.root.tips()) == \
        ["AB3.1_Mus_musculus", "NC_1.1_Homo_sapiens_mito", "NC_2.1_Pan_troglodytes", "OQ4.1_Rattus"]
    assert tv.a_full.isEnabled()
    tv.canvas.select(tv.canvas.root)
    tv.a_select.trigger()
    assert w.view.sel_rows == {0, 1, 2, 3}
    # tree files open in the viewer, also through the normal open path
    nwk = tmp_path / "x.treefile"
    nwk.write_text("(a:1,b:1,c:1);")
    w.open_path(str(nwk))
    assert len([c for c in w._children if isinstance(c, TreeWindow)]) == 2
    assert w.model.nrows == 4                                     # the alignment is untouched
    bad = tmp_path / "bad.nwk"
    bad.write_text("not a tree")
    from PySide6.QtWidgets import QMessageBox
    orig = QMessageBox.warning
    QMessageBox.warning = staticmethod(lambda *a, **k: None)
    try:
        assert w.open_tree_path(str(bad)) is None
    finally:
        QMessageBox.warning = orig
    for c in list(w._children):
        c.close()
    w.model.dirty = False
    w.close()
