from neoedit.analysis import restriction as R


def test_pool_filters():
    b = R.enzyme_pool(names=["EcoRI", "BamHI"])
    assert {str(e) for e in b} == {"EcoRI", "BamHI"}
    b6 = R.enzyme_pool(commercial_only=True, min_site=6, max_site=6)
    assert all(6 <= e.size <= 6 for e in b6) and len(b6) > 50
    bl = R.enzyme_pool(commercial_only=True, min_site=6, max_site=6, blunt=True)
    assert all(e.is_blunt() for e in bl)
    assert "New England Biolabs" in R.all_suppliers() or R.all_suppliers()


def test_search_sequence():
    seq = "AAAGAATTCTTTGGATCCAAA"          # EcoRI GAATTC at 4, BamHI GGATCC at 13
    hits = R.search_sequence(seq, R.enzyme_pool(names=["EcoRI", "BamHI", "HindIII"]))
    by = {h.enzyme: h for h in hits}
    assert by["EcoRI"].positions == [5] and by["EcoRI"].overhang == "5'"
    assert by["BamHI"].n_cuts == 1
    assert "HindIII" not in by                      # no site -> not reported
    hits = R.search_sequence(seq, R.enzyme_pool(names=["EcoRI"]), min_cuts=2)
    assert hits == []


def test_search_alignment_diagnostic():
    a = "AAAGAATTCTTTGGATCCAAA"
    b = a.replace("GAATTC", "GAATTA")               # loses EcoRI
    rows = [a, b]
    summ = {s.enzyme: s for s in R.search_alignment(rows, R.enzyme_pool(names=["EcoRI", "BamHI"]))}
    assert summ["BamHI"].is_universal([0, 1]) and not summ["BamHI"].is_diagnostic([0, 1])
    assert summ["EcoRI"].is_diagnostic([0, 1]) and summ["EcoRI"].rows_cut() == [0]


def test_features_and_digest():
    seq = "AAAGAATTCTTTGGATCCAAA"
    gapped = "AAA-GAATTC-TTTGGATCCAAA"
    hits = R.search_sequence(gapped, R.enzyme_pool(names=["EcoRI"]), row=0)
    feats = R.hits_to_features(hits, {0: gapped})
    assert len(feats) == 1
    f = feats[0]
    assert "".join(c for c in gapped[f.start:f.end] if c != "-") == "GAATTC"
    frags = R.digest_fragments(seq, ["EcoRI", "BamHI"])
    assert sum(l for _s, _e, l in frags) == len(seq) and len(frags) == 3
