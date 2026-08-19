from neoedit.analysis import translate as T
from neoedit.analysis.orf_finder import find_orfs, orfs_to_gff, gap_map
from neoedit.analysis.primers import design_primers, primer_stats, primer_mismatches


def test_translate():
    assert T.translate_gapped("ATG-GCC-TAA") == "MA*"
    assert T.translate_gapped("ATGGCCTAA", to_stop=True) == "MA"
    assert T.translate_aligned("ATG---GCC") == "M-A"
    assert T.translate_aligned("ATGG--GCC") == "MXA"
    sf = T.six_frame("ATGGCCTAA")
    assert sf["+1"] == "MA*" and len(sf) == 6
    # mito table 2: AGA is stop, ATA is Met
    assert T.translate_gapped("ATAAGA", table=2) == "M*"


def test_consensus_identity():
    rows = ["ACGT", "ACGA", "AC-T"]
    assert T.consensus(rows) == "ACGT"
    mat = T.identity_matrix(rows)
    assert mat[0, 0] == 1.0 and abs(mat[0, 1] - 0.75) < 1e-9
    ent = T.column_entropy(rows)
    assert ent[0] == 0 and ent[3] > 0
    assert T.gc_content("GGCC-AT") == 4 / 6


def test_orfs_basic():
    # ATG AAA CCC GGG TAA  -> 4 aa ORF on + strand
    seq = "TTT" + "ATGAAACCCGGGTAA" + "TTT"
    orfs = find_orfs(seq, table=1, min_aa=3, both_strands=False, allow_partial=False)
    assert len(orfs) == 1
    o = orfs[0]
    assert (o.start, o.end, o.strand, o.aa) == (3, 18, 1, "MKPG*")
    assert o.start_codon == "ATG" and o.stop_codon == "TAA"


def test_orfs_minus_and_partial():
    from Bio.Seq import Seq
    plus = "ATGAAACCCGGGTAA"
    seq = "CC" + str(Seq(plus).reverse_complement()) + "CC"
    orfs = find_orfs(seq, table=1, min_aa=3, allow_partial=False)
    neg = [o for o in orfs if o.strand < 0]
    assert len(neg) == 1 and neg[0].start == 2 and neg[0].end == 17 and neg[0].aa == "MKPG*"
    # 3' partial: no stop
    orfs = find_orfs("ATGAAACCCGGGAAA", min_aa=3, both_strands=False, allow_partial=True)
    assert any(o.partial3 and not o.partial5 for o in orfs)
    # mito: ATA start, AGA stop (table 2); T-only incomplete stop -> partial3
    orfs = find_orfs("ATAAAACCCGGGAGA", table=2, min_aa=3, both_strands=False, allow_partial=False, start_mode="table")
    assert orfs and orfs[0].aa == "MKPG*" and orfs[0].start_codon == "ATA"
    orfs = find_orfs("ATAAAACCCGGGAGA", table=2, min_aa=3, both_strands=False, start_mode="atg")
    assert not [o for o in orfs if not o.partial5]


def test_orf_nested_and_exports():
    seq = "ATGATGAAACCCGGGTAA"
    o1 = find_orfs(seq, min_aa=3, both_strands=False, allow_partial=False)
    o2 = find_orfs(seq, min_aa=3, both_strands=False, allow_partial=False, nested=True)
    assert len(o1) == 1 and len(o2) == 2
    gff = orfs_to_gff(o1, "s")
    assert "CDS\t1\t18" in gff
    assert gap_map("A-C") == [0, 2]


def test_primers():
    import random
    random.seed(1)
    tmpl = "".join(random.choice("ACGT") for _ in range(400))
    pairs = design_primers(tmpl, target=(180, 40), product_range=((100, 300),), num_return=3)
    assert pairs and pairs[0].product_size >= 100
    p = pairs[0]
    assert tmpl[p.left.start:p.left.start + p.left.length] == p.left.seq
    assert primer_mismatches(p.left.seq, [tmpl], 1, p.left.start) == [0]
    assert primer_mismatches(p.right.seq, [tmpl], -1, p.right.start) == [0]
    st = primer_stats("ACGTACGTACGTACGTACGT")
    assert st["tm"] > 0
