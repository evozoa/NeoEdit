from neoedit.model import AlignmentModel, SequenceRow


def mk():
    return AlignmentModel([SequenceRow("a", "ACGT"), SequenceRow("b", "AC-T")])


def test_insert_delete_gaps_undo():
    m = mk()
    m.insert_gaps(0, 2, 2)
    assert m.rows[0].seq == "AC--GT"
    assert m.delete_gaps(0, 2, 5) == 2  # stops at non-gap
    assert m.rows[0].seq == "ACGT"
    m.undo(); assert m.rows[0].seq == "AC--GT"
    m.undo(); assert m.rows[0].seq == "ACGT"
    m.redo(); assert m.rows[0].seq == "AC--GT"
    assert m.delete_gaps(0, 0) == 0


def test_columns():
    m = mk()
    m.insert_gap_columns(1, 1)
    assert [r.seq for r in m.rows] == ["A-CGT", "A-C-T"]
    assert m.delete_gap_columns(1, 3) == 1
    assert [r.seq for r in m.rows] == ["ACGT", "AC-T"]
    m.insert_gap_columns(4, 2)
    m.remove_gap_only_columns()
    assert [r.seq for r in m.rows] == ["ACGT", "AC-T"]


def test_overwrite_insert_pad():
    m = mk()
    m.overwrite(1, 6, "GG")
    assert m.rows[1].seq == "AC-T--GG"
    m.insert_text(0, 1, "NN")
    assert m.rows[0].seq == "ANNCGT"


def test_block_shift():
    m = AlignmentModel([SequenceRow("a", "AC--GT")])
    assert m.block_shift([0], 0, 2, 1)
    assert m.rows[0].seq == "-AC-GT"
    assert not m.block_shift([0], 1, 3, 2)  # would overwrite G
    assert m.block_shift([0], 4, 6, -1)
    assert m.rows[0].seq == "-ACGT-"


def test_rows_and_transforms():
    m = mk()
    m.reverse_complement([0]); assert m.rows[0].seq == "ACGT"
    m.reverse_complement([1]); assert m.rows[1].seq == "A-GT"
    m.remove_rows([0]); assert m.nrows == 1
    m.undo(); assert m.nrows == 2 and m.rows[0].name == "a"
    m.add_row(SequenceRow("c", "AAAA")); assert m.nrows == 3
    m.undo(); assert m.nrows == 2
    m.move_rows([1], -1); assert m.rows[0].name == "b"
    m.undo(); assert m.rows[0].name == "a"


def test_find():
    m = AlignmentModel([SequenceRow("a", "AC-GTT"), SequenceRow("b", "TTTACG")])
    assert m.find("ACG") == (0, 0, 4)
    assert m.find("ACG", 0, 1) == (1, 3, 6)
    assert m.find("ACG", 1, 4) == (0, 0, 4)  # wraps
    assert m.find("zzz") is None


def test_detect_type():
    assert mk().seq_type == "dna"
    assert AlignmentModel([SequenceRow("p", "MKVLAAG")]).seq_type == "protein"


def test_bioedit_color_tables():
    from neoedit.model import colors as C
    assert C.SCHEMES_NT["BioEdit"]["A"] == "#008000"
    assert C.SCHEMES_NT["BioEdit"]["C"] == "#0000ff"
    assert C.SCHEMES_NT["BioEdit"]["G"] == "#000000"
    assert C.SCHEMES_NT["BioEdit"]["T"] == "#ff0000"
    assert C.SCHEMES_AA["BioEdit"]["D"] == "#ff0000" and C.SCHEMES_AA["BioEdit"]["R"] == "#0000ff"
    assert list(C.SCHEMES_NT)[0] == "BioEdit" and list(C.SCHEMES_AA)[0] == "BioEdit"
    assert C._bgr_to_hex(16711680) == "#0000ff"
