"""BioEdit compatibility: binary .bio project files and BioEdit's GenBank dialect."""
import os
import struct
import pytest

from neoedit.model import bioedit_format as B
from neoedit.model import io as mio
from neoedit.model import AlignmentModel, SequenceRow

REAL = "/mnt/c/Users/User/Desktop/Sander_mitogenomesbio.bio"


def _rows():
    return [SequenceRow("seq one", "ACGTACGTAC"), SequenceRow("seq two", "ACGT--GTAC"),
            SequenceRow("third", "ACGTACGTAA")]


def test_write_read_roundtrip(tmp_path):
    p = tmp_path / "x.bio"
    rows = _rows()
    B.write_bio(str(p), rows, "dna")
    assert B.is_bio_file(str(p))
    back, info = B.read_bio(str(p))
    assert info["version"] == 7 and info["n_sequences"] == 3
    assert [r.name for r in back] == [r.name for r in rows]
    assert [r.seq for r in back] == [r.seq for r in rows]


def test_written_layout_matches_bioedit(tmp_path):
    """Header offsets and record field pattern must follow the documented layout."""
    p = tmp_path / "x.bio"
    B.write_bio(str(p), _rows(), "dna")
    buf = p.read_bytes()
    assert buf[:22] == B.MAGIC_PREFIX
    assert struct.unpack_from("<i", buf, 0x18)[0] == 3
    assert struct.unpack_from("<i", buf, 0x1C)[0] == -1
    assert struct.unpack_from("<i", buf, 0x20)[0] == -1
    offs = struct.unpack_from("<3i", buf, 0xC8)
    assert offs[0] == 0xC8 + 12
    # first record: title, sequence, type, then the tail fields
    pos = offs[0]
    lens = []
    for _ in range(15):
        (n,) = struct.unpack_from("<i", buf, pos)
        lens.append(n); pos += 4 + n
    assert lens[:3] == [len("seq one"), 10, 3]
    assert lens[3:] == [1] + [0] * 11


def test_io_layer_detects_and_saves(tmp_path):
    p = tmp_path / "aln.bio"
    m = AlignmentModel(_rows())
    mio.save(m, str(p), "bio")
    assert mio.guess_format(str(p)) == "bio"
    m2 = mio.load(str(p))
    assert m2.nrows == 3 and m2.format == "bio"
    assert [r.seq for r in m2.rows] == [r.seq for r in _rows()]


def test_genbank_bioedit_dialect(tmp_path):
    m = AlignmentModel([SequenceRow("My long title with spaces", "ACGTACGTAAACGTACGTAA",
                                    description="a description")])
    text = mio.dumps(m, "genbank-bioedit")
    assert "LOCUS" in text and "TITLE       My long title with spaces" in text
    assert "ORIGIN" in text and text.rstrip().endswith("//")
    # residues in blocks of five, fifty per line
    body = [l for l in text.splitlines() if l[:6].strip().isdigit()]
    assert body and body[0].split()[1] == "ACGTA"
    p = tmp_path / "b.gb"
    p.write_text(text)
    assert B.genbank_titles(str(p)) == ["My long title with spaces"]
    m2 = mio.load(str(p), "genbank-bioedit")
    assert m2.rows[0].name == "My long title with spaces"


def test_genbank_titles_ignores_reference_title(tmp_path):
    gb = ("LOCUS       X                 10 bp\n"
          "DEFINITION  something\n"
          "TITLE       Real title\n"
          "REFERENCE   1  (bases 1 to 10)\n"
          "  AUTHORS   Someone\n"
          "  TITLE     A paper title that is not the sequence title\n"
          "ORIGIN\n"
          "        1 acgtacgtaa\n"
          "//\n")
    p = tmp_path / "r.gb"; p.write_text(gb)
    assert B.genbank_titles(str(p)) == ["Real title"]


@pytest.mark.skipif(not os.path.exists(REAL), reason="no real BioEdit .bio file available")
def test_real_bioedit_file_roundtrip(tmp_path):
    rows, info = B.read_bio(REAL)
    assert info["n_sequences"] == len(rows) == 71
    assert all(len(r.seq) == 16944 for r in rows)
    assert rows[0].name.startswith("OR552089.1") and rows[0].description == "DNA"
    # our writer reproduces BioEdit's own bytes exactly
    out = tmp_path / "rt.bio"
    B.write_bio(str(out), rows, "dna")
    assert out.read_bytes() == open(REAL, "rb").read()
