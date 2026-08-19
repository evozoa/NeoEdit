import os
import pytest
from neoedit.model import AlignmentModel, SequenceRow
from neoedit.model import io as mio


@pytest.mark.parametrize("fmt", ["fasta", "clustal", "phylip-relaxed", "nexus", "stockholm", "bio"])
def test_roundtrip(tmp_path, fmt):
    m = AlignmentModel([SequenceRow("seq1", "ACGT-ACGT"), SequenceRow("seq2", "ACGTTACGT")])
    p = tmp_path / f"x.{fmt}"
    mio.save(m, str(p), fmt)
    m2 = mio.load(str(p), fmt)
    assert [r.seq.upper() for r in m2.rows] == [r.seq for r in m.rows]
    assert [r.name for r in m2.rows] == ["seq1", "seq2"]


def test_guess(tmp_path):
    p = tmp_path / "a.txt"
    p.write_text(">x\nACGT\n")
    assert mio.guess_format(str(p)) == "fasta"
    p.write_text("CLUSTAL W\n\nx ACGT\n")
    assert mio.guess_format(str(p)) == "clustal"
    assert mio.guess_format("foo.nex") == "nexus"


def test_genbank(tmp_path):
    m = AlignmentModel([SequenceRow("seq1", "ACGTACGTAA")])
    p = tmp_path / "x.gb"
    mio.save(m, str(p), "genbank")
    m2 = mio.load(str(p))
    assert m2.rows[0].seq == "ACGTACGTAA"
