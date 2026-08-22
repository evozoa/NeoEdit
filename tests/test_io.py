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


def test_names_are_full_deflines_and_roundtrip(tmp_path):
    """Rows are named as the source shows them (id + description); the accession stays the key;
    writing FASTA reproduces the header without duplicating the description."""
    import os
    from neoedit.model import io as mio
    fa = tmp_path / "a.fasta"
    fa.write_text(">NC_012920.1 Homo sapiens mitochondrion, complete genome\nACGTACGT\n>seq2\nACGT\n")
    m = mio.load(str(fa))
    assert m.rows[0].name == "NC_012920.1 Homo sapiens mitochondrion, complete genome"
    assert m.rows[0].accession == "NC_012920.1" and m.seqid(0) == "NC_012920.1"
    assert m.rows[1].name == "seq2" and m.rows[1].accession == "seq2"
    out = mio.dumps(m, "fasta")
    assert out.startswith(">NC_012920.1 Homo sapiens mitochondrion, complete genome\n")
    assert out.count("Homo sapiens") == 1
    m2 = mio.loads(out, "fasta")
    assert m2.rows[0].name == m.rows[0].name
    gb = os.path.join(os.path.dirname(__file__), "..", "examples", "mito", "NC_012920_MDP.gb")
    if os.path.exists(gb):
        g = mio.load(gb)
        assert g.rows[0].name == "NC_012920.1 Homo sapiens mitochondrion, complete genome"
        assert g.seqid(0) == "NC_012920.1"
