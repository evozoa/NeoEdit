"""Load / save alignments via Biopython, plus a best-effort BioEdit .bio reader."""
from __future__ import annotations

import io
import os
import re
from typing import Optional

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .alignment import AlignmentModel, SequenceRow, GAP_CHARS

# (label, biopython format name, extensions, is_alignment_format)
FORMATS = [
    ("FASTA", "fasta", [".fasta", ".fas", ".fa", ".fna", ".faa", ".seq", ".txt"], False),
    ("Clustal", "clustal", [".aln", ".clustal", ".clw"], True),
    ("PHYLIP (sequential)", "phylip-sequential", [".phy", ".phylip"], True),
    ("PHYLIP (interleaved)", "phylip", [".phy"], True),
    ("PHYLIP (relaxed)", "phylip-relaxed", [".phy"], True),
    ("NEXUS", "nexus", [".nex", ".nexus", ".nxs"], True),
    ("Stockholm", "stockholm", [".sto", ".stk"], True),
    ("GenBank", "genbank", [".gb", ".gbk", ".genbank"], False),
    ("EMBL", "embl", [".embl"], False),
    ("MSF", "msf", [".msf"], True),
    ("BioEdit (.bio)", "bio", [".bio"], True),
]

EXT_TO_FORMAT = {}
for _label, _fmt, _exts, _ in FORMATS:
    for e in _exts:
        EXT_TO_FORMAT.setdefault(e, _fmt)


def guess_format(path: str, text_head: str | None = None) -> str:
    ext = os.path.splitext(path)[1].lower()
    if ext in EXT_TO_FORMAT and ext != ".txt":
        return EXT_TO_FORMAT[ext]
    head = text_head
    if head is None:
        try:
            with open(path, "r", errors="replace") as fh:
                head = fh.read(4096)
        except OSError:
            head = ""
    h = head.lstrip()
    if h.startswith(">"):
        return "fasta"
    if h.upper().startswith("CLUSTAL"):
        return "clustal"
    if h.upper().startswith("#NEXUS"):
        return "nexus"
    if h.startswith("# STOCKHOLM"):
        return "stockholm"
    if h.startswith("LOCUS"):
        return "genbank"
    if h.startswith("ID "):
        return "embl"
    if re.match(r"^\s*\d+\s+\d+\s*$", h.splitlines()[0] if h else ""):
        return "phylip-relaxed"
    if "BioEdit" in head or "Sequence alignment" in head:
        return "bio"
    return "fasta"


def _rows_from_records(records) -> list[SequenceRow]:
    rows = []
    for rec in records:
        name = rec.id if rec.id and rec.id != "<unknown id>" else (rec.name or "seq")
        desc = rec.description or ""
        if desc.startswith(name):
            desc = desc[len(name):].strip()
        rows.append(SequenceRow(name=name, seq=str(rec.seq), description=desc, id=rec.id or ""))
    return rows


def load(path: str, fmt: Optional[str] = None) -> AlignmentModel:
    fmt = fmt or guess_format(path)
    if fmt == "bio":
        rows = _read_bio(path)
    else:
        with open(path, "r", errors="replace") as fh:
            try:
                records = list(SeqIO.parse(fh, fmt))
            except Exception:
                fh.seek(0)
                aln = AlignIO.read(fh, fmt)
                records = list(aln)
        rows = _rows_from_records(records)
    model = AlignmentModel(rows)
    model.path = path
    model.format = fmt
    model.dirty = False
    # hoist GenBank/EMBL features into the model
    if fmt in ("genbank", "embl"):
        from .alignment import Feature
        with open(path, "r", errors="replace") as fh:
            for ri, rec in enumerate(SeqIO.parse(fh, fmt)):
                for f in rec.features:
                    if f.type == "source":
                        continue
                    label = f.qualifiers.get("gene", f.qualifiers.get("product", f.qualifiers.get("label", [f.type])))[0]
                    model.features.append(Feature(ri, int(f.location.start), int(f.location.end),
                                                  f.location.strand or 1, f.type, str(label)))
    return model


def loads(text: str, fmt: str) -> AlignmentModel:
    fh = io.StringIO(text)
    records = list(SeqIO.parse(fh, fmt))
    m = AlignmentModel(_rows_from_records(records))
    m.format = fmt
    m.dirty = False
    return m


def to_records(model: AlignmentModel, rows=None) -> list[SeqRecord]:
    idx = list(rows) if rows is not None else range(model.nrows)
    recs = []
    for i in idx:
        r = model.rows[i]
        recs.append(SeqRecord(Seq(r.seq), id=r.name, name=r.name[:16], description=r.description))
    return recs


def dumps(model: AlignmentModel, fmt: str, rows=None) -> str:
    recs = to_records(model, rows)
    _label, _, _, is_aln = next((f for f in FORMATS if f[1] == fmt), (fmt, fmt, [], False))
    out = io.StringIO()
    if fmt == "bio":
        return _write_bio(model)
    if is_aln:
        w = max((len(r.seq) for r in recs), default=0)
        for r in recs:
            if len(r.seq) < w:
                r.seq = Seq(str(r.seq) + "-" * (w - len(r.seq)))
        aln = MultipleSeqAlignment(recs)
        if fmt == "nexus":
            mt = "dna" if model.seq_type in ("dna", "rna") else "protein"
            for r in aln:
                r.annotations["molecule_type"] = mt.upper()
        AlignIO.write(aln, out, fmt)
    else:
        if fmt in ("genbank", "embl"):
            mt = "DNA" if model.seq_type in ("dna", "rna") else "protein"
            for r in recs:
                r.annotations["molecule_type"] = mt
                r.seq = Seq("".join(c for c in str(r.seq) if c not in GAP_CHARS))
        SeqIO.write(recs, out, fmt)
    return out.getvalue()


def save(model: AlignmentModel, path: str, fmt: Optional[str] = None):
    fmt = fmt or model.format or guess_format(path)
    text = dumps(model, fmt)
    with open(path, "w") as fh:
        fh.write(text)
    model.path = path
    model.format = fmt
    model.dirty = False


# ------------------------------------------------------------- BioEdit .bio
# BioEdit's native format is a text file with a header and then a "sequence
# data" block. The exact layout varies across versions; this reader looks for
# the block of "name  sequence" lines (possibly with a number in between) and
# falls back to treating anything after a line that ends with "data" as rows.

def _read_bio(path: str) -> list[SequenceRow]:
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    rows: dict[str, SequenceRow] = {}
    order: list[str] = []
    in_data = False
    for ln in lines:
        if not in_data:
            if re.search(r"sequence\s*data|^data$|^seqs?\s*$", ln, re.I):
                in_data = True
            continue
        if not ln.strip():
            continue
        m = re.match(r"^(\S.*?)\s{2,}(?:\d+\s+)?([A-Za-z\-\.\~\*\?]+)\s*$", ln)
        if not m:
            m = re.match(r"^(\S+)\s+([A-Za-z\-\.\~\*\?]+)\s*$", ln)
        if not m:
            continue
        name, seq = m.group(1).strip(), m.group(2)
        if name not in rows:
            rows[name] = SequenceRow(name, "")
            order.append(name)
        rows[name].seq += seq
    if not rows:
        # last resort: maybe it's actually FASTA
        with open(path, "r", errors="replace") as fh:
            try:
                return _rows_from_records(SeqIO.parse(fh, "fasta"))
            except Exception:
                return []
    return [rows[n] for n in order]


def _write_bio(model: AlignmentModel) -> str:
    out = ["NeoEdit alignment file", "", "sequence data"]
    w = max((len(r.name) for r in model.rows), default=4)
    for r in model.rows:
        out.append(f"{r.name.ljust(w)}  {r.seq}")
    return "\n".join(out) + "\n"
