"""BioEdit compatibility: the binary .bio "BioEdit Project File" and BioEdit's GenBank dialect.

The .bio layout follows the description in the BioEdit help file (v7.x, "BioEdit
Project File Format"):

    0x00  magic  "**BioEdit Project File**"   (version 1)
                 "**BioEdit Project File02"   (version 2, current)
    0x18  int32  number of sequences
    0x1C  int32  index of the mask sequence      (-1 / 0xFFFFFFFF if none)
    0x20  int32  index of the numbering mask     (   "     "        "   )
    0xC8  int32[n]  file offset of each sequence structure

Each sequence structure is a run of length-prefixed fields, every field being an
int32 byte count followed by that many bytes: title, sequence, sequence type, then
the GenBank-derived fields, and (v5+) annotations, grouping, consensus, lock state
and positional flags. Values are little-endian (BioEdit is a Windows/x86 program).

The exact order and meaning of the trailing fields is not published, so the reader
takes the documented leading fields (title, sequence, type) and ignores the rest;
unreadable files fall back to a text scan. Round-tripping through NeoEdit is
covered by tests, but the writer has not been verified against BioEdit itself.
"""
from __future__ import annotations

import re
import struct
from dataclasses import dataclass

from .alignment import SequenceRow

MAGIC_PREFIX = b"**BioEdit Project File"
MAGIC_V1 = MAGIC_PREFIX + b"**"          # first version
MAGIC_V2 = MAGIC_PREFIX + b"02"          # version 2 (per the manual)
MAGIC_V7 = MAGIC_PREFIX + b"07"          # written by BioEdit 7.x
DEFAULT_MAGIC = MAGIC_V7

# Byte layout of one sequence record as written by BioEdit 7.x, verified against a
# real .bio file: title, sequence and type are length-prefixed fields, followed by a
# one-byte field, eleven empty fields, and a fixed 31-byte trailer holding the
# lock/group/consensus/flag state (all "unset" here).
RECORD_TAIL_FIELDS = [b"\x00"] + [b""] * 11
RECORD_TRAILER = (b"\x00\x00\x00" + b"\xff\xff\xff\xff" + b"\x00" * 8
                  + b"\xff\xff\xff\xff" + b"\x00" * 8 + b"\xff\xff\xff\xff")
FILE_TRAILER = b"\x00" * 20      # BioEdit 7.x pads the end of the file
HEADER_SIZE = 0xC8
OFF_NSEQ = 0x18
OFF_MASK = 0x1C
OFF_NUMMASK = 0x20

# sequence type codes used in the "type" field
TYPE_DNA, TYPE_RNA, TYPE_PROTEIN, TYPE_UNKNOWN = "DNA", "RNA", "Protein", "Unknown"
_RESIDUE_RE = re.compile(rb"^[A-Za-z\-\.\~\*\?\s]+$")


def is_bio_file(path: str) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(len(MAGIC_PREFIX)) == MAGIC_PREFIX
    except OSError:
        return False


def _read_field(buf: bytes, pos: int) -> tuple[bytes, int]:
    """Read one length-prefixed field; returns (data, next_pos)."""
    if pos + 4 > len(buf):
        raise ValueError("truncated field header")
    (n,) = struct.unpack_from("<i", buf, pos)
    pos += 4
    if n < 0 or pos + n > len(buf):
        raise ValueError(f"bad field length {n} at {pos - 4}")
    return buf[pos:pos + n], pos + n


def _looks_like_sequence(b: bytes) -> bool:
    if not b:
        return False
    sample = b[:200]
    return bool(_RESIDUE_RE.match(sample))


def read_bio(path: str) -> tuple[list[SequenceRow], dict]:
    """Parse a BioEdit Project file. Returns (rows, info)."""
    with open(path, "rb") as fh:
        buf = fh.read()
    if not buf.startswith(MAGIC_PREFIX):
        raise ValueError("not a BioEdit Project file")
    tag = buf[len(MAGIC_PREFIX):len(MAGIC_PREFIX) + 2]
    version = int(tag) if tag.isdigit() else 1
    (nseq,) = struct.unpack_from("<i", buf, OFF_NSEQ)
    (mask_idx,) = struct.unpack_from("<i", buf, OFF_MASK)
    (num_idx,) = struct.unpack_from("<i", buf, OFF_NUMMASK)
    info = {"version": version, "n_sequences": nseq, "mask_index": mask_idx, "numbering_index": num_idx}
    if not (0 < nseq <= 1_000_000):
        raise ValueError(f"implausible sequence count {nseq}")
    offsets = list(struct.unpack_from(f"<{nseq}i", buf, HEADER_SIZE))
    rows: list[SequenceRow] = []
    for i, off in enumerate(offsets):
        if not (HEADER_SIZE <= off < len(buf)):
            continue
        try:
            title, pos = _read_field(buf, off)
            seq, pos = _read_field(buf, pos)
            stype = b""
            try:
                stype, pos = _read_field(buf, pos)
            except ValueError:
                pass
        except ValueError:
            continue
        name = title.decode("latin-1", "replace").strip() or f"sequence_{i + 1}"
        s = seq.decode("latin-1", "replace")
        if not _looks_like_sequence(seq) and _looks_like_sequence(title):
            name, s = f"sequence_{i + 1}", title.decode("latin-1", "replace")
        s = re.sub(r"\s", "", s)
        rows.append(SequenceRow(name=name, seq=s, description=stype.decode("latin-1", "replace").strip()))
    if not rows:
        raise ValueError("no sequences could be read")
    info["read"] = len(rows)
    return rows, info


def _seq_type_code(seq_type: str) -> str:
    return {"dna": TYPE_DNA, "rna": TYPE_RNA, "protein": TYPE_PROTEIN}.get(seq_type, TYPE_UNKNOWN)


def write_bio(path: str, rows: list[SequenceRow], seq_type: str = "dna",
              mask_index: int = -1, numbering_index: int = -1, magic: bytes = DEFAULT_MAGIC):
    """Write a BioEdit Project file laid out like the ones BioEdit 7.x writes.

    `numbering_index` is BioEdit's numbering mask: the row whose own residue numbering
    is used to report positions. NeoEdit stores its pinned reference row there, which
    is the same idea (see the "Using Masks" section of the BioEdit manual).
    """
    n = len(rows)
    header = bytearray(HEADER_SIZE)
    header[0:len(magic)] = magic
    struct.pack_into("<i", header, OFF_NSEQ, n)
    struct.pack_into("<i", header, OFF_MASK, mask_index)
    struct.pack_into("<i", header, OFF_NUMMASK, numbering_index)
    body = bytearray()
    offsets = []
    base = HEADER_SIZE + 4 * n
    tcode = _seq_type_code(seq_type).encode("latin-1")
    for r in rows:
        offsets.append(base + len(body))
        for field in [r.name.encode("latin-1", "replace"),
                      r.seq.encode("latin-1", "replace"),
                      tcode] + RECORD_TAIL_FIELDS:
            body += struct.pack("<i", len(field)) + field
        body += RECORD_TRAILER
    with open(path, "wb") as fh:
        fh.write(header)
        fh.write(struct.pack(f"<{n}i", *offsets) if n else b"")
        fh.write(body)
        fh.write(FILE_TRAILER)


# --------------------------------------------------------------- GenBank dialect
TITLE_RE = re.compile(r"^TITLE\s{2,}(.+)$")
LOCUS_RE = re.compile(r"^LOCUS\s+(\S+)")


def genbank_titles(path: str) -> list[str]:
    """BioEdit writes a non-standard top-level TITLE field holding the sequence title.

    Returns one title per record (empty string where absent). REFERENCE sub-fields
    named TITLE are ignored — only top-level TITLE lines count.
    """
    titles: list[str] = []
    cur = ""
    in_reference = False
    started = False
    with open(path, "r", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("LOCUS"):
                if started:
                    titles.append(cur)
                cur = ""; started = True; in_reference = False
                continue
            if not started:
                continue
            if ln.startswith("REFERENCE"):
                in_reference = True
            elif ln[:1] not in (" ", "\t", "") and not ln.startswith("//"):
                in_reference = False
            if not in_reference:
                m = TITLE_RE.match(ln)
                if m:
                    cur = m.group(1).strip()
    if started:
        titles.append(cur)
    return titles


def read_genbank_bioedit(path: str) -> list[SequenceRow]:
    """Tolerant reader for BioEdit-style GenBank (top-level TITLE, blank line after
    ORIGIN, residues in blocks of five). Also copes with standard GenBank layout."""
    rows: list[SequenceRow] = []
    locus = definition = title = ""
    seq_parts: list[str] = []
    in_seq = False

    def flush():
        nonlocal locus, definition, title, seq_parts, in_seq
        if locus or seq_parts or title:
            name = title or definition or locus or f"sequence_{len(rows) + 1}"
            rows.append(SequenceRow(name=name.strip(), seq="".join(seq_parts),
                                    description="" if name.strip() == definition.strip() else definition.strip()))
        locus = definition = title = ""
        seq_parts = []
        in_seq = False

    with open(path, "r", errors="replace") as fh:
        for ln in fh:
            line = ln.rstrip("\n")
            if line.startswith("//"):
                flush(); continue
            if line.startswith("LOCUS"):
                if locus or seq_parts:
                    flush()
                locus = line[12:].split()[0] if len(line) > 12 else ""
                continue
            if line.startswith("DEFINITION"):
                definition = line[12:].strip(); continue
            m = TITLE_RE.match(line)
            if m:
                title = m.group(1).strip(); continue
            if line.startswith("ORIGIN"):
                in_seq = True; continue
            if in_seq and line.strip():
                body = re.sub(r"^\s*\d+", "", line)
                seq_parts.append(re.sub(r"[^A-Za-z\-\.\~\*]", "", body))
    flush()
    return [r for r in rows if r.seq or r.name]


def write_genbank_bioedit(rows: list[SequenceRow], seq_type: str = "dna",
                          line_width: int = 50, block: int = 5) -> str:
    """GenBank in BioEdit's own style: LOCUS/DEFINITION/TITLE/ORIGIN, residues in
    blocks of five, fifty per line — the layout shown in the BioEdit manual."""
    unit = "amino acids" if seq_type == "protein" else "bp"
    out = []
    for r in rows:
        seq = r.seq.replace("-", "").replace(".", "").replace("~", "")
        short = re.sub(r"\s+", "_", r.name)[:10]
        out.append(f"LOCUS       {short:<16}{len(seq)} {unit}")
        out.append(f"DEFINITION  {r.description or r.name}")
        out.append(f"TITLE       {r.name}")
        out.append("ORIGIN")
        out.append("")
        for i in range(0, len(seq), line_width):
            chunk = seq[i:i + line_width]
            groups = " ".join(chunk[j:j + block] for j in range(0, len(chunk), block))
            out.append(f"{i + 1:>6}  {groups}")
        out.append("//")
    return "\n".join(out) + "\n"
