"""Build NC_012920_MDP.gb: the human mitochondrial reference (rCRS) with the eight
described mitochondrial-derived peptides added as CDS features.

MDPs are translated by cytoplasmic ribosomes, so each MDP feature carries
/transl_table=1 while the canonical mitochondrial genes keep /transl_table=2.
Coordinates were derived by locating each UniProt peptide in the rCRS six-frame
standard-code translation (see the assertion below, which re-verifies them).

Run:  python examples/mito/make_rcrs_mdp.py
"""
from __future__ import annotations

import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "NC_012920.gb")
OUT = os.path.join(HERE, "NC_012920_MDP.gb")

# name, start, end (1-based inclusive, rCRS), strand, host gene, UniProt, peptide
MDPS = [
    ("MOTS-c", 1343, 1390, +1, "MT-RNR1", "A0A0C5B5G6", "MRWQEMGYIFYPRKLR"),
    ("SHLP3", 1706, 1819, -1, "MT-RNR2", "A0A3G1DJQ2", "MLGYNFSSFPCGTISIAPGFNFYRLYFIWVNGLAKVVW"),
    ("SHLP2", 2091, 2168, -1, "MT-RNR2", "A0A3G1DIU6", "MGVKFFTLSTRFFPSVQRAVPLWTNS"),
    ("SHLP4", 2445, 2522, -1, "MT-RNR2", "A0A3G1DJK2", "MLEVMFLVNRRGKICRVPFTFFNLSL"),
    ("SHLP1", 2488, 2559, -1, "MT-RNR2", "A0A3G1DJL7", "MCHWAGGASNTGDARGDVFGKQAG"),
    ("humanin", 2633, 2704, +1, "MT-RNR2", "Q8IVG9", "MAPRGFSCLLLLTSEIDLPVKRRA"),
    ("SHLP5", 2783, 2854, -1, "MT-RNR2", "A0A3G1DJL1", "MYCSEVGFCSEVAPTEIFNAGLVV"),
    ("SHLP6", 2990, 3049, +1, "MT-RNR2", "A0A3G1DJN1", "MLDQDIPMVQPLLKVRLFND"),
]

NOTE = ("mitochondrial-derived peptide (MDP); encoded in mtDNA but translated by "
        "cytoplasmic ribosomes, so the standard genetic code (transl_table=1) applies, "
        "unlike the surrounding mitochondrially-translated genes")


def main():
    if not os.path.exists(SRC):
        sys.exit(f"missing {SRC}; fetch NC_012920.1 in GenBank format first")
    rec = SeqIO.read(SRC, "genbank")
    added = 0
    for name, start, end, strand, host, acc, pep in MDPS:
        loc = SimpleLocation(start - 1, end, strand=strand)
        nt = loc.extract(rec.seq)
        translated = str(Seq(str(nt)).translate(table=1))
        translated = translated[:-1] if translated.endswith("*") else translated
        assert translated == pep, f"{name}: rCRS gives {translated!r}, expected {pep!r}"
        rec.features.append(SeqFeature(loc, type="CDS", qualifiers={
            "gene": [name],
            "product": [name],
            "codon_start": ["1"],
            "transl_table": ["1"],
            "translation": [pep],
            "db_xref": [f"UniProtKB/Swiss-Prot:{acc}"],
            "note": [NOTE, f"located within {host}", f"{len(pep)} aa"],
        }))
        added += 1
    rec.features.sort(key=lambda f: (int(f.location.start), -int(f.location.end)))
    rec.annotations["comment"] = ((rec.annotations.get("comment", "") + "\n\n") +
                                  f"{added} mitochondrial-derived peptides added for NeoEdit "
                                  f"(CDS features with /transl_table=1; canonical genes keep "
                                  f"/transl_table=2).").strip()
    SeqIO.write(rec, OUT, "genbank")
    print(f"wrote {OUT}: {len(rec.features)} features ({added} MDPs)")


if __name__ == "__main__":
    main()
