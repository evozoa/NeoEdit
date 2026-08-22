"""UCSC Genome Browser REST API (api.genome.ucsc.edu): assemblies, sequence, gene-model tracks.

Gene models come as genePred rows from one of the gene tracks; they are converted into the
Ensembl-style feature dicts that `ensembl.build_record` already turns into a GenBank record,
so UCSC imports get the same gene / mRNA / CDS features (UCSC exonFrames -> GFF phase)."""
from __future__ import annotations

import json
import re
import urllib.parse
from dataclasses import dataclass

from Bio.SeqRecord import SeqRecord

from .http import http_get, RemoteError
from . import ensembl as E

API = "https://api.genome.ucsc.edu"
# tried in this order when the user asks for "automatic" gene models; not every assembly has every track
GENE_TRACKS = ["ncbiRefSeq", "refGene", "ensGene", "knownGene", "augustusGene"]
TRACK_LABELS = [("Automatic (RefSeq → Ensembl → GENCODE → Augustus)", "auto"),
                ("NCBI RefSeq (ncbiRefSeq)", "ncbiRefSeq"), ("RefSeq genes (refGene)", "refGene"),
                ("Ensembl genes (ensGene)", "ensGene"), ("GENCODE / UCSC genes (knownGene)", "knownGene"),
                ("Augustus predictions (augustusGene)", "augustusGene")]
# tracks whose search hits name a gene span, most trustworthy first
_SEARCH_PRIORITY = ["hgnc", "refGene", "ncbiRefSeqCurated", "ncbiRefSeq", "knownGene", "ensGene", "mane", "augustusGene"]


@dataclass
class Genome:
    id: str                      # hg38, danRer11 …
    organism: str                # Human
    scientific_name: str         # Homo sapiens
    description: str             # Dec. 2013 (GRCh38/hg38)
    active: bool = True
    order: int = 0

    def label(self) -> str:
        return f"{self.id} — {self.organism} ({self.description})"


@dataclass
class Hit:
    name: str
    chrom: str
    start: int                   # 1-based inclusive
    end: int
    track: str
    description: str = ""

    def region(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


class UCSCClient:
    def __init__(self, fetch=http_get, api: str = API):
        self._fetch = fetch
        self.api = api.rstrip("/")

    def _get(self, path: str, **params):
        q = ";".join(f"{k}={urllib.parse.quote(str(v), safe='')}" for k, v in params.items() if v is not None)
        url = f"{self.api}{path}" + (f"?{q}" if q else "")
        body = self._fetch(url, headers={"Accept": "application/json"}, timeout=120)
        try:
            d = json.loads(body)
        except ValueError:
            raise RemoteError(f"UCSC returned non-JSON for {path}: {body[:120]!r}") from None
        if isinstance(d, dict) and d.get("error"):
            raise RemoteError(f"UCSC: {d['error']}")
        return d

    # ------------------------------------------------------------ API
    def genomes(self) -> list[Genome]:
        d = self._get("/list/ucscGenomes").get("ucscGenomes", {})
        out = []
        for gid, g in d.items():
            out.append(Genome(gid, g.get("organism", ""), g.get("scientificName", ""), g.get("description", ""),
                              bool(g.get("active", 1)), int(g.get("orderKey", 0) or 0)))
        out.sort(key=lambda g: (g.order, g.id))
        return out

    def chromosomes(self, genome: str) -> dict[str, int]:
        return dict(self._get("/list/chromosomes", genome=genome).get("chromosomes", {}))

    def sequence(self, genome: str, chrom: str, start0: int, end0: int) -> str:
        """DNA for [start0, end0) (0-based half-open), soft-masked as served by UCSC."""
        if end0 <= start0:
            raise RemoteError("Empty region.")
        out = []
        step = 8_000_000
        for s in range(start0, end0, step):
            d = self._get("/getData/sequence", genome=genome, chrom=chrom, start=s, end=min(end0, s + step))
            out.append(d.get("dna", ""))
        seq = "".join(out)
        if not seq:
            raise RemoteError(f"UCSC returned no sequence for {genome} {chrom}:{start0 + 1}-{end0}.")
        return seq

    def track(self, genome: str, track: str, chrom: str, start0: int, end0: int) -> list[dict]:
        d = self._get("/getData/track", genome=genome, track=track, chrom=chrom, start=start0, end=end0)
        items = d.get(track, [])
        if isinstance(items, dict):                       # some composites key by chromosome
            items = items.get(chrom, []) or [x for v in items.values() for x in v]
        return [x for x in items if "txStart" in x]

    def gene_models(self, genome: str, chrom: str, start0: int, end0: int, track: str = "auto") -> tuple[str, list[dict]]:
        """genePred rows overlapping the window from `track`, or the first available track for 'auto'."""
        tracks = GENE_TRACKS if track == "auto" else [track]
        last = None
        for t in tracks:
            try:
                rows = self.track(genome, t, chrom, start0, end0)
            except RemoteError as e:
                last = e
                continue
            if rows or track != "auto":
                return t, rows
        if track != "auto" and last is not None:
            raise last
        return "", []

    def search(self, genome: str, term: str) -> Hit:
        """Gene symbol / transcript id -> best positional hit (1-based inclusive)."""
        term = (term or "").strip()
        if not term:
            raise RemoteError("Enter a gene symbol or accession.")
        d = self._get("/search", genome=genome, search=term)
        groups = d.get("positionMatches", []) or []
        hits: list[tuple[int, int, Hit]] = []
        low = term.lower()
        for g in groups:
            tname = g.get("trackName") or g.get("name") or ""
            prio = _SEARCH_PRIORITY.index(tname) if tname in _SEARCH_PRIORITY else len(_SEARCH_PRIORITY)
            for m in g.get("matches", []):
                pos = m.get("position", "")
                mm = re.match(r"^([\w.-]+):([\d,]+)-([\d,]+)$", pos)
                if not mm:
                    continue
                pn = (m.get("posName") or "").strip()
                pnl = pn.lower()
                if pnl == low:
                    rank = 0
                elif pnl.startswith(low + " (") or pnl.split(" (")[0].lower() == low:
                    rank = 1
                elif pnl.split(".")[0] == low.split(".")[0]:
                    rank = 2
                else:
                    rank = 3
                hits.append((rank, prio, Hit(pn or term, mm.group(1), int(mm.group(2).replace(",", "")),
                                             int(mm.group(3).replace(",", "")), tname, m.get("description") or "")))
        if not hits:
            raise RemoteError(f"UCSC: no match for '{term}' in {genome}.")
        hits.sort(key=lambda h: (h[0], h[1]))
        return hits[0][2]

    # ------------------------------------------------------------ high level
    def fetch_genomic(self, genome: str, chrom: str, start: int, end: int, strand: int = 1, annotate: bool = True,
                      track: str = "auto", keep_mask: bool = False, name: str | None = None,
                      species: str = "", hit: Hit | None = None) -> SeqRecord:
        """Region (1-based inclusive) as an annotated SeqRecord. Returns the record; the track
        actually used is stored in rec.annotations['ucsc_track']."""
        seq = self.sequence(genome, chrom, start - 1, end)
        if not keep_mask:
            seq = seq.upper()
        used, feats = "", []
        if annotate:
            used, rows = self.gene_models(genome, chrom, start - 1, end, track)
            feats = genepred_to_features(rows)
        gene = None
        if hit is not None:
            gene = E.Lookup(id=hit.name, object_type="Gene", species=species, seq_region_name=hit.chrom,
                            start=hit.start, end=hit.end, strand=strand, display_name=hit.name,
                            biotype="", description=hit.description)
        rec = E.build_record(seq, chrom, start, end, feats, species=species, assembly=genome, strand=strand,
                             name=name, gene=gene, masked=keep_mask)
        rec.description = rec.description.replace("Ensembl", "UCSC").strip()
        if used:
            rec.description += f" UCSC {genome} {used}"
        else:
            rec.description += f" UCSC {genome}"
        rec.annotations["comment"] = f"Imported by NeoEdit from the UCSC Genome Browser ({genome}, {chrom}:{start}-{end}" \
                                     f"{', gene models from ' + used if used else ''})"
        rec.annotations["ucsc_track"] = used
        return rec


def strand_of_hit(hit: Hit, rows: list[dict]) -> int:
    """Strand of the gene model that best matches a search hit: a name match wins (names may
    differ in prefix, e.g. MT-ND6 vs ND6), otherwise the model overlapping the hit span most."""
    if not rows:
        return 1
    hn = hit.name.lower()
    alt = hn[3:] if hn.startswith("mt-") else hn
    best, best_key = None, None
    for r in rows:
        n2 = str(r.get("name2") or "").lower(); n1 = str(r.get("name") or "").lower()
        name_ok = hn in (n2, n1, n1.split(".")[0]) or alt == n2
        ov = min(int(r["txEnd"]), hit.end) - max(int(r["txStart"]), hit.start - 1)
        key = (1 if name_ok else 0, ov)
        if best_key is None or key > best_key:
            best, best_key = r, key
    return -1 if best.get("strand") == "-" else 1


# ---------------------------------------------------------------- genePred -> Ensembl-style features
def _ints(csv: str) -> list[int]:
    return [int(x) for x in str(csv).strip().rstrip(",").split(",") if x.strip() != ""]


def genepred_to_features(rows: list[dict]) -> list[dict]:
    """Convert genePred rows (0-based half-open) to the 1-based dicts `ensembl.build_record` expects:
    gene (grouped by name2/strand), transcript, exon and cds items; UCSC exonFrames become GFF phases."""
    genes: dict[tuple[str, str], dict] = {}
    feats: list[dict] = []
    for t in rows:
        strand = -1 if t.get("strand") == "-" else 1
        tname = str(t.get("name") or "")
        gname = str(t.get("name2") or "") or tname.split(".")[0]
        ts, te = int(t["txStart"]), int(t["txEnd"])
        cs, ce = int(t.get("cdsStart", ts)), int(t.get("cdsEnd", ts))
        coding = ce > cs
        gkey = (gname, t.get("strand", "+"))
        g = genes.get(gkey)
        if g is None:
            g = genes[gkey] = {"feature_type": "gene", "id": gname, "external_name": gname, "start": ts + 1, "end": te,
                               "strand": strand, "biotype": "protein_coding" if coding else "ncRNA", "description": ""}
            feats.append(g)
        else:
            g["start"] = min(g["start"], ts + 1); g["end"] = max(g["end"], te)
            if coding:
                g["biotype"] = "protein_coding"
        tid = tname
        tver = None
        m = re.match(r"^(.*)\.(\d+)$", tname)
        if m:
            tid, tver = m.group(1), m.group(2)
        tr = {"feature_type": "transcript", "id": tid, "Parent": gname, "start": ts + 1, "end": te, "strand": strand,
              "biotype": "protein_coding" if coding else "ncRNA", "external_name": tname, "is_canonical": 0}
        if tver:
            tr["version"] = tver
        feats.append(tr)
        starts, ends = _ints(t.get("exonStarts", "")), _ints(t.get("exonEnds", ""))
        frames = _ints(t.get("exonFrames", "")) if t.get("exonFrames") else []
        for i, (es, ee) in enumerate(zip(starts, ends)):
            feats.append({"feature_type": "exon", "id": f"{tname}.e{i + 1}", "Parent": tid, "start": es + 1, "end": ee,
                          "strand": strand, "rank": i + 1})
            if coding:
                a, b = max(es, cs), min(ee, ce)
                if b > a:
                    fr = frames[i] if i < len(frames) else -1
                    phase = (3 - fr) % 3 if fr >= 0 else 0
                    feats.append({"feature_type": "cds", "id": tid, "Parent": tid, "start": a + 1, "end": b,
                                  "strand": strand, "phase": phase,
                                  "protein_id": tname if re.match(r"^[NXY]P_", tname) else ""})
    return feats
