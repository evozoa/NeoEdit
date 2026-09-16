"""NCBI Entrez E-utilities: esearch / esummary / efetch for nuccore and protein.

No API key is required; NCBI asks clients to identify themselves (tool + e-mail) and to
stay under 3 requests/s (10/s with a key), which `NCBIClient` enforces."""
from __future__ import annotations

import json
import re
import time
import urllib.parse
from dataclasses import dataclass, field

from .http import http_get, RemoteError, write_download

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATABASES = [("Nucleotide (nuccore)", "nuccore"), ("Protein", "protein")]
# (label, value) -> efetch rettype per database
FORMATS = [("GenBank (with features)", "gb"), ("FASTA (sequence only)", "fasta")]
_RETTYPE = {("nuccore", "gb"): "gbwithparts", ("nuccore", "fasta"): "fasta",
            ("protein", "gb"): "gp", ("protein", "fasta"): "fasta"}
_EXT = {"gb": ".gb", "fasta": ".fasta"}

# Fields an imported sequence can be named from (key, label). All come from esummary: the
# organism, title and length, plus the record's source qualifiers (isolate, voucher, ...).
NAME_FIELDS = [
    ("accession", "Accession.version"),
    ("accession_base", "Accession (no version)"),
    ("organism", "Organism (full name)"),
    ("genus", "Genus"),
    ("species", "Species (epithet)"),
    ("isolate", "Isolate"),
    ("strain", "Strain"),
    ("specimen_voucher", "Specimen voucher"),
    ("clone", "Clone"),
    ("haplotype", "Haplotype"),
    ("cultivar", "Cultivar"),
    ("serotype", "Serotype"),
    ("segment", "Segment"),
    ("host", "Host"),
    ("country", "Country"),
    ("location", "Location (full)"),
    ("lat_lon", "Latitude / longitude"),
    ("collection_date", "Collection date"),
    ("collected_by", "Collected by"),
    ("genome", "Organelle / genome"),
    ("length", "Length"),
    ("taxid", "Taxonomy ID"),
    ("title", "Definition line (title)"),
]
NAME_LABELS = dict(NAME_FIELDS)
# NCBI's own header: what NeoEdit shows when nothing is customised
DEFAULT_NAME_FIELDS = ["accession", "title"]
NAME_SEPARATORS = [("Space", " "), ("Underscore  _", "_"), ("Pipe  |", "|"), ("Hyphen  -", "-")]
_SOURCE_KEYS = {"isolate", "strain", "specimen_voucher", "clone", "haplotype", "cultivar", "serotype",
                "segment", "host", "lat_lon", "collection_date", "collected_by"}


@dataclass
class Summary:
    uid: str
    accession: str
    title: str
    length: int
    organism: str
    moltype: str = ""
    topology: str = ""
    extra: dict = field(default_factory=dict)
    source: dict = field(default_factory=dict)     # source qualifiers: {"isolate": "...", ...}
    taxid: str = ""
    genome: str = ""                               # "mitochondrion", "chloroplast", "genomic", ...


def parse_ids(text: str) -> list[str]:
    """Split user input into unique accession/GI tokens (whitespace, commas, semicolons)."""
    out, seen = [], set()
    for tok in re.split(r"[\s,;]+", text or ""):
        tok = tok.strip().strip(">")
        if tok and tok not in seen:
            seen.add(tok); out.append(tok)
    return out


class NCBIClient:
    def __init__(self, email: str = "", api_key: str = "", fetch=http_get):
        self.email = (email or "").strip()
        self.api_key = (api_key or "").strip()
        self._fetch = fetch
        self._last = 0.0

    # ------------------------------------------------------------ plumbing
    def _params(self, **kw) -> dict:
        p = {"tool": "neoedit"}
        if self.email:
            p["email"] = self.email
        if self.api_key:
            p["api_key"] = self.api_key
        p.update({k: v for k, v in kw.items() if v not in (None, "")})
        return p

    def _throttle(self):
        gap = 0.11 if self.api_key else 0.35
        dt = time.monotonic() - self._last
        if dt < gap:
            time.sleep(gap - dt)
        self._last = time.monotonic()

    def _call(self, endpoint: str, post: bool = False, **kw) -> str:
        self._throttle()
        params = self._params(**kw)
        q = urllib.parse.urlencode(params)
        url = EUTILS + endpoint
        if post:
            body = self._fetch(url, data=q.encode(), headers={"Content-Type": "application/x-www-form-urlencoded"})
        else:
            body = self._fetch(url + "?" + q)
        text = body.decode("utf-8", "replace")
        t = text.lstrip()
        if t.startswith("Error") or t.startswith("<!DOCTYPE") or t.startswith("<html"):
            # efetch answers HTTP 200 with a text error for unknown ids
            raise RemoteError(f"NCBI: {_unspace(t[:300])}")
        return text

    # ------------------------------------------------------------ API
    def search(self, db: str, term: str, retmax: int = 200) -> tuple[int, list[str]]:
        """esearch: returns (total hit count, list of UIDs up to retmax)."""
        txt = self._call("esearch.fcgi", db=db, term=term, retmax=retmax, retmode="json")
        d = json.loads(txt).get("esearchresult", {})
        if "ERROR" in d:
            raise RemoteError(f"NCBI search: {d['ERROR']}")
        return int(d.get("count", 0)), list(d.get("idlist", []))

    def summaries(self, db: str, uids: list[str]) -> list[Summary]:
        out: list[Summary] = []
        for i in range(0, len(uids), 200):
            chunk = uids[i:i + 200]
            txt = self._call("esummary.fcgi", post=len(chunk) > 20, db=db, id=",".join(chunk), retmode="json")
            res = json.loads(txt).get("result", {})
            for uid in res.get("uids", []):
                r = res.get(uid, {})
                if "error" in r:
                    continue
                out.append(Summary(uid=uid, accession=r.get("accessionversion") or r.get("caption", uid),
                                   title=r.get("title", ""), length=int(r.get("slen") or 0),
                                   organism=r.get("organism", ""), moltype=r.get("moltype", ""),
                                   topology=r.get("topology", ""),
                                   extra={k: r[k] for k in ("biomol", "sourcedb", "completeness", "geneticcode", "updatedate")
                                          if k in r},
                                   source=parse_source(r.get("subtype", ""), r.get("subname", "")),
                                   taxid=str(r.get("taxid") or ""), genome=r.get("genome", "")))
        return out

    def namer(self, db: str, ids: list[str], fields: list[str], sep: str = " ",
              spaces: bool = False):
        """A function record id -> sequence name built from `fields` (see `format_name`), from
        the esummaries of `ids`; it returns None for records it has no summary for."""
        by_acc: dict[str, Summary] = {}
        for s in self.summaries(db, ids):
            by_acc[s.accession] = s
            by_acc.setdefault(s.accession.split(".")[0], s)
            by_acc.setdefault(s.uid, s)

        def name(record_id: str) -> str | None:
            key = record_id.split(":")[0]
            s = by_acc.get(key) or by_acc.get(key.split(".")[0])
            return format_name(name_fields(s, record_id), fields, sep, spaces) if s else None
        return name

    def fetch(self, db: str, ids: list[str], fmt: str = "gb", seq_start: int | None = None,
              seq_stop: int | None = None, strand: int | None = None) -> str:
        """efetch the records as text. A sub-range/strand only applies to a single id."""
        if not ids:
            raise RemoteError("No accession given.")
        rettype = _RETTYPE.get((db, fmt))
        if rettype is None:
            raise RemoteError(f"Unsupported database/format: {db}/{fmt}")
        kw = dict(db=db, id=",".join(ids), rettype=rettype, retmode="text")
        if len(ids) == 1:
            if seq_start:
                kw["seq_start"] = int(seq_start)
            if seq_stop:
                kw["seq_stop"] = int(seq_stop)
            if strand in (1, 2):
                kw["strand"] = strand
        text = self._call("efetch.fcgi", post=len(ids) > 20, **kw)
        if not text.strip():
            raise RemoteError("NCBI returned an empty record — check the accession and database (nucleotide vs protein).")
        return text

    def download(self, db: str, ids: list[str], fmt: str, out_dir: str, **range_kw) -> tuple[str, str]:
        """Fetch and save to out_dir. Returns (path, text)."""
        text = self.fetch(db, ids, fmt, **range_kw)
        stem = ids[0] if len(ids) == 1 else f"{ids[0]}_and_{len(ids) - 1}_more"
        rng = range_kw.get("seq_start"), range_kw.get("seq_stop")
        if len(ids) == 1 and any(rng):
            stem += f"_{rng[0] or 1}-{rng[1] or 'end'}"
        path = write_download(text, out_dir, stem, _EXT[fmt])
        return path, text


def parse_source(subtype: str, subname: str) -> dict:
    """esummary's parallel "subtype"/"subname" lists ("isolate|country", "X1|Peru") -> dict."""
    keys = [k.strip() for k in (subtype or "").split("|")]
    vals = [v.strip() for v in (subname or "").split("|")]
    if not subtype or len(keys) != len(vals):
        return {}
    return {k: v for k, v in zip(keys, vals) if k and v}


def name_fields(s: Summary, record_id: str = "") -> dict[str, str]:
    """The NAME_FIELDS values of one record. `record_id` (the id in the downloaded file) wins over
    the summary's accession, so a sub-range keeps its ":577-647" suffix."""
    acc = record_id or s.accession
    org = " ".join(s.organism.split())
    words = org.split()
    genus = species = ""
    # a Latin binomial ("Pimephales promelas", "Pimephales sp."), not "Severe acute respiratory ..."
    if (len(words) >= 2 and words[0][:1].isupper() and words[0].isalpha() and words[1].islower()
            and "virus" not in org.lower() and "phage" not in org.lower()):
        genus, species = words[0], words[1]
    location = s.source.get("geo_loc_name") or s.source.get("country", "")
    f = {
        "accession": acc,
        "accession_base": re.sub(r"\.\d+", "", acc, count=1),
        "organism": org,
        "genus": genus,
        "species": species,
        "country": location.split(":")[0].strip(),
        "location": location,
        "genome": s.genome if s.genome not in ("", "genomic") else "",
        "length": str(s.length) if s.length else "",
        "taxid": s.taxid,
        "title": s.title,
    }
    f.update({k: v for k, v in s.source.items() if k in _SOURCE_KEYS})
    return f


def format_name(values: dict[str, str], fields: list[str], sep: str = " ", spaces: bool = False) -> str:
    """Join the `fields` that have a value, in order. `spaces`: also replace whitespace inside the
    values with `sep`. Falls back to the accession when nothing is left."""
    parts = []
    for k in fields:
        v = " ".join((values.get(k) or "").split())
        if not v:
            continue
        if spaces and sep.strip():
            v = v.replace(" ", sep)
        parts.append(v)
    return (sep or " ").join(parts) or values.get("accession", "")


def _unspace(msg: str) -> str:
    """efetch letter-spaces its errors ('F a i l e d  t o  …'): collapse that back into words."""
    msg = re.sub(r"[\r\n]+", " ", msg).strip()
    m = re.match(r"^(Error:?)\s*((?:\S\s){4,}.*)$", msg)
    if not m:
        return msg
    words = re.split(r"\s{2,}", m.group(2).strip())
    return m.group(1) + " " + " ".join(w.replace(" ", "") for w in words)


def count_records(text: str, fmt: str) -> int:
    if fmt == "fasta":
        return sum(1 for ln in text.splitlines() if ln.startswith(">"))
    return sum(1 for ln in text.splitlines() if ln.startswith("LOCUS"))
