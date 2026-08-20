"""Network importers (Qt-free): NCBI Entrez E-utilities and the Ensembl REST API.

Both clients fetch records as text, write them to a local *downloads* folder as
GenBank/FASTA, and hand the path back to the UI, which opens it through the normal
file path — so features, the gene view, circular topology and "Open recent" all work
exactly as for a file on disk.
"""
from .http import RemoteError, http_get
from . import ncbi, ensembl

__all__ = ["RemoteError", "http_get", "ncbi", "ensembl"]
