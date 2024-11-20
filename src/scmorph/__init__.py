from importlib import metadata

from packaging.version import Version
from scanpy import read_h5ad, write

from scmorph import datasets, logging, pl, pp, qc, tl
from scmorph.io import read, read_cellprofiler_batches, read_cellprofiler_csv, read_sql

try:
    md = metadata.metadata(__name__)
    __version__ = md.get("version", "")
    __author__ = md.get("Author", "")
    __maintainer__ = md.get("Maintainer-email", "")
except ImportError:
    md = None  # type: ignore[assignment]
