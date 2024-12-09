from pathlib import Path

from anndata import AnnData
from scanpy import read
from scanpy.readwrite import _check_datafile_present_and_download

from scmorph.utils import _doc_params

HERE = Path(__file__).parent


_backed_docstring = """\
backed
    If ``'r'``, load AnnData in ``backed`` mode instead
    of fully loading it into memory (`memory` mode). If you want to modify
    backed attributes of the AnnData object, you need to choose ``'r+'``.
"""


def rohban2017_minimal_csv() -> Path:
    """Provides a minimal csv file in CellProfiler format, data from :cite:t:`Rohban2017`"""
    filename = HERE / "rohban2017_CellProfiler_minimal.csv"
    backup = "https://figshare.com/ndownloader/files/50656098"
    _check_datafile_present_and_download(filename, backup_url=backup)
    return filename


@_doc_params(backed=_backed_docstring)
def rohban2017_minimal(backed=None) -> AnnData:
    """\
    Load a subset of a multi-plate experiment by :cite:t:`Rohban2017` with ~12,000 cells

    Parameters
    ----------
    {backed}
    """
    filename = HERE / "rohban2017_subset.h5ad"
    backup = "https://figshare.com/ndownloader/files/50656878"
    return read(filename, backup_url=backup, backed=backed)


@_doc_params(backed=_backed_docstring)
def rohban2017(backed=None) -> AnnData:
    """\
    Load a large multi-plate experiment by :cite:t:`Rohban2017` with ~1.2M cells

    Parameters
    ----------
    {backed}
    """
    filename = HERE / "rohban2017_subset.h5ad"
    backup = "https://figshare.com/ndownloader/files/50650236"
    return read(filename, backup_url=backup, backed=backed)


@_doc_params(backed=_backed_docstring)
def rohban2017_imageQC(backed=None) -> AnnData:
    """\
    Load image-level data for a multi-plate experiment by :cite:t:`Rohban2017`

    Parameters
    ----------
    {backed}
    """
    filename = HERE / "rohban2017_imageQC.h5ad"
    backup = "https://figshare.com/ndownloader/files/50651907"
    return read(filename, backup_url=backup, backed=backed)
