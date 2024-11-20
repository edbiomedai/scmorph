from pathlib import Path

from anndata import AnnData
from scanpy import read
from scanpy.readwrite import _check_datafile_present_and_download

HERE = Path(__file__).parent


# def multi_plate_experiment() -> AnnData:
#     """Load a multi-plate experiment of ~1,600 cells"""
#     filename = HERE / "test_data.h5ad"
#     backup = "https://github.com/edbiomedai/scmorph/raw/launch_prep/src/scmorph/datasets/test_data.h5ad"
#     return read(filename, backup_url=backup)


def rohban2018_minimal_csv() -> Path:
    """Provides a minimal csv file in CellProfiler format, data from Rohban et al. 2018"""
    filename = HERE / "rohban2018_CellProfiler_minimal.csv"
    backup = "https://figshare.com/ndownloader/files/50656098"
    _check_datafile_present_and_download(filename, backup_url=backup)
    return filename


def rohban2018_minimal(**kwargs) -> AnnData:
    """Load a subset of a multi-plate experiment by Rohban with ~12,000 cells"""
    filename = HERE / "rohban2018_subset.h5ad"
    backup = "https://figshare.com/ndownloader/files/50656878"
    return read(filename, backup_url=backup, **kwargs)


def rohban2018(**kwargs) -> AnnData:
    """Load a large multi-plate experiment by Rohban with ~1.2M cells"""
    filename = HERE / "rohban2018_subset.h5ad"
    backup = "https://figshare.com/ndownloader/files/50650236"
    return read(filename, backup_url=backup, **kwargs)


def rohban2018_imageQC(**kwargs) -> AnnData:
    """Load image-level data for a multi-plate experiment by Rohban"""
    filename = HERE / "rohban2018_imageQC.h5ad"
    backup = "https://figshare.com/ndownloader/files/50651907"
    return read(filename, backup_url=backup, **kwargs)


# def image_qc_data() -> AnnData:
#     """Load image QC data"""
#     filename = HERE / "test_image_data.h5ad"
#     backup = "https://github.com/edbiomedai/scmorph/raw/launch_prep/src/scmorph/datasets/test_image_data.h5ad"
#     return read(filename, backup_url=backup)
