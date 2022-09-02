from pathlib import Path

from anndata import AnnData
from scanpy import read

HERE = Path(__file__).parent


def multi_plate_experiment() -> AnnData:
    """Load a multi-plate experiment of ~1,600 cells"""
    filename = HERE / "test_data.h5ad"
    backup = "https://github.com/edbiomedai/scmorph/raw/launch_prep/src/scmorph/datasets/test_data.h5ad"
    return read(filename, backup_url=backup)


def image_qc_data() -> AnnData:
    """Load image QC data"""
    filename = HERE / "test_image_data.h5ad"
    backup = "https://github.com/edbiomedai/scmorph/raw/launch_prep/src/scmorph/datasets/test_image_data.h5ad"
    return read(filename, backup_url=backup)
