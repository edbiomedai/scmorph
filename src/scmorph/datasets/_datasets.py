from pathlib import Path

from anndata import AnnData
from scanpy import read

HERE = Path(__file__).parent


def multi_plate_experiment() -> AnnData:
    """Load a multi-plate experiment of ~1,600 cells"""
    filename = HERE / "test_data.h5ad"
    backup = "https://uoe-my.sharepoint.com/:u:/g/personal/s2221912_ed_ac_uk/Eb8nsk2jHZhLrUDrMow5hykBYn28Q5AjPtxIe-0t6u1sWQ?download=1"
    return read(filename, backup_url=backup)


def image_qc_data() -> AnnData:
    """Load image QC data"""
    filename = HERE / "test_image_data.h5ad"
    backup = "https://uoe-my.sharepoint.com/:u:/g/personal/s2221912_ed_ac_uk/EZWSp1ys8exJvxZ-10SaRPMBqXEyIZYKa_PEuhPsqEMSWQ?download=1"
    return read(filename, backup_url=backup)
