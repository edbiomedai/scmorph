import numpy as np
import pandas as pd
import pytest

import scmorph as sm

data_nrows_no_na = 12286
data_nrows_qc = 11126


@pytest.fixture
def adata():
    adata = sm.datasets.rohban2017_minimal()
    sm.pp.drop_na(adata)
    return adata


@pytest.fixture
def adata_var_filtered(adata):
    pass_var = np.empty(len(adata.var), dtype=bool)
    for i, feat in enumerate(adata.var_names):
        pass_var[i] = False if np.var(adata[:, feat].X) < 1e-5 else True
    return adata[:, pass_var].copy()


@pytest.fixture
def adata_imageQC():
    adata = sm.datasets.rohban2017_imageQC()
    adata.obs.index = pd.Series(np.arange(adata.shape[0])).astype(str)
    return adata


def test_filter_outliers(adata_var_filtered):
    # Outlier detection requires features with non-zero variance
    assert (
        sm.qc.filter_outliers(adata_var_filtered, n_obs=1000, outliers=0.1).shape[0]
        == data_nrows_qc
    )


def test_count_cells_per_group(adata):
    sm.qc.count_cells_per_group(
        adata, ["Image_Metadata_Plate", "Image_Metadata_Well", "Image_Metadata_Site"]
    )
    assert adata.obs["CellsPerGroup"].describe().loc["50%"] == 77.0


def test_qc_images_by_dissimilarity(adata, adata_imageQC):
    sm.qc.qc_images_by_dissimilarity(adata, adata_imageQC, threshold=0.2)
    assert adata.obs["PassQC"].value_counts().loc["True"] == 11584
