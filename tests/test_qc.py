import warnings

import numpy as np
import pandas as pd
import pytest

import scmorph as sm

data_nrows_no_na = 12286
data_nrows_qc = 11126


@pytest.fixture
def adata_var_filtered(adata_no_na):
    adata = adata_no_na
    pass_var = np.empty(len(adata.var), dtype=bool)
    for i, feat in enumerate(adata.var_names):
        pass_var[i] = False if np.var(adata[:, feat].X) < 1e-5 else True
    return adata[:, pass_var].copy()


@pytest.fixture
def adata_imageQC():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sm.datasets.rohban2017_imageQC()
    adata.obs.index = pd.Series(np.arange(adata.shape[0])).astype(str)
    return adata


def test_filter_outliers(adata_var_filtered):
    # Outlier detection requires features with non-zero variance
    assert (
        sm.qc.filter_outliers(adata_var_filtered, n_obs=1000, outliers=0.1).shape[0]
        == data_nrows_qc
    )


def test_count_cells_per_group(adata_no_na):
    adata = adata_no_na
    sm.qc.count_cells_per_group(
        adata, ["Image_Metadata_Plate", "Image_Metadata_Well", "Image_Metadata_Site"]
    )
    assert adata.obs["cells_per_group"].describe().loc["50%"] == 76.0


def test_qc_images_by_dissimilarity(adata_no_na, adata_imageQC):
    adata = adata_no_na
    sm.qc.qc_images_by_dissimilarity(adata, adata_imageQC, threshold=0.2)
    assert adata.obs["PassQC"].value_counts().loc["True"] == 11584
