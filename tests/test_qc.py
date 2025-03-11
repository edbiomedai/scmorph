import numpy as np
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


def test_filter_outliers(adata_var_filtered):
    # Outlier detection requires features with non-zero variance
    assert (
        sm.qc.filter_outliers(adata_var_filtered, n_obs=1000, outliers=0.1).shape[0]
        == data_nrows_qc
    )
