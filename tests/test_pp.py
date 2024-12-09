import pandas as pd
import pytest

import scmorph as sm

data_nrows_no_na = 12286
data_nrows = 12352


@pytest.fixture
def adata():
    return sm.datasets.rohban2017_minimal()


@pytest.fixture
def adata_treat():
    adata = sm.datasets.rohban2017_minimal()
    adata.obs["TARGETGENE"] = adata.obs["TARGETGENE"].astype(str).replace("nan", "DMSO")
    return adata


def test_drop_na(adata):
    assert sm.pp.drop_na(adata, inplace=False).shape == (data_nrows_no_na, 1687)


def test_pca(adata):
    sm.pp.drop_na(adata)
    sm.pp.scale(adata)
    sm.pp.pca(adata, n_comps=2)
    assert adata.obsm["X_pca"].shape == (data_nrows_no_na, 2)


def test_aggregate_modes(adata):
    modes = ["mean", "median", "std", "var", "sem", "mad", "mad_scaled"]
    for m in modes:
        agg = sm.pp.aggregate(
            adata, method=m, group_keys=["Image_Metadata_Plate"], well_key="Image_Metadata_Well", progress=False
        )
        assert agg.shape == (20, adata.shape[1])


def test_aggregate_mahalanobis(adata_treat):
    agg = sm.pp.aggregate_mahalanobis(adata_treat, treatment_key="TARGETGENE", well_key="Image_Metadata_Well")
    assert agg.shape == (1,)


def test_aggregate_pc(adata_treat):
    sm.pp.drop_na(adata_treat)
    agg = sm.pp.aggregate_pc(adata_treat, treatment_key="TARGETGENE")
    assert agg.shape == (2,)


@pytest.mark.filterwarnings("ignore:Precision loss")
def test_aggregate_tstat(adata_treat):
    agg = sm.pp.aggregate_ttest(adata_treat, treatment_key="TARGETGENE")
    assert agg[0].shape == (adata_treat.shape[1], 1)  # one treatment in test data
    assert 0 <= agg[1].max().max() <= 1  # test p-values are valid


@pytest.mark.filterwarnings("ignore:Precision loss")
def test_aggregate_ttest_summary(adata_treat):
    agg = sm.pp.aggregate_ttest(adata_treat, treatment_key="TARGETGENE")[0]
    t = sm.pp.tstat_distance(agg)
    assert t.shape == (1,)


def test_scale_by_plate(adata):
    sm.pp.scale_by_batch(adata, batch_key="Image_Metadata_Plate")
    X = pd.concat([adata.obs["Image_Metadata_Plate"], adata[:, 0].to_df()], axis=1)
    assert all(X.groupby("Image_Metadata_Plate").mean() < 1e-7)
