import numpy as np
import pandas as pd
import pytest

import scmorph as sm

data_nrows_no_na = 12286
data_nrows = 12352


@pytest.fixture
def adata():
    return sm.datasets.rohban2017_minimal()


@pytest.fixture
def adata_treat(adata_no_na):
    adata = adata_no_na
    adata.obs["TARGETGENE"] = adata.obs["TARGETGENE"].astype(str).replace("nan", "DMSO")
    return adata


@pytest.fixture
def adata_treat_mini(adata_treat):
    return adata_treat[:200, :10]


def test_drop_na(adata_no_na):
    assert adata_no_na.shape == (data_nrows_no_na, 1687)


def test_pca(pca_result):
    assert pca_result.obsm["X_pca"].shape == (data_nrows_no_na, 2)


def test_aggregate_modes(adata):
    modes = ["mean", "median", "std", "var", "sem", "mad", "mad_scaled"]
    for m in modes:
        agg = sm.pp.aggregate(
            adata,
            method=m,
            group_keys=["Image_Metadata_Well", "Image_Metadata_Plate"],
            progress=False,
        )
        assert agg.shape == (20, adata.shape[1])


def test_aggregate_mahalanobis(adata_treat):
    agg = sm.pp.aggregate_mahalanobis(
        adata_treat,
        treatment_key="TARGETGENE",
        control="DMSO",
        group_key="Image_Metadata_Well",
    )
    assert agg.shape == (1,)


def test_single_cell_mahalanobis(adata_fixed_values):
    agg = sm.pp.aggregate_mahalanobis(
        adata_fixed_values,
        treatment_key="Treatment",
        control="DMSO",
        cov_from_single_cell=True,
    )
    assert "treated" in agg.index
    assert np.allclose(agg["treated"], 1.391, atol=1e-3)
    assert agg.shape == (1,)


def test_aggregate_pc(adata_treat):
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


def test_scale_by_plate_with_control(adata_treat):
    sm.pp.scale_by_batch(
        adata_treat,
        batch_key="Image_Metadata_Plate",
        treatment_key="TARGETGENE",
        control="DMSO",
    )
    X = pd.concat(
        [
            adata_treat.obs[["Image_Metadata_Plate", "TARGETGENE"]],
            adata_treat[:, 0].to_df(),
        ],
        axis=1,
    )
    ref_vals = (
        X.groupby(["Image_Metadata_Plate", "TARGETGENE"])
        .mean()
        .reset_index()
        .query("TARGETGENE == 'DMSO'")["Nuclei_AreaShape_Area"]
        .values
    )
    assert np.allclose(ref_vals, 0, atol=1e-7)


def test_batch_effect_removal(adata_treat):
    # test with copy

    adata_corr = sm.pp.remove_batch_effects(
        adata_treat,
        batch_key="Image_Metadata_Plate",
        treatment_key="TARGETGENE",
        control="DMSO",
        copy=True,
    )

    # test without copy
    sm.pp.remove_batch_effects(
        adata_treat,
        batch_key="Image_Metadata_Plate",
        treatment_key="TARGETGENE",
        control="DMSO",
    )

    assert np.allclose(adata_corr[:, 0].X, adata_treat[:, 0].X)

    agg = sm.pp.aggregate(adata_treat, ["Image_Metadata_Plate", "TARGETGENE"], "mean")
    feat_values = agg[agg.obs["TARGETGENE"] == "DMSO", 0].X.flatten()
    assert np.allclose(feat_values - np.mean(feat_values), 0, atol=0.001)


def test_bio_batch_effect_removal(adata_treat):
    adata_treat.obs["cell_line"] = "A cell line"
    adata_treat.obs.loc[:"199", "cell_line"] = "Another cell line"

    sm.pp.remove_batch_effects(
        adata_treat,
        batch_key="Image_Metadata_Plate",
        bio_key="cell_line",
        treatment_key="TARGETGENE",
        control="DMSO",
    )

    assert "A cell line" in adata_treat.uns["batch_effects"].columns
    assert "Another cell line" in adata_treat.uns["batch_effects"].columns

    # Check for bio effect retention
    agg = (
        pd.merge(
            adata_treat.obs["cell_line"],
            adata_treat[adata_treat.obs["TARGETGENE"] == "DMSO", "Nuclei_AreaShape_Area"].to_df(),
            left_index=True,
            right_index=True,
        )
        .groupby("cell_line")
        .mean()
        .diff()
    )

    measured_effect = adata_treat.uns["batch_effects"].loc[
        "Nuclei_AreaShape_Area", "Another cell line"
    ]
    recovered_effect = agg.loc["Another cell line"]

    assert np.allclose(measured_effect, recovered_effect, atol=0.0001)


def test_chatterjee_correlation():
    np.random.seed(2025)
    x1 = np.random.normal(size=100)
    x2 = x1 * 10 + np.random.normal(scale=0.5, size=100)
    x3 = np.random.normal(size=2)
    x2d = np.random.normal(size=(100, 2))

    corr = sm.pp.correlation.xim(x1, x2)
    assert np.all(corr > 0.9)

    corr = sm.pp.correlation.xim(x2d)
    assert corr[0, 1] < 0.1

    with pytest.raises(ValueError):
        sm.pp.correlation.xim(x1, x2d)

    with pytest.raises(ValueError):
        sm.pp.correlation.xim(x1)

    with pytest.raises(ValueError):
        sm.pp.correlation.xim(x1, x3)


@pytest.mark.filterwarnings("ignore:RuntimeWarning")
def test_correlation(adata_treat_mini):
    corr = sm.pp.correlation.corr(adata_treat_mini.X, method="pearson")
    assert corr.shape == (adata_treat_mini.shape[1], adata_treat_mini.shape[1])
    assert np.all(corr <= 1) and np.all(corr >= -1)

    corr = sm.pp.correlation.corr(adata_treat_mini.X, method="spearman")
    assert corr.shape == (adata_treat_mini.shape[1], adata_treat_mini.shape[1])
    assert np.all(corr <= 1) and np.all(corr >= -1)


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_feature_correlations(adata):
    raw_shape = adata.shape
    sm.pp.select_features(adata)
    assert adata.shape[1] < raw_shape[1]


def test_kruskal_test(adata):
    sm.pp.kruskal_test(adata, "Image_Metadata_Plate")
    assert "kruskal_test" in adata.uns_keys()
    assert adata.uns["kruskal_test"]["Image_Metadata_Plate"].shape == (
        adata.shape[1],
        3,
    )

    adata_filtered = sm.pp.kruskal_filter(adata, "Image_Metadata_Plate", copy=True)
    assert adata_filtered.shape[1] < adata.shape[1]
