import pytest

import scmorph as sm


@pytest.fixture
def adata_treat():
    adata = sm.datasets.rohban2017_minimal()
    adata.obs["TARGETGENE"] = adata.obs["TARGETGENE"].astype(str).replace("nan", "UNTREATED")
    sm.pp.drop_na(adata)
    sm.pp.scale_by_batch(
        adata,
        batch_key="Image_Metadata_Plate",
        treatment_key="TARGETGENE",
        control="UNTREATED",
    )
    return adata


def test_get_ks(adata_treat):
    adata = adata_treat
    ref_ks, treat_ks = sm.tl.get_ks(
        adata,
        treatment_key="TARGETGENE",
        control="UNTREATED",
        well_key="Image_Metadata_Well",
        batch_key="Image_Metadata_Plate",
        control_wells=None,
    )
    assert ref_ks.shape == (10, 7)
    assert treat_ks.shape == (5, 8)
    assert treat_ks["is_significant_0.05"].all()
    assert not ref_ks["is_significant_0.05"].any()
