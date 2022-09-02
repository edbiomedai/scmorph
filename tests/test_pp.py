import scmorph as sm

data_nrows = 1659


def test_pca():
    adata = sm.datasets.multi_plate_experiment()
    sm.pp.scale(adata)
    sm.pp.pca(adata, n_comps=2)
    assert adata.obsm["X_pca"].shape == (data_nrows, 2)
