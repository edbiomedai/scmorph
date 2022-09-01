import scmorph as sm

filename = "data/test_data.csv"
n_headers = 1
data_nrows = 1659


def test_pca():
    adata = sm.read_cellprofiler(filename, n_headers=n_headers)
    sm.pp.scale(adata)
    sm.pp.pca(adata, n_comps=2)
    assert adata.obsm["X_pca"].shape == (data_nrows, 2)
