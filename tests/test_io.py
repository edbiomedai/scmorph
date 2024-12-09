import tempfile

import pandas as pd
import pytest
from anndata import AnnData

import scmorph as sm

tmpfile = tempfile.NamedTemporaryFile("w+b", suffix=".h5ad").name

data_nrows = 12352
raw_cols = 1703
n_headers = 2
meta_cols, feature_cols = 16, 1687


@pytest.fixture
def rohban():
    return sm.datasets.rohban2017_minimal()


@pytest.fixture
def rohban_minimal_csv_file():
    return sm.datasets._datasets.rohban2017_minimal_csv()


def test_parse_csv_header(rohban_minimal_csv_file):
    header = sm.io.io._parse_csv_headers(rohban_minimal_csv_file, n_headers=n_headers, sanitize=True, sep=",")
    assert isinstance(header, list) and len(header) == raw_cols


def test_parse_csv(rohban_minimal_csv_file):
    df = sm.io.io._parse_csv(rohban_minimal_csv_file, n_headers=n_headers, sep=",")
    assert isinstance(df, pd.DataFrame) and df.shape == (4, raw_cols)


def test_split_feature_names(rohban):
    features = sm.io.io.split_feature_names(rohban.var.index, feature_delim="_")
    assert rohban.var.shape == features.shape
    assert isinstance(features, pd.DataFrame)
    assert features.shape == (feature_cols, 6)


def test_split_meta(rohban_minimal_csv_file):
    df = sm.io.io._parse_csv(rohban_minimal_csv_file, n_headers=n_headers)
    meta, X = sm.io.io._split_meta(df, meta_cols=None)
    assert isinstance(meta, pd.DataFrame) and meta.shape == (4, meta_cols)


def test_make_AnnData(rohban):
    df = pd.concat([rohban.obs, rohban.to_df()], axis=1)
    adata = sm.io.make_AnnData(df, feature_delim="_")
    assert isinstance(adata, AnnData)


def test_read_cellprofiler(rohban_minimal_csv_file):
    adata = sm.read_cellprofiler_csv(rohban_minimal_csv_file, n_headers=n_headers)
    assert isinstance(adata, AnnData)
    assert adata.shape == (4, feature_cols)


def test_read_h5ad(rohban):
    rohban.write(tmpfile)
    rohban.file.close()
    adata = sm.read_h5ad(tmpfile, backed="r+")
    assert isinstance(adata, AnnData) and adata.shape == (data_nrows, feature_cols)
    adata.file.close()
