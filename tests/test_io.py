import os
import tempfile

import pandas as pd
from anndata import AnnData

import scmorph as sm

filename = "src/scmorph/datasets/test_data.csv"
n_headers = 1
tmpfile = tempfile.NamedTemporaryFile("w+b", suffix=".h5ad").name

data_nrows = 1659
raw_cols = 1083
meta_cols, feature_cols = 9, 920


def test_parse_csv_header():
    header = sm._io.io._parse_csv_headers(filename, n_headers=n_headers, sanitize=True, sep=",")
    assert isinstance(header, list) and len(header) == raw_cols


def test_parse_csv():
    df = sm._io.io._parse_csv(filename, n_headers=n_headers, sep=",")
    assert isinstance(df, pd.DataFrame) and df.shape == (data_nrows, raw_cols)


def test_split_feature_names():
    df = sm._io.io._parse_csv(filename, n_headers=n_headers)
    features = df.columns
    featData = sm._io.io._split_feature_names(features, feature_delim="_")
    assert isinstance(featData, pd.DataFrame) and featData.shape == (raw_cols, 6)


def test_split_meta():
    df = sm._io.io._parse_csv(filename, n_headers=n_headers)
    meta, X = sm._io.io._split_meta(df, meta_cols=None)
    assert isinstance(meta, pd.DataFrame) and meta.shape == (data_nrows, meta_cols)


def test_make_AnnData():
    df = sm._io.io._parse_csv(filename, n_headers=n_headers)
    adata = sm._io.io._make_AnnData(df, feature_delim="_")
    assert isinstance(adata, AnnData)


def test_read_cellprofiler():
    adata = sm.read_cellprofiler(filename, n_headers=n_headers)
    assert isinstance(adata, AnnData)


def test_writing_h5ad():
    adata = sm.datasets.multi_plate_experiment()
    adata.write(tmpfile)
    assert os.path.isfile(tmpfile)


def test_read_h5ad():
    adata = sm.datasets.multi_plate_experiment()
    adata.write(tmpfile)
    adata.file.close()
    del adata  # just to make sure
    adata = sm.read_h5ad(tmpfile, backed="r+")
    assert isinstance(adata, AnnData) and adata.shape == (data_nrows, feature_cols)
    adata.file.close()
