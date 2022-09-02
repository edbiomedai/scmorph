import functools
import glob
import os
import re
import sys
from typing import Any, List, Tuple, Union

import numpy as np
import pandas as pd
import pyarrow
from anndata import AnnData
from scanpy import read_h5ad

from scmorph._logging import get_logger

if sys.version_info >= (3, 9):
    Pattern = re.Pattern
elif sys.version_info >= (3, 8):
    from typing import Pattern
else:
    from typing.re import Pattern

"""
Functions to generate AnnData objects from a morphological datasets in .csv format

Key components:
- read_cellprofiler : convert .csv files produced by CellProfiler into AnnData objects
- _make_AnnData : lower-level wrapper to create AnnData objects from pd.DataFrames stored in memory

Note that some of these functions consume a lot of memory when reading in large .csv files.
For the creation of file-backed AnnData files it is thus advisable to use environments
with a lot of available RAM. Alternatively, `read_cellprofiler_from_path` automatically
discovers csv files and iteratively creates an AnnData object from the files. This is
memory-friendly but relies on you having one file per well.
"""

log = get_logger()


def _parse_csv_headers(
    filename: Union[str, List["str"]],
    n_headers: int = 1,
    sanitize: bool = True,
    sep: str = ",",
) -> List["str"]:
    """
    Parses csv file headers with multiple headers.

    Parameters
    ----------
    filename : str
            Path to .csv file. If list is given, will return header of first file.

    n_headers : int
            1-indexed row number of last header

    sanitize : bool
            Remove everything after last dot of headers?

    sep : char
            Column deliminator

    Returns
    ----------
    list with merged header names and length equal to number
    of columns in csv file
    """
    if isinstance(filename, list):
        filename = filename[0]

    df = pd.read_csv(filename, header=None, nrows=n_headers, sep=sep)

    if sanitize:
        reg = r"(.*)\.[^.]*$"
        df.replace(reg, r"\1", regex=True, inplace=True)

    out = df.agg("_".join, axis=0).tolist()
    return out


def _parse_csv(
    path: Union[str, List["str"]], n_headers: int = 1, sep: str = ","
) -> pd.DataFrame:
    """
    Parses csv files with multiple headers.

    Parameters
    ----------
    path : str
            Path to .csv file. If list is given, will append files vertically and use
            header of first file

    n_headers : int
            1-indexed row number of last header

    sep : char
            Column deliminator

    Returns
    ----------
    pd.Dataframe

    Note
    ----------
    Depending on the size of the input matrix, this function can take a lot of memory
    """

    # get header information
    head = _parse_csv_headers(path, n_headers, sep=sep)

    _read_csv = functools.partial(
        pd.read_csv,
        sep=sep,
        names=head,
        skiprows=n_headers,
        header=None,
        engine="pyarrow",
    )

    path_is_list = isinstance(path, list)

    # read in data including header
    if not path_is_list:
        df = _read_csv(path)
    else:
        df = pd.concat([_read_csv(f) for f in path], axis=0)

    return df


def _split_feature_names(
    features: Union[pd.Series, List["str"]], feature_delim: str = "_"
) -> pd.DataFrame:
    """
    Split feature names into pd.DataFrame

    Parameters
    ----------
    features : pd.Series or list
            Feature names

    feature_delim : str
            Character delimiting feature names

    Returns
    ----------
    pd.DataFrame of feature names split into columns
    """

    features = pd.Series(features)

    df = features.str.split(feature_delim, expand=True)  # split feature names
    df.index = features
    df.columns = ["feature_" + str(i) for i in df.columns]  # ensure str columns
    return df


def _split_meta(
    df: pd.DataFrame, meta_cols: Union[List["str"], None]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split pd.DataFrame into two, keeping metadata and measurements separate

    Parameters
    ----------
    df : pd.DataFrame
            DataFrame with metadata and measurement columns

    meta_cols : list, optional
            Names of metadata columns. None for automatic detection. Default: None

    Returns
    ----------
    Tuple of pd.DataFrames where first element represents metadata and second element is measurements.

    Note
    ----------
    Note that `df` is modified in-place to stay memory efficient
    """
    meta_cols = _match_meta(df.columns, meta_cols)

    meta = df.loc[:, meta_cols]
    df.drop(columns=meta_cols, inplace=True)

    return (meta, df)


def _make_AnnData(
    df: pd.DataFrame,
    meta_cols: Union[List["str"], None] = None,
    feature_delim: str = "_",
) -> AnnData:
    """
    Make annotated data matrix from pd.DataFrame

    Parameters
    ----------
    df : pd.DataFrame
            Phenotypic measurements, e.g. derived from CellProfiler

    meta_cols : list, optional
        Names of metadata columns. None for automatic detection. Default: None

    feature_delim : str
            Character delimiting feature names

    Returns
    ----------
    Annotated data matrix
    """

    meta, X = _split_meta(df, meta_cols=meta_cols)
    dropcols = _match_drop(df.columns)
    if dropcols:
        log.warning(
            "Non-continous and rotation-variant features are not currently supported and "
            + "will be discarded! The following features are dropped:\n"
            + "\n".join(dropcols)
        )
        df.drop(columns=dropcols, inplace=True)

    featData = _split_feature_names(X.columns, feature_delim=feature_delim)

    ad = AnnData(
        X=X.to_numpy(),
        obs=meta.to_dict("list"),  # avoid conversion bugs
        var=featData,
        dtype="float32",  # avoid future warning
    )

    return ad


def _find_files(path: Union[str, List["str"]], suffix: str = ".csv") -> List["str"]:
    """
    Find single-cell csv files recursively

    Parameters
    ----------
    path : str
            Path to input directory. If a path to a matching file is given, will return that path.

    suffix : str
            File suffix to match. Default: .csv
    Returns
    -------
    List of str
        Matching files
    """
    # check input modes
    if isinstance(path, str):
        if path.endswith("f{suffix}"):
            return [path]
        elif os.path.isdir(path):
            path = [path]
        else:
            raise ValueError(f"{path} is neither a {suffix} nor a directory")

    path = [os.path.abspath(p) for p in path]

    # recursively find csv files in path
    files = [glob.glob(f"{p}/**/*{suffix}", recursive=True) for p in path]
    return np.hstack(files).tolist()


def read_cellprofiler(
    filename: str,
    n_headers: int = 1,
    meta_cols: Union[List["str"], None] = None,
    feature_delim: str = "_",
    sep: str = ",",
) -> AnnData:
    """
    Read a matrix from a .csv file created with CellProfiler

    Parameters
    ----------
    filename : str
            Path to .csv file

    n_headers : int, optional
            Number of header rows. Default: 2

    meta_cols: list, optional
            Names of metadata columns. None for automatic detection. Default: None

    feature_delim : str, optional
            Feature deliminator. Default: "_"

    sep : str, optional
            Column deliminator. Default: ","

    Returns
    ----------
    adata : :class:`~AnnData`

    Note
    ----------
    Depending on the size of the input matrix, this function can take a lot of memory.
    If needed, try exporting CellProfiler in batches of smaller csv files and read them in using
    :func:`read_cellprofiler_batch`.
    """
    # TODO: think about having temporary file-backing to lower memory usage
    df = _parse_csv(filename, n_headers, sep=sep)

    # convert df to AnnData objects
    ad = _make_AnnData(df, meta_cols=meta_cols, feature_delim=feature_delim)

    return ad


def read_cellprofiler_batches(
    path: str,
    output_file: str,
    file_pattern: str = "Nuclei.csv",
    n_headers: int = 1,
    meta_cols: Union[List["str"], None] = None,
    sep: str = ",",
) -> AnnData:
    """
    Read CellProfiler data from directories

    Parameters
    ----------
    path : str
            Path to a directory containing .csv files

    output_file : str
            Path to output file, will create a .h5ad file. This is needed
            to prevent large memory allocations.

    file_pattern : str
            Pattern to match .csv files. Default: "Nuclei.csv"

    n_headers : int, optional
            Number of header rows. Default: 2

    meta_cols: list, optional
            Names of metadata columns. None for automatic detection. Default: None

    feature_delim : str, optional
            Feature deliminator. Default: "_"

    sep : str, optional
            Column deliminator. Default: ","

    progress : bool, optional
            Show progress bar. Default: True

    Returns
    ----------
    adata : :class:`~AnnData`
    """
    import anndata as ad
    import h5py
    from anndata.experimental import write_elem
    from tqdm import tqdm

    tqdm = functools.partial(tqdm, unit=" files", dynamic_ncols=True, mininterval=1)

    files = _find_files(path, suffix=file_pattern)

    if len(files) == 0:
        raise ValueError(f"No files ending in {file_pattern} found in {path}")

    sample_file = files[0]

    log.info(f"Found {len(files)} files")
    log.info("Reading in all metadata...")

    # read in obs metadata
    obs = pd.concat(
        [
            read_meta(f, n_headers=n_headers, meta_cols=meta_cols, sep=sep)
            for f in tqdm(files)
        ]
    ).reset_index(drop=True)
    obs.fillna("", inplace=True)
    obs.index = obs.index.astype(str)

    # extract var metadata from first file
    var = read_cellprofiler(sample_file, sep=sep, n_headers=n_headers).var
    var.fillna("", inplace=True)
    var.index = var.index.astype(str)

    # create output file with metadata and empty X
    log.info("Creating intermediary output file...")
    with h5py.File(output_file, "w") as target:
        target.create_dataset(
            "X",
            (obs.shape[0], var.shape[0]),
            dtype="float32",
            chunks=(min(10000, obs.shape[0]), min(10, var.shape[0])),
        )
        write_elem(target, "obs", obs)
        write_elem(target, "var", var)

    # read in created output file
    adata = ad.read(output_file, backed="r+")

    log.info("Converting all data. This may take a while...")
    counter = 0
    for f in tqdm(files):
        cur_X = read_X(f, meta_cols=meta_cols, n_headers=n_headers, sep=sep)
        adata[counter : counter + cur_X.shape[0], :].X = cur_X
        counter += cur_X.shape[0]
    return adata


def _meta_terms() -> Pattern[str]:
    filters = [
        "^Image_",
        "FileName",
        "Object.?Number",
        "Image.?Number",
        "^URL",
        "Metadata",
    ]
    filt = "|".join(filters)
    return re.compile(filt)


def _drop_terms() -> Pattern[str]:
    filters = [
        "Phase",
        "NumberOfNeighbors",
        "Orientation",
        "Extent",
        "BoundingBox",
        "SpatialMoment",
        "CentralMoment",
        "NormalizedMoment",
        "InertiaMoment",  # TODO: double-check whether to kick this out
        "Location",
        "[XY]$",
    ]
    filt = "|".join(filters)
    return re.compile(filt)


def _match_meta(
    header: List[str], meta_cols: Union[None, List[str]] = None
) -> List[str]:
    import re

    re_meta = _meta_terms()
    meta_cols = [col for col in header if re.search(re_meta, col)]
    return meta_cols


def _match_drop(
    header: List[str], meta_cols: Union[None, List[str]] = None
) -> List[str]:
    import re

    re_drop = _drop_terms()
    drop_cols = [col for col in header if re.search(re_drop, col)]
    return drop_cols


def _read_csv_columns(
    path: str,
    columns: List[str],
    column_names: List[str],
    n_headers: int = 1,
    sep: str = ",",
) -> pyarrow.Table:
    """
    Read specific columns from a .csv file given then column names

    Parameters
    ----------
    path : str
        Path to csv file
    columns : List[str]
        Names of columns to include
    column_names : List[str]
        Column names
    n_headers : int, optional
        Number of headers, by default 1
    sep : str, optional
        Column deliminiator, by default ","

    Returns
    -------
    Pyarrow table
    """
    from pyarrow import csv

    parseopts = csv.ParseOptions(delimiter=sep)
    readopts = csv.ReadOptions(skip_rows=n_headers, column_names=column_names)
    convopts = csv.ConvertOptions(include_columns=columns)
    tab = csv.read_csv(path, readopts, parseopts, convopts)
    return tab


def read_meta(
    path: str,
    meta_cols: Union[List[str], None] = None,
    n_headers: int = 1,
    sep: str = ",",
) -> pd.DataFrame:
    """
    Read metadata from a .csv file

    Parameters
    ----------
    path : str
            Path to .csv file

    meta_cols: list, optional
            Names of metadata columns. None for automatic detection. Default: None

    Returns
    ----------
    meta : :class:`~pandas.DataFrame`
    """
    header = _parse_csv_headers(path, n_headers=n_headers, sep=sep)
    meta_cols = _match_meta(header, meta_cols)
    df = _read_csv_columns(
        path=path, columns=meta_cols, column_names=header, sep=sep, n_headers=n_headers
    )
    df = df.to_pandas()
    return df


def read_X(
    path: str,
    meta_cols: Union[List[str], None] = None,
    n_headers: int = 1,
    sep: str = ",",
) -> np.array:
    """
    Read X from a .csv file

    Parameters
    ----------
    path : str
            Path to .csv file

    meta_cols: list, optional
            Names of metadata columns. None for automatic detection. Default: None

    Returns
    ----------
    X : :class:`~numpy.array`
    """
    header = _parse_csv_headers(path, n_headers=n_headers, sep=sep)
    meta_cols = _match_meta(header, meta_cols)
    drop_cols = _match_drop(header, meta_cols)
    columns = [col for col in header if col not in meta_cols + drop_cols]
    tab = _read_csv_columns(
        path=path, columns=columns, column_names=header, sep=sep, n_headers=n_headers
    )
    arr = np.array(tab, dtype="float32").T
    return arr


def read(filename: str, **kwargs: Any) -> AnnData:
    """
    Read csv or h5ad files.

    This function wraps read_cellprofiler and read_h5ad and uses to appropriate one depending on file ending.

    Parameters
    ----------
    filename : str
            Path to .csv or h5ad file

    kwargs : Any, optional
            Other parameters passed to `read_cellprofiler` or `read_h5ad`

    Returns
    -------
    adata : :class:`~AnnData`
    """
    _, fileending = os.path.splitext(filename)
    if fileending == ".csv":
        return read_cellprofiler(filename, **kwargs)
    elif fileending == ".h5ad":
        return read_h5ad(filename, **kwargs)
    else:
        raise ValueError(f"File ending {fileending} not supported")
