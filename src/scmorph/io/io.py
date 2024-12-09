import functools
import glob
import os
import re
from typing import Any

import numpy as np
import pandas as pd
import pyarrow
from anndata import AnnData
from scanpy import read_h5ad

from scmorph.logging import get_logger

log = get_logger()


def _parse_csv_headers(
    filename: str | list[str],
    n_headers: int = 1,
    sanitize: bool = True,
    sep: str = ",",
) -> list[str]:
    """
    Parses csv file headers with multiple headers.

    Parameters
    ----------
    filename
        Path to .csv file. If list is given, will return header of first file.

    n_headers
        1-indexed row number of last header

    sanitize
        Remove everything after last dot of headers?

    sep
        Column deliminator

    Returns
    -------
    list with merged header names and length equal to number
    of columns in csv file
    """
    if isinstance(filename, list):
        filename = filename[0]

    df = pd.read_csv(filename, header=None, nrows=n_headers, sep=sep)

    if sanitize:
        reg = r"(.*)\.[^.]*$"
        df.replace(reg, r"\1", regex=True, inplace=True)

    return df.agg("_".join, axis=0).tolist()


def _parse_csv(
    path: str | list[str],
    n_headers: int = 1,
    sep: str = ",",
    backup_url: str | None = None,
) -> pd.DataFrame:
    """
    Parses csv files with multiple headers.

    Parameters
    ----------
    path
        Path to .csv file. If list is given, will append files vertically and use
        header of first file

    n_headers
        1-indexed row number of last header

    sep
        Column deliminator

    backup_url
        Retrieve the file from an URL if not present on disk.

    Note
    ----------
    Depending on the size of the input matrix, this function can take a lot of memory
    """
    if isinstance(path, list):
        for p in path:
            _cache_file(p, backup_url=backup_url)
    else:
        _cache_file(path, backup_url=backup_url)

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


def split_feature_names(features: pd.Series | list[str], feature_delim: str = "_") -> pd.DataFrame:
    """
    Split feature names into a :class:`~pandas.DataFrame`

    Parameters
    ----------
    features
        Feature names

    feature_delim
        Character delimiting feature names

    Returns
    -------
    Feature names split into columns of a DataFrame
    """
    features = pd.Series(features)

    df = features.str.split(feature_delim, expand=True)  # split feature names
    df.index = features
    df.columns = [f"feature_{str(i)}" for i in df.columns]
    return df


def _split_meta(df: pd.DataFrame, meta_cols: list[str] | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split DataFrame into two, keeping metadata and measurements separate

    Parameters
    ----------
    df
        DataFrame with metadata and measurement columns

    meta_cols
        Names of metadata columns. `None` for automatic detection.

    Returns
    -------
    Tuple of DataFrames where first element represents metadata and second element is measurements.

    Note
    ----------
    Note that `df` is modified in-place to stay memory efficient
    """
    meta_cols = _match_meta(df.columns, meta_cols)

    meta = df.loc[:, meta_cols]
    df.drop(columns=meta_cols, inplace=True)

    return (meta, df)


def make_AnnData(
    df: pd.DataFrame,
    meta_cols: list[str] | None = None,
    feature_delim: str = "_",
) -> AnnData:
    """
    Make annotated data matrix from :class:`~pandas.DataFrame`

    Parameters
    ----------
    df
        Phenotypic measurements, e.g. derived from CellProfiler

    meta_cols
        Names of metadata columns. `None` for automatic detection.

    feature_delim
            Character delimiting feature names

    Returns
    -------
    Annotated data matrix
    """
    meta, X = _split_meta(df, meta_cols=meta_cols)
    dropcols = _match_drop(df.columns)
    if dropcols:
        log.warning(
            "Non-continous and rotation-variant features are not currently supported and "
            "will be discarded! The following features are dropped:\n%s",
            "\n".join(dropcols),
        )
        df.drop(columns=dropcols, inplace=True)

    featData = split_feature_names(X.columns, feature_delim=feature_delim)

    ad = AnnData(
        X=X.to_numpy(),
        obs=meta.to_dict("list"),  # avoid conversion bugs
        var=featData,
    )

    return ad


def _find_files(path: str | list[str], suffix: str = ".csv") -> list[str]:
    """
    Find single-cell csv files recursively

    Parameters
    ----------
    path
        Path to input directory. If a path to a matching file is given, will return that path.

    suffix
        File suffix to match.

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


def read_cellprofiler_csv(
    filename: str,
    n_headers: int = 1,
    meta_cols: list[str] | None = None,
    feature_delim: str = "_",
    sep: str = ",",
    backup_url: str | None = None,
) -> AnnData:
    """
    Read a matrix from a .csv file created with CellProfiler

    Note
    -----
    Depending on the CellProfiler version, you may have one or two header rows.
    Before using this function, do check the csv file to see if feature names
    are contained only in the first row or split over two and set `n_headers` accordingly.

    Parameters
    ----------
    filename
        Path to .csv file

    n_headers
        Number of header rows.

    meta_cols: list
        Names of metadata columns. `None` for automatic detection.

    feature_delim
        Feature deliminator.

    backup_url
        Retrieve the file from an URL if not present on disk.

    sep
        Column deliminator.

    Returns
    -------
    The read in AnnData object

    Note
    ----------
    Depending on the size of the input matrix, this function can take a lot of memory.
    If needed, try exporting CellProfiler in batches of smaller csv files and read them in using
    :func:`.io.read_cellprofiler_batches`.
    """
    # TODO: think about having temporary file-backing to lower memory usage
    df = _parse_csv(filename, n_headers, sep=sep, backup_url=backup_url)

    return make_AnnData(df, meta_cols=meta_cols, feature_delim=feature_delim)


def read_cellprofiler_batches(
    path: str,
    output_file: str,
    file_pattern: str = "Nuclei.csv",
    n_headers: int = 1,
    meta_cols: list[str] | None = None,
    sep: str = ",",
) -> AnnData:
    """
    Read CellProfiler data from directories

    Parameters
    ----------
    path
        Path to a directory containing .csv files

    output_file
        Path to output file, will create a .h5ad file. This is needed
        to prevent large memory allocations.

    file_pattern
        re.Pattern to match .csv files.

    n_headers
        Number of header rows.

    meta_cols
        Names of metadata columns. `None` for automatic detection.

    feature_delim
        Feature deliminator.

    sep
        Column deliminator.

    progress
        Show progress bar.

    Returns
    -------
    adata
    """
    import anndata as ad
    import h5py
    from tqdm import tqdm

    try:
        from anndata.io import write_elem
    except ImportError:
        from anndata.experimental import write_elem

    tqdm = functools.partial(tqdm, unit=" files", dynamic_ncols=True, mininterval=1)

    files = _find_files(path, suffix=file_pattern)

    if len(files) == 0:
        raise ValueError(f"No files ending in {file_pattern} found in {path}")

    sample_file = files[0]

    log.info("Found %s files", len(files))
    log.info("Reading in all metadata...")

    # read in obs metadata
    obs = pd.concat([read_meta(f, n_headers=n_headers, meta_cols=meta_cols, sep=sep) for f in tqdm(files)]).reset_index(
        drop=True
    )
    obs.fillna("", inplace=True)
    obs.index = obs.index.astype(str)

    # extract var metadata from first file
    var = read_cellprofiler_csv(sample_file, sep=sep, n_headers=n_headers).var
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


def _meta_terms() -> re.Pattern[str]:
    filters = [
        "^Image_",
        "FileName",
        "Object.?Number",
        "Image.?Number",
        "^URL",
        "Metadata",
        "Parent_(Cells|Nuclei)",
        "Children",
        "TableNumber",
    ]
    filt = "|".join(filters)
    return re.compile(filt)


def _drop_terms() -> re.Pattern[str]:
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
        "[XYZ]$",
    ]
    filt = "|".join(filters)
    return re.compile(filt)


def _match_meta(header: list[str], meta_cols: None | list[str] = None) -> list[str]:
    re_meta = _meta_terms()
    meta_cols = [col for col in header if re.search(re_meta, col)]
    return meta_cols


def _match_drop(header: list[str]) -> list[str]:
    re_drop = _drop_terms()
    return [col for col in header if re.search(re_drop, col)]


def _read_csv_columns(
    path: str,
    columns: list[str],
    column_names: list[str],
    n_headers: int = 1,
    sep: str = ",",
) -> pyarrow.Table:
    """
    Read specific columns from a .csv file given then column names

    Parameters
    ----------
    path
        Path to csv file
    columns
        Names of columns to include
    column_names
        Column names
    n_headers
        Number of headers
    sep
        Column deliminiator

    Returns
    -------
    Table with the requested columns
    """
    from pyarrow import csv

    parseopts = csv.ParseOptions(delimiter=sep)
    readopts = csv.ReadOptions(skip_rows=n_headers, column_names=column_names)
    convopts = csv.ConvertOptions(include_columns=columns)
    return csv.read_csv(path, readopts, parseopts, convopts)


def _cache_file(path: str, backup_url: str | None) -> None:
    """
    Check if a file is present, download if not.

    Parameters
    ----------
    path
        Path to file
    backup_url
        URL to backup file

    Note
    ----------
    Uses scanpy's functionality for file caching.
    All rights for these functions lie with the license holders:
    https://github.com/scverse/scanpy/blob/master/LICENSE
    """
    from scanpy.readwrite import _check_datafile_present_and_download

    log.debug("Checking if file is present, downloading if not...")
    is_present = _check_datafile_present_and_download(path, backup_url=backup_url)

    if not is_present:
        log.error("%s not found, and could not be downloaded.", path)


def read_meta(
    path: str,
    meta_cols: list[str] | None = None,
    n_headers: int = 1,
    sep: str = ",",
) -> pd.DataFrame:
    """
    Read metadata from a .csv file

    Parameters
    ----------
    path
        Path to .csv file

    meta_cols: list
        Names of metadata columns. `None` for automatic detection.

    Returns
    -------
    The read in metadata table.
    """
    header = _parse_csv_headers(path, n_headers=n_headers, sep=sep)
    meta_cols = _match_meta(header, meta_cols)
    df = _read_csv_columns(path=path, columns=meta_cols, column_names=header, sep=sep, n_headers=n_headers)
    df = df.to_pandas()
    return df


def read_X(
    path: str,
    meta_cols: list[str] | None = None,
    n_headers: int = 1,
    sep: str = ",",
) -> np.ndarray:
    """
    Read X from a .csv file

    Parameters
    ----------
    path
        Path to .csv file

    meta_cols: list
        Names of metadata columns. `None` for automatic detection.

    Returns
    -------
    The read in feature matrix.
    """
    header = _parse_csv_headers(path, n_headers=n_headers, sep=sep)
    meta_cols = _match_meta(header, meta_cols)
    drop_cols = _match_drop(header)
    columns = [col for col in header if col not in meta_cols + drop_cols]
    tab = _read_csv_columns(path=path, columns=columns, column_names=header, sep=sep, n_headers=n_headers)
    return np.array(tab, dtype="float32").T


def read_sql(filename: str, backup_url: str | None = None) -> AnnData:
    """
    Read sql files.

    Parameters
    ----------
    filename
        Path to .sql file

    backup_url
        URL to backup file.

    Returns
    -------
    adata
    """
    import sqlite3

    _cache_file(filename, backup_url=backup_url)

    known_tables = ["Image", "Nuclei", "Cytoplasm", "Cells"]
    conn = sqlite3.connect(filename)
    c = conn.cursor()
    c.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = [i[0] for i in c.fetchall()]

    # sanity check
    if unknown_tables := list(set(tables) - set(known_tables)):
        log.warning("Unknown tables found in SQL database: %s", unknown_tables)

    keep_tables = list(set(tables) & set(known_tables))

    data = {  # nosec
        i: pd.read_sql_query(f"SELECT * from {i}", conn) for i in keep_tables
    }

    meta_present = False
    if "Image" in data:
        meta = data.pop("Image")
        meta_present = True

    X = pd.concat(data, axis=1, ignore_index=False)
    # remove table names from index, this information is in features
    if isinstance(X.columns, pd.MultiIndex):
        X.columns = X.columns.droplevel(0)
    adata = make_AnnData(X)

    if meta_present:
        # keep only subset of metadata table, to avoid creating large memory and disk-size overhead
        # for information that is likely not needed. If you do need other information contained
        # in the Image table, please open a GitHub issue.
        log.info("Adding metadata to AnnData object...")

        meta_regex_keep = r"^Metadata|TableNumber|Count_Cells|Count_Cytoplasm|Count_Nuclei"
        meta_cols_keep = [i for i in meta.columns if re.match(meta_regex_keep, i)]
        log.info(
            "Metadata found in SQL database, adding to AnnData object. Will only keep the following columns: %s",
            ", ".join(meta_cols_keep),
        )

        adata.obs = pd.merge(adata.obs, meta[meta_cols_keep])

    return adata


def read(filename: str, **kwargs: Any) -> AnnData:
    """
    Read csv, h5ad or sql files.

    This function wraps :func:`~.read_cellprofiler_csv`, :func:`~.read_sql`
    and :func:`~scanpy.read_h5ad` and uses to appropriate one depending on file ending.
    For details, see the respective functions.

    Parameters
    ----------
    filename
        Path to .csv or h5ad file

    kwargs
        Other parameters passed to the corresponding functions.

    Returns
    -------
    The read in AnnData object, possibly after conversion.
    """
    _, fileending = os.path.splitext(filename)
    if fileending == ".csv":
        return read_cellprofiler_csv(filename, **kwargs)
    elif fileending == ".h5ad":
        return read_h5ad(filename, **kwargs)
    elif fileending in [".sql", ".sqlite"]:
        return read_sql(filename, **kwargs)
    else:
        raise ValueError(f"File ending {fileending} not supported")
