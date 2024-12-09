import re
from collections.abc import Callable, Iterable, Sequence
from functools import partial
from inspect import signature
from textwrap import dedent
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

from scmorph.logging import get_logger


def _infer_names(target: str, options: Iterable[str]) -> Sequence[str]:
    logger = get_logger()

    if target == "batch":
        reg = re.compile("batch|plate")
    elif target in {"well", "group"}:
        reg = re.compile("well$")
    elif target == "treatment":
        reg = re.compile("treatment", re.IGNORECASE)
    elif target == "site":
        reg = re.compile("site$")
    else:
        raise ValueError("type must be one of 'batch', 'well', 'treatment', 'site'")
    res = [x for x in options if reg.search(x)]
    if len(res) > 1:
        logger.warning(
            f"Found multiple {target} columns, are these duplicates?\n"
            + "Will just use the first one, if that is not desired please specify "
            + f"the correct column name using the {target}_key argument.\n"
            + f"Problematic columns were: {', '.join(res)}"
        )
        res = [res[0]]
    return res


def _grouped_obs_fun(
    adata: AnnData,
    group_key: str | list[str],
    fun: Callable[..., Any],
    layer: str | None = None,
    progress: bool = True,
) -> pd.DataFrame:
    """
    Grouped operations on anndata objects

    Slightly adapted from https://github.com/scverse/scanpy/issues/181#issuecomment-534867254
    All copyright lies with Isaac Virshup.
    """
    from tqdm import tqdm

    def getX(adata: AnnData, layer: None | str) -> np.ndarray:
        return adata.X if layer is None else adata.layers[layer]

    grouped = adata.obs.groupby(group_key, observed=True)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names,
    )
    items = tqdm(grouped.indices.items(), unit=" groups") if progress else grouped.indices.items()

    for group, idx in items:
        X = getX(adata[idx], layer)
        out[group] = np.array(fun(X))

    return out


def grouped_op(
    adata: AnnData,
    group_key: str | list[str],
    operation: str,
    layer: str | None = None,
    progress: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Grouped operations on anndata objects

    Parameters
    ----------
    adata
        AnnData object
    group_key
        Column name in `obs` metadata to group by
    operation
        What operation to perform, one of "mean", "logmean", "median", "std",
        "var", "sem", "mad" and "mad_scaled".
    layer
        Which layer ("X" or "X_pca", for example) to aggregate
    progress
         Whether to show a progress bar

    Returns
    -------
    pd.DataFrame
        Data averaged per group in `group_key`
    """
    if operation == "mean":
        fun = partial(np.mean, axis=0, dtype=np.float64, **kwargs)
    elif operation == "logmean":

        def fun(x):  # type: ignore
            return np.mean(np.log1p(x), axis=0, dtype=np.float64, **kwargs)

    elif operation == "median":
        fun = partial(np.median, axis=0, **kwargs)
    elif operation == "std":
        fun = partial(np.std, axis=0, dtype=np.float64, **kwargs)
    elif operation == "var":
        fun = partial(np.var, axis=0, dtype=np.float64, **kwargs)
    elif operation == "sem":
        from scipy.stats import sem

        fun = partial(sem, axis=0, **kwargs)
    elif operation == "mad":
        from scipy.stats import median_abs_deviation as mad

        fun = partial(mad, axis=0, **kwargs)
    elif operation == "mad_scaled":
        from scipy.stats import median_abs_deviation as mad

        f1 = partial(mad, axis=0, scale="normal", **kwargs)
        f2 = partial(np.median, axis=0, **kwargs)

        # Chung 2008, following https://github.com/cytomining/pycytominer/blob/master/pycytominer/operations/transform.py#L175
        def fun(x: np.array) -> np.array:  # type: ignore
            return f2(x) / (f1(x) + 1e-18)

    else:
        raise ValueError("Operation must be one of 'mean', 'median', 'std', 'var', 'sem', 'mad', 'mad_scaled'")

    return _grouped_obs_fun(adata, group_key, fun=fun, layer=layer, progress=progress)


def group_obs_fun_inplace(
    adata: AnnData,
    group_key: str | list[str],
    fun: Callable[..., Any],
    progress: bool = True,
) -> AnnData:
    """
    Alter adata.X inplace by performing fun in each group

    Parameters
    ----------
    adata
        Annotated data matrix object

    group_key
        obs keys to group by

    fun
        Function that takes array and returns array of equal size.
        The function may either only take an array, or the array and the group key.
        In the latter case, the group key must be the second argument!

    Returns
    -------
    AnnData
        Annotated data matrix object after the operation
    """
    from tqdm import tqdm

    grouped = adata.obs.groupby(group_key, observed=True)

    takes_group = len(signature(fun).parameters) > 1

    items = grouped.indices.items()
    items = tqdm(items, unit=" groups") if progress else items

    for group, idx in items:
        X = adata[idx].X
        adata[idx].X = fun(X, group) if takes_group else fun(X)

    return adata


def _get_group_keys(
    adata: AnnData,
    treatment_key: str | None,
    group_key: str | None | list[str],
) -> tuple[list[str], list[str]]:
    # inferring treatment names if necessary
    if isinstance(treatment_key, str):
        if treatment_key == "infer":
            treatment_col = _infer_names("treatment", adata.obs.columns)
        else:
            treatment_col = [treatment_key]
    elif treatment_key is None:
        treatment_col = []
    else:
        treatment_col = treatment_key

    if isinstance(group_key, str):
        if group_key == "infer":
            group_key = _infer_names("group", adata.obs.columns)  # type: ignore
        else:
            group_key = [group_key]
    elif group_key is None:
        group_key = []
    # end inference

    group_keys = [*treatment_col, *group_key]
    group_keys = [x for x in group_keys if x]  # remove None's
    return group_keys, treatment_col  # type: ignore


def get_grouped_op(
    adata: AnnData,
    group_key: list[str],
    operation: str,
    as_anndata: bool = False,
    layer: str | None = None,
    store: bool = True,
    progress: bool = True,
) -> pd.DataFrame | AnnData:
    """
    Retrieve from cache or compute a grouped operation

    Parameters
    ----------
    adata
        AnnData object
    group_key
        Column name in `obs` metadata to group by
    operation
        What operation to perform, one of "mean", "logmean", "median", "std",
        "var", "sem", "mad" and "mad_scaled".
    as_anndata
        Whether to return an AnnData object
    layer
        Which layer to retrieve data from
    store
        Whether to retrieve from/save to cache the result
    progress
        Whether to show a progress bar

    Returns
    -------
    pd.DataFrame
        Result of grouped operation
    """
    keys_tuple = tuple(group_key)
    stored_present = False

    if store:
        if "grouped_ops" not in adata.uns:
            adata.uns["grouped_ops"] = {}

        if keys_tuple not in adata.uns["grouped_ops"]:
            adata.uns["grouped_ops"][keys_tuple] = {}

        if operation in adata.uns["grouped_ops"][keys_tuple]:
            stored_present = True
            res = adata.uns["grouped_ops"][keys_tuple][operation]

    if not stored_present:
        res = grouped_op(
            adata,
            group_key=group_key,
            operation=operation,
            layer=layer,
            progress=progress,
        )

        if store:
            adata.uns["grouped_ops"][keys_tuple][operation] = res

    return grouped_op_to_anndata(res, group_key) if as_anndata else res


def grouped_op_to_anndata(df: pd.DataFrame, group_key: list[str]) -> AnnData:
    """
    Convert a result from a grouped operation into AnnData

    Parameters
    ----------
    df
            Result from grouped operation
    group_key
            Keys used for grouping

    Returns
    -------
    AnnData
            Converted object
    """
    if len(group_key) == 1:
        obs = pd.DataFrame(df.columns, index=df.columns, columns=group_key)
    else:
        obs = pd.DataFrame.from_records(df.columns, columns=group_key)
    obs.index = obs.index.astype(str)
    X = df.T
    X.index = obs.index
    return AnnData(X=X, obs=obs)


def anndata_to_df(adata: AnnData) -> pd.DataFrame:
    """Convert an AnnData object to a pandas DataFrame, keeping .obs"""
    return pd.concat([adata.obs, adata.to_df()], axis=1)


def _doc_params(**kwds):
    """Document common parameters in docstrings."""

    def dec(obj):
        obj.__orig_doc__ = obj.__doc__
        obj.__doc__ = dedent(obj.__doc__).format_map(kwds)
        return obj

    return dec
