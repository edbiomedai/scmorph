import re
from functools import partial
from inspect import signature
from typing import Any, Callable, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from scmorph.logging import get_logger


def _infer_names(type: str, options: Iterable[str]) -> Sequence[str]:
    logger = get_logger()

    if type == "batch":
        reg = re.compile("batch|plate")
    elif type in ["well", "group"]:
        reg = re.compile("well$")
    elif type == "treatment":
        reg = re.compile("treatment", re.IGNORECASE)
    elif type == "site":
        reg = re.compile("site$")
    else:
        raise ValueError("type must be one of 'batch', 'well', 'treatment', 'site'")
    res = [x for x in options if reg.search(x)]
    if len(res) > 1:
        logger.warning(
            f"Found multiple {type} columns, are these duplicates?\n"
            + "Will just use the first one, if that is not desired please specify "
            + f"the correct column name using the {type}_key argument.\n"
            + f"Problematic columns were: {', '.join(res)}"
        )
        res = [res[0]]
    return res


def grouped_obs_fun(
    adata: AnnData,
    group_key: Union[str, List[str]],
    fun: Callable[..., Any],
    layer: Optional[str] = None,
) -> pd.DataFrame:
    """
    Slightly adapted from https://github.com/scverse/scanpy/issues/181#issuecomment-534867254
    All copyright lies with Isaac Virshup.
    """

    def getX(adata: AnnData, layer: Union[None, str]) -> np.array:
        return adata.X if layer is None else adata.layers[layer]

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names,
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx], layer)
        out[group] = np.ravel(fun(X))
    return out


def grouped_op(
    adata: AnnData,
    group_key: Union[str, List[str]],
    operation: str,
    layer: Optional[str] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    if operation == "mean":
        fun = partial(np.mean, axis=0, dtype=np.float64, **kwargs)
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

        f1 = partial(mad, axis=0, **kwargs)
        f2 = partial(np.median, axis=0, **kwargs)

        def fun(x: np.array) -> np.array:  # type: ignore
            return f2(x) / f1(x) * 1.4826  # Chung 2008

    else:
        raise ValueError(
            "Operation must be one of 'mean', 'median', 'std', 'var', 'sem'"
        )
    return grouped_obs_fun(adata, group_key, fun=fun, layer=layer)


def group_obs_fun_inplace(
    adata: AnnData,
    group_key: Union[str, List[str]],
    fun: Callable[..., Any],
) -> AnnData:
    """
    Alter adata.X inplace by performing fun in each group

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix object

    group_key : Union[str, List[str]]
        obs keys to group by

    fun : Callable
        Function that takes array and returns array of equal size.
        The function may either only take an array, or the array and the group key.
        In the latter case, the group key must be the second argument!

    Returns
    -------
    AnnData
        The annotated data matrix object after the operation
    """
    grouped = adata.obs.groupby(group_key)

    takes_group = len(signature(fun).parameters) > 1

    for group, idx in grouped.indices.items():
        X = adata[idx].X
        adata[idx].X = fun(X, group) if takes_group else fun(X)

    return adata


def get_group_keys(
    adata: AnnData,
    treatment_key: Union[str, List[str]],
    group_key: Optional[Union[str, List[str]]],
) -> Tuple[List[str], List[str]]:
    # convert to list
    if not isinstance(group_key, list):
        group_col = [group_key]

    # inferring treatment names if necessary
    if isinstance(treatment_key, str):
        if treatment_key == "infer":
            treatment_col = _infer_names("treatment", adata.obs.columns)
        else:
            treatment_col = [treatment_key]

    if isinstance(group_key, str):
        if group_key == "infer":
            group_col = _infer_names("group", adata.obs.columns)  # type: ignore
        else:
            group_col = [group_key]
    # end inference

    group_keys = [*treatment_col, *group_col]
    group_keys = [x for x in group_keys if x]  # remove None's
    return group_keys, treatment_col  # type: ignore


def get_grouped_op(
    adata: AnnData,
    group_key: List[str],
    operation: str,
    as_anndata: bool = False,
    layer: Optional[str] = None,
    store: bool = True,
) -> pd.DataFrame:
    """
    Retrieve from cache or compute a grouped operation

    Parameters
    ----------
    adata : AnnData
        AnnData object
    group_key : List[str]
        Grouping keys
    operation : str
        Operation, e.g. mean or median
    as_anndata : bool,
        Whether to return an AnnData object, by default False
    layer : Optional[str], optional
        Which layer to retrieve data from, by default None
    store : bool, optional
        Whether to retrieve from/save to cache the result, by default True

    Returns
    -------
    pd.DataFrame
        Result of grouped operation
    """

    keys_tuple = tuple(group_key)
    stored_present = False

    if store:
        if "grouped_ops" not in adata.uns:
            adata.uns["grouped_ops"] = dict()

        if keys_tuple not in adata.uns["grouped_ops"]:
            adata.uns["grouped_ops"][keys_tuple] = dict()

        if operation in adata.uns["grouped_ops"][keys_tuple]:
            stored_present = True
            res = adata.uns["grouped_ops"][keys_tuple][operation]

    if not stored_present:
        res = grouped_op(adata, group_key, operation, layer)

        if store:
            adata.uns["grouped_ops"][keys_tuple][operation] = res

    if as_anndata:
        return grouped_op_to_anndata(res, group_key)

    return res


def grouped_op_to_anndata(df: pd.DataFrame, group_key: List[str]) -> AnnData:
    """
    Convert a result from a grouped operation into AnnData

    Parameters
    ----------
    df : pd.DataFrame
            Result from grouped operation
    group_key : List[str]
            Keys used for grouping

    Returns
    -------
    AnnData
            Converted object
    """
    obs = pd.DataFrame.from_records(df.columns, columns=group_key)
    obs.index = obs.index.astype(str)
    X = df.T
    X.index = obs.index
    return AnnData(X=X, obs=obs)


def anndata_to_df(adata: AnnData) -> pd.DataFrame:
    """
    Convert an AnnData object to a pandas DataFrame, keeping .obs
    """
    return pd.concat([adata.obs, adata.to_df()], axis=1)
