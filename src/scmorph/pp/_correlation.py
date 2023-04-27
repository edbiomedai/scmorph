from typing import Iterator, Tuple, Union

import numpy as np
from numba import jit
from scipy.stats import spearmanr


def _iter_cols(
    X: np.array, upper: bool = False, return_arr: bool = False
) -> Union[Iterator[Tuple[int, int]], Iterator[np.array]]:
    """Iterate over columns of matrix, return indices. Can also return view of array"""
    cols = range(X.shape[1])
    for col1 in cols:
        for col2 in cols:
            if upper and col1 >= col2:
                continue
            if not return_arr:
                yield (col1, col2)
            else:
                yield X[:, [col1, col2]]


def xim(X: np.array, Y: Union[None, np.array] = None, M: int = 5) -> np.array:
    """
    Compute the XIM (revised Chatterjee) correlation coefficient

    Parameters
    ----------
    X, Y: 1D or 2D array_like, Y is optional
            One or two 1-D or 2-D arrays containing multiple variables and observations.
            When these are 1-D, each represents a vector of observations of a single variable.
            In the 2-D case, each row is assumed to contain an observation.
            Both arrays need to have the same length.

    M : int
            Number of right nearest neighbors

    Returns
    ----------
    xim : :class:~`float`
            Value of XIM

    Note
    ----------
    This code is an implementation of [Lin21]_.
    Implementation by Jesko Wagner.
    """

    @jit(nopython=True)
    def _comp_coef(xorder: np.array, yrank: np.array, M: int) -> float:
        """helper for XIM correlation"""
        n = yrank.size
        yrank = yrank[xorder]
        coef_sum = 0
        for m in range(1, M + 1):
            coef_sum_temp = np.sum(np.minimum(yrank[: (n - m)], yrank[m:n]) + 1)
            coef_sum_temp = coef_sum_temp + np.sum(yrank[(n - m) : n] + 1)
            coef_sum = coef_sum + coef_sum_temp
        return -2 + 6 * coef_sum / ((n + 1) * (n * M + M * (M + 1) / 4))

    if len(X.shape) > 1:
        if Y is not None:
            raise ValueError("Y must not be provided if X is a 2D array")
    elif Y is None:
        raise ValueError("Y must be provided if X is a 1D array")
    elif Y.shape != X.shape:
        raise ValueError("X and Y must have the same shape")
    else:
        X = np.column_stack((X, Y))

    # pseudorandom tie breaker based on adding random noise
    X = X + np.random.uniform(-0.000001, 0.000001, X.shape)
    orders = X.argsort(axis=0)
    ranks = orders.argsort(axis=0)

    cormat_1d = np.array(
        [
            _comp_coef(orders[:, i], ranks[:, j], M=M)
            for i, j in _iter_cols(X, upper=False)
        ]
    )

    cormat_2d = np.eye(X.shape[1])
    cormat_2d[np.triu_indices(cormat_2d.shape[0], k=1)] = cormat_1d
    cormat_2d[np.tril_indices(cormat_2d.shape[0], k=-1)] = cormat_1d.T

    return cormat_2d


def corr(
    X: np.array, Y: Union[np.array, None] = None, method: str = "pearson", M: int = 5
) -> float:
    """
    Compute pairwise correlations

    Parameters
    ----------
    X, Y : 1D or 2D array_like, Y is optional
            One or two 1-D or 2-D arrays containing multiple variables and observations.
            When these are 1-D, each represents a vector of observations of a single variable.
            In the 2-D case, each row is assumed to contain an observation.
            Both arrays need to have the same length.

    method : str
            One of "pearson", "spearman", or "chatterjee" ([Lin21]_), by default "pearson"

    M : int
            Number of right nearest neighbors to use for Chatterjee correlation.

    Returns
    -------
    corr : :class:~`float`
            Correlation coefficient
    """

    if method == "pearson":
        result = np.corrcoef(X, Y, False)
    elif method == "spearman":
        result = spearmanr(X, Y)[0]
    elif method == "chatterjee":
        result = xim(X, Y, M=M)
    else:
        raise ValueError(
            f"Method must be one of 'pearson', 'spearman', or 'chatterjee'. Received {method}"
        )
    return result
