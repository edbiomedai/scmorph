from collections.abc import Iterator

import numpy as np
from numba import jit
from scipy.stats import spearmanr


def _iter_cols(
    X: np.ndarray, upper: bool = False, return_arr: bool = False
) -> Iterator[tuple[int, int]] | Iterator[np.ndarray]:
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


def xim(X: np.ndarray, Y: None | np.ndarray = None, M: int = 5) -> np.ndarray:
    """
    Compute the XIM (revised Chatterjee) correlation coefficient

    Parameters
    ----------
    X, Y: np.ndarray
            One or two 1-D or 2-D arrays containing multiple variables and observations.
            When these are 1-D, each represents a vector of observations of a single variable.
            In the 2-D case, each row is assumed to contain an observation.
            Both arrays need to have the same length.

    M
            Number of right nearest neighbors

    Returns
    -------
    xim
            Value of XIM

    Note
    ----------
    This code is an implementation of [:cite:p:`LinHan2021`]_.
    Implementation by Jesko Wagner.
    """

    @jit(nopython=True)
    def _comp_coef(xorder: np.ndarray, yrank: np.ndarray, M: int) -> float:
        """Helper for XIM correlation"""
        n = yrank.size
        yrank = yrank[xorder]
        coef_sum = 0
        for m in range(1, M + 1):
            coef_sum_temp = np.sum(np.minimum(yrank[: (n - m)], yrank[m:n]) + 1)
            coef_sum_temp = coef_sum_temp + np.sum(yrank[(n - m) : n] + 1)
            coef_sum = coef_sum + coef_sum_temp
        return float(-2 + 6 * coef_sum / ((n + 1) * (n * M + M * (M + 1) / 4)))

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
    rng = np.random.default_rng(seed=2024)
    X = X + rng.uniform(-0.000001, 0.000001, X.shape)
    orders = X.argsort(axis=0)
    ranks = orders.argsort(axis=0)

    # compute 1D correlation vector of pairwise comparisons
    cormat_1d = np.array([_comp_coef(orders[:, i], ranks[:, j], M=M) for i, j in _iter_cols(X, upper=False)])

    # reshape to 2D correlation matrix
    cormat_2d = cormat_1d.reshape(X.shape[1], X.shape[1])

    # force xim=1 for diagonal elements
    np.fill_diagonal(cormat_2d, 1)

    return cormat_2d


def corr(X: np.ndarray, Y: np.ndarray | None = None, method: str = "pearson", M: int = 5) -> np.ndarray:
    """
    Compute pairwise correlations

    Parameters
    ----------
    X, Y
            One or two 1-D or 2-D arrays containing multiple variables and observations.
            When these are 1-D, each represents a vector of observations of a single variable.
            In the 2-D case, each row is assumed to contain an observation.
            Both arrays need to have the same length.

    method
            One of "pearson", "spearman", or "chatterjee" ([:cite:p:`LinHan2021`]_)

    M
            Number of right nearest neighbors to use for Chatterjee correlation.

    Returns
    -------
    corr
            Correlation coefficient
    """
    if method == "pearson":
        result = np.corrcoef(X, Y, False)
    elif method == "spearman":
        result = spearmanr(X, Y)[0]
    elif method == "chatterjee":
        result = xim(X, Y, M=M)
    else:
        raise ValueError(f"Method must be one of 'pearson', 'spearman', or 'chatterjee'. Received {method}")
    return result
