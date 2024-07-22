from typing import Tuple, Union

import numpy as np
import pandas as pd
from formulaic import Formula
from numpy.lib.stride_tricks import sliding_window_view


def running_median(x: np.ndarray, window: int = 3) -> np.ndarray:
    # follows R's runmed function with endrule="constant"
    if window % 2 == 0:
        raise ValueError("Window must be odd")
    if window == 1:
        return x

    # Pad array to ensure len(m) = len(x)
    k2 = window // 2
    xp = np.pad(x, (k2, k2), mode="edge")
    w = sliding_window_view(xp, window)
    m = np.median(w, axis=1)

    # Emulate R's endrule="constant"
    # I.e. fill edges with the last correct medians
    m[-k2:] = m[-k2 - 1]
    m[:k2] = m[k2]
    return m


# 1. Compute quantiles
def quantile(
    x: np.ndarray, q: np.ndarray, axis: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    # method = "inverted_cdf" to be consistent with qsmooth's implementation
    # this could be adapted by defaulting back to "linear" which
    # may give better interpolation
    quant = np.quantile(x, q, axis=axis, method="inverted_cdf")
    return q, quant


# 2. Quantile regression
def quantile_regression(
    Z: Union[pd.DataFrame, np.ndarray], q: np.ndarray, axis: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    if not isinstance(Z, pd.DataFrame):
        Z = pd.DataFrame(Z, columns=["Z"])

    if not isinstance(q, np.ndarray):
        q = np.array(q)
    # Assumes Z is categorical
    m = np.array(Formula("0 + C(Z)").get_model_matrix(Z))
    mt = m.T
    qbetas = (np.linalg.inv(mt @ m).T @ mt @ q.T).T
    qhat = qbetas @ mt
    qbar = np.mean(q.T, axis=axis)
    return qhat, qbar


# 3. Compute residual components
def residuals(
    f: np.ndarray, fhat: np.ndarray, fbar: float, axis: int = 0
) -> Tuple[float, float, float]:
    SST = np.sum(np.power(f.T - fbar, 2), axis=axis)
    SSB = np.sum(np.power(fhat.T - fbar, 2), axis=axis)
    SSE = SST - SSB
    return SST, SSB, SSE


# 4. Compute weights
def weights(res: Tuple[float, float, float], window: float = 0.05) -> np.ndarray:
    SST, SSB, _ = res

    # Compute rough weights
    roughWeights = 1 - SSB / SST
    roughWeights = np.where(roughWeights < 1e-6, 1, roughWeights)

    # Compute smooth weights
    k = np.floor(window * len(roughWeights))  # type: ignore
    k = k if k % 2 == 1 else k + 1
    k = int(k)
    smoothWeights = running_median(roughWeights, window=k)

    return smoothWeights


# 5. Correct
def apply_correction(fbar: np.ndarray, fhat: np.ndarray, w: float) -> np.ndarray:
    return (np.array(w * fbar + (1 - w)) * fhat.T).T


# 6. Do all
# This implementation differs from qsmooth in that we
# use computed quantiles instead sorting the matrix
# and using that as quantiles. This is because
# we will want to branch out into the multiple-instance case
# where we might have an unequal number of instances per sample
# and therefore sorting and recomputing quantiles in the original
# way would not be possible.
# One can still get the original qsmooth behavior by using
# q = np.linspace(0, 1, X.shape[0])
def quantile_norm(
    X: np.ndarray,
    Z: Union[pd.DataFrame, np.ndarray],
    q: np.ndarray,
    window: float = 0.05,
    axis: int = 0,
) -> np.ndarray:
    quantiles = quantile(X, q, axis=axis)[1]
    qhat, qbar = quantile_regression(Z, quantiles, axis=axis)
    res = residuals(quantiles, qhat, qbar, axis=axis)
    w = weights(res, window=window)
    return apply_correction(qbar, qhat, w)
