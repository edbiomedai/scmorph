from typing import Optional, Tuple

import numpy as np
from anndata import AnnData
from scanpy.pp import subsample

from ._correlation import compute_corr


def corr_features(adata: AnnData, method: str = "pearson", M: int = 5) -> AnnData:
    """
    correlate features and save in `.varm` slot

    Parameters
    ----------
    adata : AnnData
            The (annotated) data matrix of shape `n_obs` × `n_vars`.
            Rows correspond to cells and columns to genes.

    method : str, optional
            One of "pearson", "spearman" and "chatterjee", by default "pearson"

    M : int, optional
            Number of right nearest neighbors to use for chatterjee correlation.

    Returns
    -------
    AnnData
        AnnData object with feature correlations saved in `.varm` slot
    """
    adata.varm[method] = compute_corr(adata.X, method=method, M=M)
    return adata.varm[method]


def select_features(
    adata: AnnData,
    cor_cutoff: Tuple[float, float, float] = (0.8, 0.8, 0.8),
    fraction: Optional[float] = None,
    n_obs: Optional[int] = None,
) -> AnnData:
    """
    Select features by feature correlations. Allows measuring correlation
    on a subset of cells to speed up computations. See `fraction` and `n_obs`
    for details.

    Parameters
    ----------
    adata : AnnData
            The (annotated) data matrix of shape `n_obs` × `n_vars`.
            Rows correspond to cells and columns to genes.

    cor_cutoff : tuple, optional
            Cutoff beyond which features with a correlation coefficient
            higher than it are removed. Expects three values for Pearson,
            Spearman and Chatterjee correlation coefficient cutoffs, respectively.

    fraction: float, optional
            Subsample to this `fraction` of the number of observations.

    n_obs : int, optional
            Subsample to this number of observations.

    Returns
    -------
    AnnData
        AnnData object with feature correlations saved in `.varm` slot
        and feature selection saved in `.var` slot.

    """

    # sampling
    adata_ss = subsample(adata, fraction=fraction, n_obs=n_obs, copy=True)

    # variance filter
    pass_var = np.empty(len(adata.var), dtype=bool)

    for i, feat in enumerate(adata_ss.var_names):
        pass_var[i] = False if np.var(adata_ss[:, feat].X) < 1e-5 else True
    adata.var["qc_pass_var"] = pass_var
    adata[:, adata.var["qc_pass_var"]]

    # correlation filter
    for m in ["pearson", "spearman", "chatterjee"]:
        compute_corr(adata_ss.X, method=m, M=5)
        adata.varm[m] = adata_ss.varm[m]

    # TODO implement logic to subset features based on correlation
    # see: https://github.com/cytomining/pycytominer/blob/02d05223effad9c12239ba331a8ce709153b18a4/pycytominer/operations/correlation_threshold.py#L80
