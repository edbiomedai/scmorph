import numpy as np
import pandas as pd
from anndata import AnnData
from scanpy.pp import subsample

from .correlation import corr


def corr_features(adata: AnnData, method: str = "pearson", M: int = 5) -> AnnData:
    """
    correlate features and save in `.varm` slot

    Parameters
    ----------
    adata
            The (annotated) data matrix of shape `n_obs` × `n_vars`.
            Rows correspond to cells and columns to genes.

    method
            One of "pearson", "spearman" and "chatterjee" ([:cite:p:`LinHan2021`]_)

    M
            Number of right nearest neighbors to use for Chatterjee correlation.

    Returns
    -------
    adata
        Feature correlations saved in `.varm` slot
    """
    adata.varm[method] = corr(adata.X, method=method, M=M)
    return adata.varm[method]


def _corr_wide_to_long(adata: AnnData, method: str) -> pd.DataFrame:
    x = adata.varm[method].copy()
    n_target = int(len(x) ** 2 / 2)
    xdf = pd.DataFrame(np.tril(x, -1), index=adata.var.index, columns=adata.var.index)
    return xdf.stack().iloc[:n_target]


def _corr_filter(adata: AnnData, method: str, cor_cutoff: float) -> list[str]:
    """Filter pairwise feature correlations, discard features with highest correlation to other features"""
    pairwise_corr = _corr_wide_to_long(adata, method)
    candidate_pairs = pairwise_corr.loc[pairwise_corr.abs() > cor_cutoff].index.to_frame()
    candidate_pairs.columns = ["feature_1", "feature_2"]
    candidate_singlets = np.unique(candidate_pairs.to_numpy().flatten())

    tot_corr = np.abs(adata[:, candidate_singlets].varm[method]).sum(axis=0)
    feat_corrs = pd.Series(tot_corr, index=adata.var.index)

    candidate_pairs["corr_1"] = feat_corrs.loc[candidate_pairs["feature_1"]].values
    candidate_pairs["corr_2"] = feat_corrs.loc[candidate_pairs["feature_2"]].values

    def pick_higher_cor(row: pd.Series) -> str:
        if row["corr_1"] < row["corr_2"]:
            return row["feature_1"]
        else:
            return row["feature_2"]

    exclude = list(set(candidate_pairs.apply(pick_higher_cor, axis=1)))
    return exclude


def select_features(
    adata: AnnData,
    method: str = "pearson",
    cor_cutoff: float = 0.9,
    fraction: float | None = None,
    n_obs: int | None = None,
    copy: bool = False,
) -> AnnData | None:
    """
    Feature selection

    Select features by feature correlations. Allows measuring correlation
    on a subset of cells to speed up computations. See ``fraction`` and ``n_obs``
    for details.

    Parameters
    ----------
    adata
            The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
            Rows correspond to cells and columns to genes.

    method
            Which correlation coefficient to use for filtering.
            One of "pearson", "spearman" and "chatterjee" ([:cite:p:`LinHan2021`]_)

    cor_cutoff
            Cutoff beyond which features with a correlation coefficient
            higher than it are removed. Must be between 0 and 1.

    fraction: float
            Subsample to this ``fraction`` of the number of observations.

    n_obs
            Subsample to this number of observations.

    copy
            Whether to return a copy or modify ``adata`` inplace
            (i.e. operate inplace)

    Returns
    -------
    adata
        Feature correlations saved in ``.varm`` slot
        and feature selection saved in ``.var`` slot.

    """
    # sampling
    if fraction or n_obs:
        adata_ss = subsample(adata, fraction=fraction, n_obs=n_obs, copy=True)
    else:
        adata_ss = adata

    # variance filter
    pass_var = np.empty(len(adata.var), dtype=bool)

    for i, feat in enumerate(adata_ss.var_names):
        pass_var[i] = False if np.var(adata_ss[:, feat].X) < 1e-5 else True

    adata.var["qc_pass_var"] = pass_var
    adata_ss.var["qc_pass_var"] = pass_var

    # correlation filter
    corr_features(adata_ss, method=method, M=5)
    adata.varm[method] = adata_ss.varm[method]
    exclude = _corr_filter(adata_ss, method=method, cor_cutoff=cor_cutoff)
    keep = np.invert(adata.var.index.isin(exclude))

    if not copy:
        adata._inplace_subset_var(keep)
        return None
    else:
        return adata[:, keep].copy()
