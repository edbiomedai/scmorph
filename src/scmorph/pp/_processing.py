from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData
from scanpy._utils import AnyRandom

neighbors = sc.pp.neighbors
neighbors.__doc__ = "| Copied from scanpy [Wolf18]_." + neighbors.__doc__
umap = sc.tl.umap
umap.__doc__ = "| Copied from scanpy [Wolf18]_." + umap.__doc__


def pca(
    adata: AnnData,
    n_comps: Optional[int] = None,
    whiten: bool = False,
    zero_center: Optional[bool] = True,
    svd_solver: str = "arpack",
    random_state: AnyRandom = 0,
    dtype: str = "float32",
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> AnnData:
    adata = adata.copy() if copy else adata
    sc.pp.pca(
        adata,
        n_comps=n_comps,
        zero_center=zero_center,
        svd_solver=svd_solver,
        random_state=random_state,
        dtype=dtype,
        copy=False,
        chunked=chunked,
        chunk_size=chunk_size,
    )

    if whiten:
        adata.obsm["X_pca"] /= np.sqrt(adata.uns["pca"]["variance"])
    return adata if copy else None


pca.__doc__ = (
    "| Copied from scanpy [Wolf18]_ with added whitening.\n" + sc.pp.pca.__doc__
)


def scale(adata: AnnData, chunked: bool = False) -> None:
    """
    Scale data to unit variance per feature while maintaining a low memory footprint

    Parameters
    ----------
    adata :class:`~anndata.AnnData`
            Object as returned by :func:`scmorph.read_cellprofiler`. Represents AnnData object populated with data.

    chunked: bool
            Whether to save memory by processing in chunks. This is slower but less memory intensive.

    Returns
    ----------
    adata : :class:`~anndata.AnnData`
    """

    if not chunked:
        from sklearn.preprocessing import StandardScaler

        scaler = StandardScaler(copy=False)
        if adata.isbacked:
            for i in range(adata.shape[1]):
                scaler.fit_transform(adata[:, i].X)
        else:
            scaler.fit_transform(adata.X)

    else:
        # process one feature at a time
        def scaler(x: np.array) -> None:
            x -= x.mean()
            x /= x.std()

        np.apply_along_axis(scaler, 0, adata.X)


def scale_by_batch(adata: AnnData, batch_key: str, chunked: bool = False) -> None:
    """
    Scale data to unit variance per batch

    Parameters
    ----------
    adata :class:`~anndata.AnnData`
            Object as returned by :func:`scmorph.read_cellprofiler`. Represents AnnData object populated with data.

    batch_key : str
            Name of the column in the AnnData object that contains the batch information.

    chunked : bool
            Whether to save memory by processing in chunks. This is slower but less memory intensive.

    Returns
    -------
    adata : :class:`~anndata.AnnData`

    Note
    -------
    Initial idea taken from https://github.com/scverse/scanpy/issues/2142#issuecomment-1041591406
    """

    for _, idx in adata.obs.groupby(batch_key).indices.items():
        scale(adata[idx, :], chunked=chunked)


def drop_na(
    adata: AnnData,
    feature_threshold: float = 0.9,
    cell_threshold: float = 0,
    inplace: bool = True,
) -> AnnData:
    """
    Drop features with many NAs, then drop cells with any NAs (or infinite values)

    Parameters
    ----------
    adata :class:`~anndata.AnnData`
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
            Rows correspond to cells and columns to genes.

    feature_threshold : float
        Features whose fraction of cells with NA is higher than this will be discarded.

    cell_threshold : float
        Cells whose fraction of features with NA is higher than this will be discarded.

    inplace : bool
        Whether to drop the features and/or cells inplace.
    Returns
    -------
    adata : :class:`~anndata.AnnData`
    """

    isna = np.bitwise_or(np.isinf(adata.X), np.isnan(adata.X))

    if isna.sum() > 0:
        # filter columns where most entries are NaN
        col_mask = isna.sum(axis=0) <= adata.shape[0] * feature_threshold

        # filter cells with any NAs
        row_mask = isna[:, col_mask].sum(axis=1) <= adata.shape[1] * cell_threshold

        if inplace:
            adata._inplace_subset_obs(row_mask)
            adata._inplace_subset_var(col_mask)
        else:
            return adata[row_mask, col_mask].copy()
    else:
        if not inplace:
            return adata
