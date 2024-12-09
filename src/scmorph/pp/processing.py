from typing import Any

import numpy as np
import scanpy as sc
from anndata import AnnData


def neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    n_pcs: int | None = None,
    use_rep: str = "X_pca",
    copy: bool = False,
    **kwargs: Any,
) -> AnnData | None:
    """
    Compute a neighborhood graph of observations using the PCA representation.

    This function wraps the higher-level function :func:`~scanpy.pp.neighbors` (:cite:p:`Wolf18`.)

    Parameters
    ----------
    adata
        AnnData object

    n_neighbors
        The size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation. Larger values result in more
        global views of the manifold, while smaller values result in more local
        data being preserved. In general values should be in the range 2 to 100.
        If `knn` is `True`, number of nearest neighbors to be searched. If `knn`
        is `False`, a Gaussian kernel width is set to the distance of the
        `n_neighbors` neighbor.

    n_pcs
        Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is not None`.

    use_rep
        Use the indicated representation. 'X' or any key for .obsm is valid.
        If None, the representation is chosen automatically:
        For .n_vars < N_PCS, .X is used, otherwise `X_pca` is used.
        If `X_pca` is not present, it`s computed with default parameters or n_pcs if present.

    copy
        Return a copy instead of writing to adata.

    kwargs
        Additional arguments passed to :func:`~scanpy.pp.neighbors`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obsp['distances' | key_added+'_distances']` : :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Distance matrix of the nearest neighbors search. Each row (cell) has `n_neighbors`-1 non-zero entries. These are the distances to their `n_neighbors`-1 nearest neighbors (excluding the cell itself).
    `adata.obsp['connectivities' | key_added+'_connectivities']` : :class:`scipy.sparse._csr.csr_matrix` (dtype `float`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    `adata.uns['neighbors' | key_added]` : :class:`dict`
        neighbors parameters.
    """
    adata = adata.copy() if copy else adata
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep, copy=False, **kwargs)
    return adata if copy else None


def umap(adata: AnnData, **kwargs) -> AnnData | None:
    """
    Embed the neighborhood graph using UMAP :cite:p:`McInnes2018`.

    This function wraps the higher-level function :func:`~scanpy.tl.umap` (:cite:p:`Wolf18`.)

    Parameters
    ----------
    adata
        AnnData object

    kwargs
        Additional arguments passed to :func:`~scanpy.tl.umap`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obsm['X_umap' | key_added]` : :class:`numpy.ndarray` (dtype `float`)
        UMAP coordinates of data.
    `adata.uns['umap' | key_added]` : :class:`dict`
        UMAP parameters.
    """
    return sc.tl.umap(adata, **kwargs)


def pca(
    adata: AnnData,
    n_comps: int | None = None,
    whiten: bool = False,
    *,
    copy: bool = False,
    zero_center: bool | None = True,
    random_state: int | None = 0,
    **kwargs: Any,
) -> AnnData | None:
    """
    Principal component analysis :cite:p:`Pedregosa2011`.

    Computes PCA coordinates, loadings and variance decomposition.
    Uses the implementation of *scikit-learn* :cite:p:`Pedregosa2011`.

    This function wraps the higher-level function :func:`~scanpy.pp.pca` (:cite:p:`Wolf18`).

    It provides addional functionality to whiten the resulting PCA coordinates, which
    may help de-correlate them.

    Parameters
    ----------
    adata
        AnnData object

    n_comps
        Number of principal components to compute. Defaults to 50, or 1 - minimum
        dimension size of selected representation.

    whiten
        Whether to whiten the data. This can help remove correlation in the resulting PC axes.

    copy
        Return a copy instead of writing to adata.

    zero_center
        If `True`, compute standard PCA from covariance matrix (strongly recommended).
        If `False`, omit zero-centering variable.

    random_state
        Change to use different initial states for the optimization.

    kwargs
        Additional arguments passed to :func:`~scanpy.pp.pca`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `.obsm['X_pca' | key_added]` : :class:`~scipy.sparse.spmatrix` | :class:`~numpy.ndarray` (shape `(adata.n_obs, n_comps)`)
        PCA representation of data.
    `.varm['PCs' | key_added]` : :class:`~numpy.ndarray` (shape `(adata.n_vars, n_comps)`)
        The principal components containing the loadings.
    `.uns['pca' | key_added]['variance_ratio']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Ratio of explained variance.
    `.uns['pca' | key_added]['variance']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Explained variance, equivalent to the eigenvalues of the
        covariance matrix.
    """
    adata = adata.copy() if copy else adata
    sc.pp.pca(
        adata,
        n_comps=n_comps,
        zero_center=zero_center,
        random_state=random_state,
        copy=False,
    )

    if whiten:
        adata.obsm["X_pca"] /= np.sqrt(adata.uns["pca"]["variance"])
    return adata if copy else None


def scale(adata: AnnData, treatment_key: str | None = None, control: str | None = None, chunked: bool = False) -> None:
    """
    Scale data to unit variance per feature while maintaining a low memory footprint (operates in-place).

    Parameters
    ----------
    adata
            Annotated data matrix.

    batch_key
            Name of the column in the AnnData object that contains the batch information.

    treatment_key
            Name of column used to delinate treatments. This is used when computing batch effects across drug-treated plates.
            In that case, we compute batch effects only on untreated cells and then apply the correction factors to all cells.
            If using, please also see `control`.

    chunked
            Whether to save memory by processing in chunks. This is slower but less memory intensive.

    """
    if not chunked:
        from sklearn.preprocessing import StandardScaler

        scaler = StandardScaler()
        if adata.isbacked:
            for i in range(adata.shape[1]):
                if treatment_key and control:
                    X_ctrl = adata[adata.obs[treatment_key] == control, i].X
                    scaler.fit(np.array(X_ctrl))
                else:
                    scaler.fit(np.array(adata[:, i].X))
                adata[:, i] = scaler.transform(np.array(adata[:, i].X))
        else:
            if treatment_key and control:
                X_ctrl = adata[adata.obs[treatment_key] == control, :].X
                scaler.fit(np.array(X_ctrl))
                adata.X = scaler.transform(np.array(adata.X))
            else:
                adata.X = scaler.fit_transform(np.array(adata.X))

    else:
        # process one feature at a time
        def scale_func(x: np.ndarray) -> None:
            if treatment_key and control:
                X_ctrl = x[adata.obs[treatment_key] == control, :]
                x_mean = X_ctrl.mean()
                x_std = X_ctrl.std()
            else:
                x_mean = x.mean()
                x_std = x.std()

            x -= x_mean
            x /= x_std

        np.apply_along_axis(scale_func, 0, adata.X)  # type: ignore


def scale_by_batch(
    adata: AnnData, batch_key: str, treatment_key: str | None = None, control: str | None = None, chunked: bool = False
) -> None:
    """
    Scale data to zero-center and unit variance per batch in-place.

    Parameters
    ----------
    adata
            Annotated data matrix.

    batch_key
            Name of the column in the AnnData object that contains the batch information.

    treatment_key
            Name of column used to delinate treatments. This is used when computing batch effects across drug-treated plates.
            In that case, we compute batch effects only on untreated cells and then apply the correction factors to all cells.
            If using, please also see `control`.

    control
            Name of control treatment. Must be valid value in `treatment_key`.

    chunked
            Whether to save memory by processing in chunks. This is slower but less memory intensive.
    """
    for _, idx in adata.obs.groupby(batch_key, observed=True).indices.items():
        scale(adata[idx, :], treatment_key=treatment_key, control=control, chunked=chunked)


def drop_na(
    adata: AnnData,
    feature_threshold: float = 0.9,
    cell_threshold: float = 0,
    inplace: bool = True,
) -> None | AnnData:
    """
    Drop features with many NAs, then drop cells with any NAs (or infinite values)

    Parameters
    ----------
    adata
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.

    feature_threshold
        Features whose fraction of cells with NA is higher than this will be discarded.

    cell_threshold
        Cells whose fraction of features with NA is higher than this will be discarded.

    inplace
        Whether to drop the features and/or cells inplace.

    Returns
    -------
    Depending on `inplace`, returns or updates `adata` with the filtered data.
    """
    X_dense = adata.X if isinstance(adata.X, np.ndarray) else np.array(adata.X)

    isna = np.bitwise_or(np.isinf(X_dense), np.isnan(X_dense))

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


def leiden(
    adata: AnnData,
    resolution: float = 1.0,
    copy: bool = False,
    **kwargs: Any,
) -> AnnData | None:
    """
    Cluster cells into subgroups :cite:p:`Traag2019`.

    Cluster cells using the Leiden algorithm :cite:p:`Traag2019`,
    an improved version of the Louvain algorithm :cite:p:`Blondel2008`.
    It has been proposed for single-cell analysis by :cite:t:`Levine2015`.

    This requires having ran :func:`~.pp.neighbors` or
    :func:`~scanpy.external.pp.bbknn` first.

    This function wraps the higher-level function :func:`~scanpy.tl.leiden` (:cite:p:`Wolf18`).

    Parameters
    ----------
    adata
        AnnData object

    resolution
        A parameter value controlling the coarseness of the clustering.
        Higher values lead to more clusters.
        Set to `None` if overriding `partition_type`
        to one that doesn't accept a `resolution_parameter`.

    copy
        Return a copy instead of writing to adata.

    kwargs
        Additional arguments passed to :func:`~scanpy.tl.leiden`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obs['leiden' | key_added]` : :class:`pandas.Series` (dtype ``category``)
        Array of dim (number of samples) that stores the subgroup id
        (``'0'``, ``'1'``, ...) for each cell.

    `adata.uns['leiden' | key_added]['params']` : :class:`dict`
        A dict with the values for the parameters `resolution`, `random_state`,
        and `n_iterations`.
    """
    return sc.tl.leiden(adata, resolution=resolution, copy=copy, **kwargs)
