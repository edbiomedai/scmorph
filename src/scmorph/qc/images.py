import numpy as np
import pandas as pd
from anndata import AnnData

import scmorph as sm


def _embed(adata):
    sm.pp.scale(adata)
    sm.pp.pca(adata)


def _get_pca(adata):
    if "X_pca" not in adata.obsm_keys():
        _embed(adata)
    return adata.obsm["X_pca"] * adata.uns["pca"]["variance_ratio"]


def kNN_dists(adata: AnnData, pcs: int = 3, neighbors: int = 10):
    """
    Compute maximum kNN distance (i.e. radius of smallest enclosing circle of kNNs)

    Parameters
    ----------
    adata
        image-level data
    pcs
        Number of PCs to use
    neighbors
        Number of image neigbors in PC

    Returns
    -------
    For each image, how far is the k-th nearest neighbor away in PC space (measured as Euclidean distance)
    """
    from sklearn.neighbors import NearestNeighbors

    a = _get_pca(adata)[:, :pcs]
    nbrs = NearestNeighbors(n_neighbors=neighbors + 1).fit(a)
    d, _ = nbrs.kneighbors(a)
    dss = d[:, 1:]
    dist = dss[:, neighbors - 1]
    return dist / np.sqrt(pcs)


def unsupervised_imageQC(qcadata: AnnData, pcs: int = 3, neighbors: int = 10):
    """
    Compute maximum kNN distance (i.e. radius of smallest enclosing circle of kNNs).
    This function will perform center-scaling and PCA transform, before computing distances.
    It also saves the image-level PCA in obsm["X_pca"].

    Parameters
    ----------
    qcadata
        image-level data
    pcs
        Number of PCs
    neighbors
        Number of image neigbors in PC

    Returns
    -------
    Image-level data with added ImageQCDistance in `obs` and `X_pca` in `obsm`. Does not operate on X in-place.
    """
    qcadatac = qcadata.copy()
    dists = kNN_dists(qcadatac, pcs, neighbors)
    qcadata.obs["ImageQCDistance"] = dists
    qcadata.obsm["X_pca"] = qcadatac.obsm["X_pca"]
    qcadata.uns["pca"] = qcadatac.uns["pca"]
    qcadata.varm["PCs"] = qcadatac.varm["PCs"]
    return qcadata


def _filter_adata_by_qc(scadata: AnnData, qcadata: AnnData):
    return (
        scadata[scadata.obs["PassQC"] == "True"].copy(),
        qcadata[qcadata.obs["PassQC"] == "True"].copy(),
    )


def qc_images_by_dissimilarity(
    adata: AnnData,
    qcadata: AnnData,
    filter: bool = True,
    threshold: float = 0.05,
    **kwargs,
) -> tuple[AnnData, AnnData]:
    """
    Perform QC of datasets using unsupervised, kNN-based distance filtering

    Parameters
    ----------
    adata
        Single-cell data
    qcadata
        Image-level data
    filter
        Whether to return filtered or unfiltered (i.e. only annotated) adatas
    threshold
        Threshold for removal

    Returns
    -------
    Tuple of single-cell and image-level adatas
    """
    merge_cols = adata.obs.columns.isin(qcadata.obs.columns)
    merge_cols = adata.obs.columns[merge_cols]
    merge_cols = merge_cols.drop(["PassQC", "ImageQCDistance"], errors="ignore")
    merge_cols = merge_cols.to_list()
    assert len(merge_cols) > 0, "No columns to merge on"
    qcadata = unsupervised_imageQC(qcadata, **kwargs)
    qcadata.obs["PassQC"] = qcadata.obs["ImageQCDistance"] < threshold
    qcadata.obs["PassQC"] = qcadata.obs["PassQC"].astype(str).astype("category")
    new = pd.merge(adata.obs, qcadata.obs, how="left", on=merge_cols)
    new.index = adata.obs.index
    adata.obs = new
    if not filter:
        return adata, qcadata
    return _filter_adata_by_qc(adata, qcadata)


def count_cells_per_group(
    adata: AnnData,
    group_keys: list[str],
    inplace: bool = True,
):
    """
    Count number of cells per group.

    Parameters
    ----------
    adata
        Annotated data matrix.
    group_keys
        List of keys to group by, typically plate, well and site or image ID.
    inplace
        Whether to modify the AnnData object in place.

    Returns
    -------
    Annotated data matrix with an additional column "cells_per_group" in `adata.obs` containing the number of cells group.
    """
    if not inplace:
        adata = adata.copy()
    if not isinstance(group_keys, list):
        group_keys = [group_keys]
    obs = adata.obs
    new_obs = (
        obs.groupby(group_keys, observed=True)
        .apply(len, include_groups=False)
        .rename("cells_per_group")
        .to_frame()
        .reset_index()
        .merge(obs, how="right")
        .set_index(obs.index)
    )
    adata.obs = new_obs
    return adata
