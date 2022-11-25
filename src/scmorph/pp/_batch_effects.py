from typing import Any, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

from scmorph._utils import _infer_names, get_grouped_op, group_obs_fun_inplace

"""Functions to remove batch effects from morphological datasets."""


def compute_batch_effects(
    adata: AnnData,
    bio_key: Optional[str] = None,
    batch_key: str = "infer",
    treatment_key: Optional[str] = None,
    control: str = "DMSO",
    progress: bool = True,
) -> Tuple[np.array, np.array]:
    """
    Compute batch effects using scone's method of deconvoluting technical from biological effects
    through linear modeling. If no biological differences are expected (e.g. because data is all from
    the same cell line), set `bio_key` to None.

    For details on scone, see [Cole19]_. For details on the results particularly view Eq. 3.

    Parameters
    ----------
    adata :class:`~anndata.AnnData`
            Annotated data matrix

    bio_key : str
            Name of column used to delineate biological entities, e.g. cell lines. Default: None

    batch_key : str
            Name of column used to delineate batch effects, e.g. plates. Will try to guess if
            no argument is given. Default: "infer"

    treatment_key: str
            Name of column used to delinate treatments. This is used when computing batch effects across drug-treated plates.
            In that case, we compute batch effects only on untreated cells and then apply the correction factors to all cells.
            If using, please also see `control`.

    control: str
            Name of control treatment. Must be valid value in `treatment_key`.

    progress: bool
            Whether to show a progress bar, by default True

    Returns
    -------
    betas : :class:`~numpy.array`
        Biological effects, i.e. how much each feature varied because of biological differences
    gammas: :class:`~numpy.array`
        Technical effects, i.e. batch effects.
    """
    from formulaic import Formula
    from patsy.contrasts import Treatment

    withbio = bio_key is not None

    if batch_key == "infer":
        batch_key = _infer_names("batch", adata.obs.columns)[0]

    # combine keys for batch and bio
    joint_keys = [i for i in [bio_key, batch_key] if i is not None]

    if treatment_key is not None:
        adata = adata[adata.obs[treatment_key] == control, :]

    data = get_grouped_op(
        adata, group_key=joint_keys, operation="mean", progress=progress
    )  # compute average feature per batch/bio group

    data = data.loc[adata.var.index]  # ensure correct feature order

    groups = pd.DataFrame(data.columns.to_list(), columns=joint_keys)

    # create design matrix
    if withbio:
        dmatrix = Formula(f"C({bio_key}) +C({batch_key})").get_model_matrix(groups)
    else:
        dmatrix = Formula(f"C({batch_key})").get_model_matrix(groups)

    # identify model parameters using linear modeling
    # rcond=None to avoid FutureWarning
    params = np.linalg.lstsq(dmatrix, data.T, rcond=None)[0]

    # extract parameters into beta and gammas
    gamma_ind = np.where(dmatrix.columns.str.contains(batch_key))[0]

    batch_groups = groups[batch_key].drop_duplicates().to_list()
    contrast_gamma = Treatment(reference=0).code_without_intercept(batch_groups)
    gammas = (contrast_gamma.matrix @ params[gamma_ind, :]).T

    betas = []
    if withbio:
        beta_ind = np.where(dmatrix.columns.str.contains(bio_key))[0]
        bio_groups = groups[bio_key].drop_duplicates().to_list()
        contrast_beta = Treatment(reference=0).code_without_intercept(bio_groups)
        betas = (contrast_beta.matrix @ params[beta_ind, :]).T

    # convert to DataFrames to store metadata
    gammas_df = pd.DataFrame(gammas, columns=sorted(batch_groups), index=data.index)
    betas_df = (
        pd.DataFrame(betas, columns=sorted(bio_groups), index=data.index)
        if len(betas) != 0
        else []
    )
    return betas_df, gammas_df


def remove_batch_effects(
    adata: AnnData,
    bio_key: Optional[str] = None,
    batch_key: str = "infer",
    copy: Optional[bool] = False,
    **kwargs: Any,
) -> AnnData:
    """
    Remove batch effects using scone's method of deconvoluting technical from biological effects
    through linear modeling. Note that this preserves biological differences and only
    removes technical effects. If no biological differences are expected (e.g. because data
    is all from the same cell line), set `bio_key` to None.

    For details on scone, see [Cole19]_.

    Parameters
    ----------
    adata :class:`~anndata.AnnData`
            Annotated data matrix

    bio_key : _type_
            Name of column used to delineate biological entities, e.g. cell lines. Default: None

    batch_key : str
            Name of column used to delineate batch effects, e.g. plates. Will try to guess if
            no argument is given. Default: "infer"

    copy : bool
            If False, will perform operation in-place, else return a modified copy of the data.

    kwargs: dict
            Other parameter passed to :func:`compute_batch_effects`

    Returns
    -------
    adata : :class:`~anndata.AnnData`
            Annotated data matrix with batch effects removed. If `copy` is False, will modify in-place and not return anything.
    """
    if copy:
        adata = adata.copy()
    betas, gammas = compute_batch_effects(
        adata, **kwargs, bio_key=bio_key, batch_key=batch_key
    )
    adata.uns["batch_effects"] = (
        pd.concat((betas, gammas), axis=1) if len(betas) > 0 else gammas
    )
    group_obs_fun_inplace(
        adata, batch_key, lambda x, group: x - gammas[group].to_numpy()
    )
    if copy:
        return adata
