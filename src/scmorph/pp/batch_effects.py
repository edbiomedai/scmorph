from typing import Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

from scmorph._utils import _infer_names, group_obs_fun_inplace, grouped_op

"""Functions to remove batch effects from morphological datasets."""


def compute_batch_effects(
    adata: AnnData,
    bio_key: Optional[str] = None,
    batch_key: str = "infer",
    treatment_key: Optional[str] = None,
    control: str = "DMSO",
) -> Tuple[np.array, np.array]:
    """
    Compute batch effects using scone's method of deconvoluting technical from biological effects
     through linear modeling. If no biological differences are expected (e.g. because data is all from
     the same cell line), set `bio_key` to None.

     For details on scone, see Cole et al., 2019: https://doi.org/10.1016/j.cels.2019.03.010

     Parameters
     ----------
     adata : AnnData
             Annotated data matrix

     bio_key : _type_, optional
             Name of column used to delineate biological entities, e.g. cell lines. Default: None

     batch_key : str, optional
             Name of column used to delineate batch effects, e.g. plates. Will try to guess if
             no argument is given. Default: "infer"

    treatment_key: str, optional
            Name of column used to delinate treatments. This is used when computing batch effects across drug-treated plates.
            In that case, we compute batch effects only on untreated cells and then apply the correction factors to all cells.
            If using, please also see `control`.

    control: str, optional
            Name of control treatment. Must be valid value in `treatment_key`.

     Returns
     -------
     Tuple[np.array, np.array]:
             Betas, and gammas of Eq. 3 in the paper. In other words, betas represent the biological effects,
             i.e. how much each feature varied because of biological differences. Gammas represent the technical
             effects, i.e. batch effects.
    """
    from formulaic import Formula
    from patsy.contrasts import Treatment

    withbio = bio_key is not None

    if batch_key == "infer":
        batch_key = _infer_names("batch", adata.obs.columns)[0]

    # combine keys for batch and bio
    joint_keys = [i for i in [bio_key, batch_key] if i is not None]

    if treatment_key is not None:
        adata = adata[adata.obs[treatment_key == control], :]

    data = grouped_op(
        adata, joint_keys, "mean"
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

    if withbio:
        beta_ind = np.where(dmatrix.columns.str.contains(bio_key))[0]
        bio_groups = groups[bio_key].drop_duplicates().to_list()
        contrast_beta = Treatment(reference=0).code_without_intercept(bio_groups)
        betas = (contrast_beta.matrix @ params[beta_ind, :]).T
    else:
        betas = []

    batch_groups = groups[batch_key].drop_duplicates().to_list()
    contrast_gamma = Treatment(reference=0).code_without_intercept(batch_groups)
    gammas = (contrast_gamma.matrix @ params[gamma_ind, :]).T

    # convert to DataFrame's to store metadata
    gammas_df = pd.DataFrame(gammas, columns=batch_groups, index=data.index)
    betas_df = (
        pd.DataFrame(betas, columns=bio_groups, index=data.index)
        if len(betas) != 0
        else []
    )
    return betas_df, gammas_df


def remove_batch_effects(
    adata: AnnData,
    bio_key: Optional[str] = None,
    batch_key: str = "infer",
    copy: Optional[bool] = False,
) -> AnnData:
    """
    Remove batch effects using scone's method of deconvoluting technical from biological effects
    through linear modeling. Note that this preserves biological differences and only
    removes technical effects. If no biological differences are expected (e.g. because data
    is all from the same cell line), set `bio_key` to None.

    For details on scone, see Cole et al., 2019: https://doi.org/10.1016/j.cels.2019.03.010

    Parameters
    ----------
    adata : AnnData
            Annotated data matrix

    bio_key : _type_, optional
            Name of column used to delineate biological entities, e.g. cell lines. Default: None

    batch_key : str, optional
            Name of column used to delineate batch effects, e.g. plates. Will try to guess if
            no argument is given. Default: "infer"

    copy: bool, optional
            If False, will perform operation in-place, else return a modified copy of the data.

    Returns
    -------
    Batch-corrected annotated data matrix. Will only return object if copy=True.
    """
    _, gammas = compute_batch_effects(adata, bio_key, batch_key)
    if copy:
        adata = adata.copy()
    group_obs_fun_inplace(
        adata, batch_key, lambda x, group: x - gammas[group].to_numpy()
    )
    if copy:
        return adata
