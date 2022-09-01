from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

from scmorph._utils import get_group_keys, get_grouped_op
from scmorph.logging import get_logger
from scmorph.pp import drop_na, pca, scale


def _split_adata_control_drugs(
    adata: AnnData, treatment_key: str, control: str, group_key: Optional[str] = None
) -> Tuple[AnnData, AnnData, List[str]]:
    """
    Split adata into control and drugs
    """
    group_keys, treatment_col = get_group_keys(adata, treatment_key, group_key)

    adata_control = adata[(adata.obs[treatment_col] == control).to_numpy(), :]
    adata_drugs = adata[~(adata.obs.index.isin(adata_control.obs.index)), :]

    return adata_control, adata_drugs, group_keys


def _pca_aggregate(
    adata: AnnData, cum_var_explained: float = 0.9
) -> Tuple[AnnData, np.array]:
    scale(adata)
    pca(adata)

    weights = adata.uns["pca"]["variance_ratio"]

    pc_cutoff = np.where(np.cumsum(weights) > cum_var_explained)[0]
    if pc_cutoff.size > 0:  # sum of variance explained is more than cum_var_explained
        weights = weights[: pc_cutoff[0]]
        adata.obsm["X_pca"] = adata.obsm["X_pca"][:, : pc_cutoff[0]]

    return adata, weights


def _pca_mahalanobis(
    joint_adata: AnnData,
    treatment_col: str,
    control: str,
    cov_include_treatment: bool = False,
) -> pd.Series:
    from scipy.spatial.distance import mahalanobis

    logger = get_logger()
    drop_na(joint_adata, feature_threshold=0, cell_threshold=1)  # drop NA columns
    joint_adata, _ = _pca_aggregate(joint_adata)

    control_idx = (joint_adata.obs[treatment_col] == control).to_numpy().squeeze()
    control_data = joint_adata.obsm["X_pca"][control_idx, :]

    control_centroid = control_data.mean(axis=0)

    drug_adata = joint_adata[np.invert(control_idx), :]
    drug_data = drug_adata.obsm["X_pca"]

    if control_data.shape[1] > 1:
        cov = np.cov(control_data, rowvar=False)

        if cov_include_treatment:
            # make a second covariance matrix for the treatment data
            # weight each matrix by number of wells in the group
            # then combine
            if drug_data.shape[0] < 2:
                logger.warning(
                    "Not enough drug replicates to compute covariance."
                    + " Using control covariance. Use cov_include_treatment=False"
                    + " to avoid computing covariance on treatments when not all drugs have replicates."
                )
            else:
                drug_cov = np.cov(drug_data, rowvar=False)

                # weigh
                drug_cov /= drug_data.shape[0]
                cov /= control_data.shape[0]

                # combine
                cov += drug_cov

        cov_inv = np.linalg.inv(cov)
    else:
        cov_inv = np.array([1])

    dists = np.apply_along_axis(
        lambda x: mahalanobis(control_centroid, x, cov_inv), axis=1, arr=drug_data
    )

    # add back information about drugs, then collapse drugs with multiple measurements
    dists = (
        pd.Series(
            dists, index=pd.Series(drug_adata.obs[treatment_col[0]]), name="mahalanobis"
        )
        .groupby(level=0)
        .mean()
    )

    return dists


def aggregate_mahalanobis(
    adata: AnnData,
    treatment_key: str = "infer",
    control: str = "DMSO",
    well_key: str = "infer",
    per_treatment: bool = False,
    cov_include_treatment: bool = False,
    cov_from_single_cell: bool = False,
) -> pd.DataFrame:
    """
    Aggregate single-cell matrix by treatment using Mahalanobis distance

    Parameters
    ----------
    adata : AnnData
            The AnnData object as read in with annmorph

    treatment_key : str
            Name of column in metadata used to define treatments

    control : str
            Name of control treatment. Must be valid value in `treatment_key`.

    well_key : str
            Name of column in metadata used to define wells. This is needed
            to define the covariance matrix for Mahalanobis distance.

    per_treatment : bool
            Whether to compute PCA and Mahalanobis distance for each treatment separately.

    cov_include_treatment : bool
            Whether to compute covariance matrix from control alone (False) or control and treatment together (True).
            If True, covariance matrices are combined through a weighted sum, where weights represent the number of
            replicates for this drug.

    cov_from_single_cell : bool
            Whether to compute covariance matrix from single cells. This computes distances directly on features
            with no prior PCA. As a result, cov_include_treatment and per_treatment will be ignore (both False).

    Returns
    ----------
    DataFrame, containing mahalanobis distances between treatments
    """
    import anndata

    group_keys, treatment_col = get_group_keys(adata, treatment_key, well_key)

    # aggregate
    agg_data = get_grouped_op(adata, group_keys, "median")
    X = agg_data.T
    # TODO: the multiindex only works if well_key was passed
    meta_cols = pd.MultiIndex.from_tuples([*X.index]).to_frame().reset_index(drop=True)
    meta_cols.columns = group_keys
    meta_cols.index = meta_cols.index.astype(str)  # avoid conversion warnings
    agg_adata = AnnData(X=X.to_numpy(), obs=meta_cols)

    # compute dists on PCs
    if not per_treatment and not cov_from_single_cell:
        dists = _pca_mahalanobis(agg_adata, treatment_col[0], control)
        return dists

    adata_control, adata_drugs, _ = _split_adata_control_drugs(
        agg_adata, treatment_col[0], control, well_key
    )

    dists = pd.Series(
        index=adata_drugs.obs[treatment_col[0]].unique(),
        name="mahalanobis",
        dtype=np.float64,
    )
    if cov_from_single_cell:
        from scipy.spatial.distance import mahalanobis

        adata_control_sc, _, _ = _split_adata_control_drugs(
            adata, treatment_col[0], control, well_key
        )
        cov = np.cov(adata_control_sc.X, rowvar=False)
        try:
            vi = np.linalg.inv(cov)
        except np.linalg.LinAlgError as e:
            if "Singular matrix" in str(e):
                logger = get_logger()
                logger.warning(
                    f"Covariance matrix estimated from single cells of {control} was not invertible."
                    + "This is likely because there are very few cells."
                    + " Falling back to estimating covariance matrix from aggregate data."
                )
                cov_from_single_cell = False
            else:
                raise

        if cov_from_single_cell:  # check that covariance matrix was invertible
            control_centroid = np.median(adata_control.X, axis=0)
            for cur_treatment in adata_drugs.obs[treatment_col[0]].unique():
                drug_idx = adata_drugs.obs[treatment_col[0]] == cur_treatment
                if sum(drug_idx) == 1:
                    drug_centroid = adata_drugs[drug_idx].X
                else:
                    drug_centroid = np.median(adata_drugs[drug_idx].X, axis=0).flatten()

                dists[cur_treatment] = mahalanobis(
                    control_centroid, drug_centroid, VI=vi
                )

    if not cov_from_single_cell:
        for cur_treatment in adata_drugs.obs[treatment_col[0]].unique():
            drug_idx = adata_drugs.obs[treatment_col[0]] == cur_treatment

            joint_adata = anndata.concat([adata_control, adata_drugs[drug_idx]])

            dists[cur_treatment] = _pca_mahalanobis(
                joint_adata, treatment_col[0], control, cov_include_treatment
            )

    return dists


def aggregate_pc(
    adata: AnnData,
    treatment_key: str = "infer",
    control: str = "DMSO",
    cum_var_explained: float = 0.9,
) -> pd.Series:
    """
    Aggregate single-cell matrix by treatment using weighted principal component distances

    Parameters
    ----------
    adata : AnnData
            The AnnData object as read in with annmorph

    treatment_key : str
            Name of column in metadata used to define treatments

    control : str
            Name of control treatment. Must be valid value in `treatment_key`.

    cum_var_explained : float
            This allows thresholding how many PCs to use during computation of distances.
            It will select the first n PCs until at least this sum of variance has been explained.
            Must be a value between 0 and 1.

    Returns
    ----------
    One Series, containing weighted principal component distances to control
    """

    group_keys, treatment_col = get_group_keys(adata, treatment_key, None)

    agg_data = get_grouped_op(adata, group_keys, "median")
    X = agg_data.T
    meta_cols = pd.DataFrame(X.index, columns=group_keys)
    meta_cols.index = meta_cols.index.astype(str)  # avoid conversion warnings
    agg_adata = AnnData(X=X.to_numpy(), obs=meta_cols)
    drop_na(agg_adata, 0, 1)  # drop NA columns
    agg_adata, weights = _pca_aggregate(agg_adata, cum_var_explained)

    # determine reference point
    agg_control = agg_adata[(agg_adata.obs[treatment_col] == control).to_numpy(), :]
    control_centroid = agg_control.obsm["X_pca"]

    # compute euclidean distances
    # square root of (Var_PC1 x (drug1_PC1 - DMSO_PC1)^2 + Var_PC2 x (drug1_PC2 - DMSO_PC2)^2 + … + Var_PCk x (drug1_PCk - DMSO_PCk)^2 )
    dist = np.sqrt(
        np.sum(
            weights * np.square(control_centroid - agg_adata.obsm["X_pca"]),
            axis=1,
        )
    )

    return pd.Series(dist, index=agg_adata.obs[treatment_col[0]], name="pc_dist")


# this implementation is slower than the below, but guaranteed to give the right results
# we noticed differences with the alternative implementation that still have to be addressed
# TODO: check differences between implementations
def aggregate_ttest(
    adata: AnnData,
    treatment_key: str = "infer",
    control: str = "DMSO",
    group_key: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate single-cell matrix by treatment using t-statistics

    Parameters
    ----------
    adata : AnnData
            The AnnData object as read in with annmorph

    treatment_key : str
            Name of column in metadata used to define treatments

    control : str
            Name of control treatment. Must be valid value in `treatment_key`.

    group_key : str
            Name of column in metadata used to define groups

    Returns
    ----------
    Two pd.DataFrames, one with the t-statistics and one with the q-values (i.e. FDR-corrected p-values)
    """
    import scipy
    from statsmodels.stats.multitest import fdrcorrection

    adata_control, adata_drugs, group_keys = _split_adata_control_drugs(
        adata, treatment_key, control, group_key
    )

    tstats = {}
    pvals = {}

    def _get_stats(control: np.array, drug: np.array) -> Tuple[np.array, np.array]:
        tstat, pval = scipy.stats.ttest_ind(control, drug, axis=0, equal_var=False)
        return tstat, pval

    for group, idx in adata_drugs.obs.groupby(group_keys).groups.items():
        cur_drug = adata_drugs[idx, :]
        tstat, pval = _get_stats(adata_control.X, cur_drug.X)
        pvals[group] = pval
        tstats[group] = tstat

    pvalsdf = pd.DataFrame(pvals, index=adata.var.index).T
    qvalsdf = pd.DataFrame(
        fdrcorrection(np.ravel(pvalsdf))[1].reshape(pvalsdf.shape),
        columns=pvalsdf.columns,
        index=pvalsdf.index,
    )

    tstatsdf = pd.DataFrame(tstats, index=adata.var.index).T

    # verify that all rows and columns are in correct order
    qvalsdf = qvalsdf.reindex_like(tstatsdf, copy=False)

    return tstatsdf.T, qvalsdf.T  # transpose to match new implementation


def tstat_distance(tstats: pd.DataFrame) -> pd.DataFrame:
    """
    Compute t-static distances

    Parameters
    ----------
    tstats : pd.DataFrame
        t-statistics computed with `aggregate_test`

    Returns
    -------
    DataFrame
        Per-group t-statistic distances
    """
    # score[j] = sqrt(t_1^2 + ... + t_i^2)
    # where i = features and j = compounds
    return tstats.pow(2).sum(axis=0).pow(0.5)


# keep for now until differences have been found
# def aggregate_ttest(
#     adata: AnnData,
#     treatment_key: str = "infer",
#     control: str = "DMSO",
# ) -> Tuple[pd.DataFrame, pd.DataFrame]:
#     """
#     Aggregate single-cell matrix by treatment using t-statistics

#     Parameters
#     ----------
#     adata : AnnData
#             The AnnData object as read in with annmorph

#     treatment_key : str
#             Name of column in metadata used to define treatments

#     control : str
#             Name of control treatment

#     Returns
#     ----------
#     Two DataFrames, one with the t-statistics and one with the q-values (i.e. FDR-corrected p-values)
#     """
#     import scipy
#     from statsmodels.stats.multitest import fdrcorrection

#     # inferring treatment names if necessary
#     if isinstance(treatment_key, str) and treatment_key == "infer":
#         treatment_col = _infer_names("treatment", adata.obs.columns)
#     else:
#         treatment_col = [treatment_key]
#     # end inference

#     means = grouped_op(adata, treatment_col, "mean")
#     vars = grouped_op(adata, treatment_col, "var")
#     Ns = adata.obs[treatment_col].squeeze().value_counts()

#     groups = pd.Index(means.columns, name=treatment_col[0])

#     means.columns = groups
#     vars.columns = groups

#     u_control = means.loc[:, control]
#     u_drugs = means.loc[:, means.columns != control]

#     var_control = vars.loc[:, control]
#     var_drugs = vars.loc[:, vars.columns != control]

#     N_control = Ns.loc[control]
#     N_drugs = Ns.loc[Ns.index != control]

#     ### This section was adapted from scipy.stats._mstats_basic.ttest_ind
#     # This adaption was implemented to compute all t-statistics and df's at once,
#     # i.e. to avoid computing variance and mean of control multiple times
#     # as would happen when using scipy.stats.ttest_ind
#     # it uses pandas' vectorized implementations of operators to efficiently
#     # compute the statistics

#     vn1 = var_control / N_control
#     vn2 = var_drugs / N_drugs
#     with np.errstate(divide="ignore", invalid="ignore"):
#         df_nom = (vn2.add(vn1, axis="rows")).pow(2)
#         df_control_denom = (vn1.pow(2)).div(N_control - 1)
#         df_drug_denom = (vn2.pow(2)).div(N_drugs - 1)
#         df = df_nom / (df_drug_denom.add(df_control_denom, axis="rows"))
#         # df = (vn1 + vn2)**2 / (vn1**2 / (n1 - 1) + vn2**2 / (n2 - 1))

#     # If df is undefined, variances are zero.
#     # It doesn't matter what df is as long as it is not NaN.
#     df = np.where(np.isnan(df), 1, df)
#     denom = vn2.add(vn1, axis="rows").pow(0.5)

#     with np.errstate(divide="ignore", invalid="ignore"):
#         t = (-u_drugs).add(u_control, axis="rows") / denom
#         # (u_control - u_drugs) / denom, but vectorized

#     t, pval = scipy.stats._stats_py._ttest_finish(df, t, "two-sided")
#     ### End of adapted section

#     qvals = pd.DataFrame(
#         fdrcorrection(np.ravel(pval))[1].reshape(pval.shape),
#         columns=pval.columns,
#         index=pval.index,
#     )
#     return t, qvals
