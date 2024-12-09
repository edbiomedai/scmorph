import numpy as np
import pandas as pd
from anndata import AnnData

from scmorph.logging import get_logger
from scmorph.pp import drop_na, pca, scale
from scmorph.utils import _get_group_keys, _infer_names, get_grouped_op


def _split_adata_control_drugs(
    adata: AnnData, treatment_key: str, control: str, group_key: str | None = None
) -> tuple[AnnData, AnnData, list[str]]:
    """Split adata into control and drugs"""
    group_keys, treatment_col = _get_group_keys(adata, treatment_key, group_key)

    adata_control = adata[(adata.obs[treatment_col] == control).to_numpy(), :]
    adata_drugs = adata[~(adata.obs.index.isin(adata_control.obs.index)), :]

    return adata_control, adata_drugs, group_keys


def _pca_aggregate(adata: AnnData, cum_var_explained: float = 0.9) -> tuple[AnnData, np.ndarray]:
    scale(adata)
    pca(adata)

    weights = adata.uns["pca"]["variance_ratio"]

    pc_cutoff = np.where(np.cumsum(weights) > cum_var_explained)[0]

    # check if sum of variance explained is more than cum_var_explained
    # and if not just take all PCs
    pc_cutoff = pc_cutoff[0] if pc_cutoff.size > 0 else len(weights)

    weights = weights[:pc_cutoff]
    adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :pc_cutoff]

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

    dists = np.apply_along_axis(lambda x: mahalanobis(control_centroid, x, cov_inv), axis=1, arr=drug_data)

    # add back information about drugs, then collapse drugs with multiple measurements
    dists = pd.Series(dists, index=pd.Series(drug_adata.obs[treatment_col]), name="mahalanobis").groupby(level=0).mean()

    return dists


def aggregate(
    adata: AnnData,
    well_key: str = "infer",
    group_keys: str | list[str] | None = None,
    method: str = "median",
    progress: bool = True,
) -> AnnData:
    """
    Aggregate single-cell measurements into well-level profiles

    Parameters
    ----------
    adata
        Annotated data matrix
    well_key
        Name of column in metadata used to define wells.
    group_keys
        Other column names to group by, e.g. plate names
    method
        Which aggregation to perform. Must be one of 'mean', 'median', 'std',
        'var', 'sem', 'mad', and 'mad_scaled' (i.e. median/mad)
    progress
        Whether to show a progress bar

    Note
    ---------
    If this function produces warnings about dividing by zero, this means that at least
    one group had a median absolute deviation of 0 for a feature. This means that this
    feature is constant in that group. However, this will produce missing values.
    Before proceeding, you should therefore use
    :func:`~.pp.drop_na` with `feature_threshold=1` and `cell_threshold=0`
    to remove features with missing values.

    Returns
    -------
    dists
            Aggregated annotated data matrix
    """
    if well_key == "infer":
        well_key = _infer_names("well", adata.obs.columns)[0]

    if not isinstance(group_keys, list):
        group_keys = [group_keys] if group_keys is not None else []

    group_keys = [well_key, *group_keys]

    return get_grouped_op(adata, group_keys, operation=method, as_anndata=True, progress=progress)


def aggregate_mahalanobis(
    adata: AnnData,
    treatment_key: str = "infer",
    control: str = "DMSO",
    well_key: str = "infer",
    per_treatment: bool = False,
    cov_include_treatment: bool = False,
    cov_from_single_cell: bool = False,
    progress: bool = False,
) -> pd.DataFrame:
    """
    Measure distance between groups using mahalanobis distance

    Parameters
    ----------
    adata
            Annotated data matrix

    treatment_key
            Name of column in metadata used to define treatments

    control
            Name of control treatment. Must be valid value in `treatment_key`.

    well_key
            Name of column in metadata used to define wells. This is needed
            to define the covariance matrix for Mahalanobis distance.

    per_treatment
            Whether to compute PCA and Mahalanobis distance for each treatment separately.

    cov_include_treatment
            Whether to compute covariance matrix from control alone (False) or control and treatment together (True).
            If True, covariance matrices are combined through a weighted sum, where weights represent the number of
            replicates for this drug.

    cov_from_single_cell
            Whether to compute covariance matrix from single cells. This computes distances directly on features
            with no prior PCA. As a result, cov_include_treatment and per_treatment will be ignored (both False).

    progress
            Whether to show a progress bar

    Returns
    -------
    dists
            Mahalanobis distances between treatments
    """
    import anndata
    from tqdm import tqdm

    group_keys, treatment_col = _get_group_keys(adata, treatment_key, well_key)
    treatment_col = treatment_col[0]

    # aggregate
    agg_adata = get_grouped_op(adata, group_keys, "median", progress=progress, as_anndata=True, store=False)

    # compute dists on PCs
    if not per_treatment and not cov_from_single_cell:
        return _pca_mahalanobis(agg_adata, treatment_col, control)

    adata_control, adata_drugs, _ = _split_adata_control_drugs(agg_adata, treatment_col, control, well_key)

    dists = pd.Series(
        index=adata_drugs.obs[treatment_col].unique(),
        name="mahalanobis",
        dtype=np.float64,
    )
    iterator = tqdm(dists.index) if progress else dists.index

    if cov_from_single_cell:
        from scipy.spatial.distance import mahalanobis

        adata_control_sc, _, _ = _split_adata_control_drugs(adata, treatment_col, control, well_key)
        cov = np.cov(adata_control_sc.X, rowvar=False)
        try:
            vi = np.linalg.inv(cov)
        except np.linalg.LinAlgError as e:
            if "Singular matrix" not in str(e):
                raise e

            logger = get_logger()
            logger.warning(
                f"Covariance matrix estimated from single cells of {control} was not invertible."
                + " This is likely because there are very few cells. Falling back to estimating covariance matrix from aggregate data."
            )

            cov_from_single_cell = False

    if cov_from_single_cell:  # check that covariance matrix was invertible
        control_centroid = np.median(adata_control.X, axis=0)

        for cur_treatment in iterator:
            drug_idx = adata_drugs.obs[treatment_col] == cur_treatment
            if sum(drug_idx) == 1:
                drug_centroid = adata_drugs[drug_idx].X.flatten()
            else:
                drug_centroid = np.median(adata_drugs[drug_idx].X, axis=0).flatten()

            dists[cur_treatment] = mahalanobis(control_centroid, drug_centroid, VI=vi)
        return dists

    # only per_treatment = True and cov_from_single_cell = True remaining
    # (or cov_from_single_cell = True but covariance matrix not invertible)
    if not cov_from_single_cell:
        for cur_treatment in iterator:
            drug_idx = adata_drugs.obs[treatment_col] == cur_treatment
            joint_adata = anndata.concat([adata_control, adata_drugs[drug_idx]])
            dists[cur_treatment] = _pca_mahalanobis(joint_adata, treatment_col, control, cov_include_treatment)[0]

    return dists


def aggregate_pc(
    adata: AnnData,
    treatment_key: str = "infer",
    control: str = "DMSO",
    cum_var_explained: float = 0.9,
    progress: bool = True,
) -> pd.Series:
    """
    Measure distance between groups using principle components weighted by variance explained

    Parameters
    ----------
    adata
            Annotated data matrix

    treatment_key
            Name of column in metadata used to define treatments

    control
            Name of control treatment. Must be valid value in `treatment_key`.

    cum_var_explained
            This allows thresholding how many PCs to use during computation of distances.
            It will select the first n PCs until at least this sum of variance has been explained.
            Must be a value between 0 and 1.

    progress
            Whether to show a progress bar

    Returns
    -------
    dists
            Weighted principal component distances to control
    """
    group_keys, treatment_col = _get_group_keys(adata, treatment_key, None)

    agg_adata = get_grouped_op(adata, group_keys, "median", progress=progress, as_anndata=True)
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


# TODO: check speedup options
def aggregate_ttest(
    adata: AnnData,
    treatment_key: str = "infer",
    control: str = "DMSO",
    group_key: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Measure per-feature distance between groups using t-statistics.

    Can be aggregated to a single distance by using :func:`scmorph.pp.tstat_distance`

    Parameters
    ----------
    adata
            Annotated data matrix

    treatment_key
            Name of column in metadata used to define treatments

    control
            Name of control treatment. Must be valid value in `treatment_key`.

    group_key
            Name of column in metadata used to define groups

    Returns
    -------
    dists
            T-statistics between groups

    qvals: :class:`~pandas.DataFrame`
            q-values (i.e. FDR-corrected p-values)
    """
    import scipy
    from statsmodels.stats.multitest import fdrcorrection

    logger = get_logger()

    adata_control, adata_drugs, group_keys = _split_adata_control_drugs(adata, treatment_key, control, group_key)

    tstats = {}
    pvals = {}

    def _get_stats(control: np.array, drug: np.array) -> tuple[np.array, np.array]:
        tstat, pval = scipy.stats.ttest_ind(control, drug, axis=0, equal_var=False)
        return tstat, pval

    for group, idx in adata_drugs.obs.groupby(group_keys, observed=True).groups.items():
        cur_drug = adata_drugs[idx, :]
        tstat, pval = _get_stats(adata_control.X, cur_drug.X)
        pvals[group] = pval
        tstats[group] = tstat

    pvalsdf = pd.DataFrame(pvals, index=adata.var.index).T

    # When feature is similar between treatment and control
    # the p-value is nan. When doing FDR, this forces all pvals
    # to nan. We replace these with 1.
    if np.isnan(pvalsdf.values).any():
        logger.warning(
            "Some p-values are NaN. This is likely due to constant features between control and treatment."
            + " These will be replaced with 1."
        )
        pvalsdf.fillna(1, inplace=True)

    qvalsdf = pd.DataFrame(
        fdrcorrection(np.ravel(pvalsdf))[1].reshape(pvalsdf.shape),
        columns=pvalsdf.columns,
        index=pvalsdf.index,
    )

    tstatsdf = pd.DataFrame(tstats, index=adata.var.index).T

    # verify that all rows and columns are in correct order
    qvalsdf = qvalsdf.reindex_like(tstatsdf, copy=False)

    return tstatsdf.T, qvalsdf.T


def tstat_distance(tstats: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize t-statistics into per group. See :func:`scmorph.pp.aggregate_ttest` for details.

    Parameters
    ----------
    tstats
        t-statistics computed with :func:`scmorph.pp.aggregate_ttest`

    Returns
    -------
    dists
        Per-group t-statistic distances
    """
    # score[j] = sqrt(t_1^2 + ... + t_i^2)
    # where i = features and j = compounds
    return tstats.pow(2).sum(axis=0).pow(0.5)
