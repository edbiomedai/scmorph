import numpy as np
import pandas as pd
from anndata import AnnData

from scmorph.utils import _clean_R_env, _load_R_functions, _None_converter


def slingshot(
    adata: AnnData,
    cluster_labels: str = "leiden",
    start_clus: int | None = None,
    end_clus: int | None = None,
    n_comps: int = 10,
) -> None:
    """
    Trajectory inference using Slingshot

    Parameters
    ----------
    adata
        AnnData object
    cluster_labels
        Column name in `obs` defining clusters
    start_clus
        Start cluster label
    end_clus
        End cluster label[s]
    n_comps
        Number of principal components to use for trajectory inference. Default 10

    Returns
    -------
    AnnData object is modified in-place with trajectory information added to `.obsm` and `.uns`
    """
    import rpy2.robjects as ro
    from rpy2.robjects import default_converter, pandas2ri
    from rpy2.robjects.conversion import localconverter

    from scmorph.logging import get_logger

    log = get_logger()
    r_slingshot = _load_R_functions("run_slingshot")
    none_cv = _None_converter()

    clusters = adata.obs[cluster_labels]

    with localconverter(default_converter + pandas2ri.converter):
        X_pca = adata.obsm["X_pca"][:, :n_comps]
        X_pca = ro.conversion.py2rpy(pd.DataFrame(X_pca))
        cluster = ro.conversion.py2rpy(clusters)

    log.info("Running Slingshot...")
    with localconverter(default_converter + none_cv):
        r_res = r_slingshot(X_pca, cluster, start_clus, end_clus)

    adata.uns["slingshot_object"] = r_res[0]
    adata.uns["slingshot_curve_coords"] = np.array(r_res[1])
    adata.obsm["slingshot_pseudotime"] = np.array(r_res[2])
    adata.obsm["slingshot_cell_assignments"] = np.array(r_res[3])
    _clean_R_env(except_obj="slingshot_object")


def _test_condiments(
    fun: str,
    adata: AnnData,
    conditions: pd.Series | np.ndarray,
    all_pairs: bool = True,
    pairwise: bool = True,
    parallel: bool = True,
    lineages: bool = True,
) -> None:
    import rpy2.robjects as ro
    from rpy2.robjects import default_converter, pandas2ri
    from rpy2.robjects.conversion import localconverter

    from scmorph import get_logger

    log = get_logger()

    if "slingshot_object" not in adata.uns.keys():
        log.error("Slingshot assignments not found in `.uns`. Please run :func:`scmorph.tl.slingshot` first.")
        raise KeyError("adata.uns['slingshot_object'] not found.")

    r_test_fun = _load_R_functions(fun)

    conditions = adata.obs[conditions]

    with localconverter(default_converter + pandas2ri.converter):
        conditions = ro.conversion.py2rpy(conditions)
        if fun != "test_common_trajectory":
            weights = ro.conversion.py2rpy(pd.DataFrame(adata.obsm["slingshot_cell_assignments"]))

    if fun == "test_common_trajectory":
        r_res = r_test_fun(adata.uns["slingshot_object"], conditions, parallel=parallel)
    elif fun == "test_differential_differentiation":
        r_res = r_test_fun(weights, conditions, all_pairs, pairwise=pairwise)
    elif fun == "test_differential_progression":
        with localconverter(default_converter + pandas2ri.converter):
            pseudotime = ro.conversion.py2rpy(pd.DataFrame(adata.obsm["slingshot_pseudotime"]))
        r_res = r_test_fun(weights, pseudotime, conditions, all_pairs, lineages=lineages)
    else:
        log.error("Function must be one of `test_differential_differentiation` or `test_differential_progression`")
        r_res = None

    # convert to pandas
    with localconverter(default_converter + pandas2ri.converter):
        res = ro.conversion.rpy2py(r_res)

    print(res)
    adata.uns[fun] = res
    _clean_R_env()


def test_common_trajectory(adata: AnnData, conditions: pd.Series | np.ndarray, parallel: bool = True) -> None:
    """
    Test for common trajectory using condiments' `topologyTest`

    Parameters
    ----------
    adata
        AnnData object
    conditions
        Column name in `obs` defining conditions
    parallel
        Use parallel processing.

    Returns
    -------
    AnnData object is modified in-place with common trajectory test results added to `.uns`
    """
    import rpy2

    from scmorph.logging import get_logger

    log = get_logger()

    # compute p values
    log.info("Testing for common trajectory...")
    try:
        _test_condiments("test_common_trajectory", adata, conditions, parallel=parallel)
    except rpy2.rinterface_lib.embedded.RRuntimeError as e:
        err = "Error in topologyTest::testCommonTrajectory(slingshot_object, conditions,"
        err += f"parallel=parallel): \n{e}\n\n"
        err += "This is often caused because the clusters were too granular. Try clustering with fewer clusters, rerun slingshot and try again."
        log.error(err)


def test_differential_progression(
    adata: AnnData,
    conditions: pd.Series | np.ndarray,
    all_pairs: bool = True,
    lineages: bool = True,
) -> None:
    """
    Test for differential progression using condiments' `progressionTest`

    Parameters
    ----------
    adata
        AnnData object
    conditions
        Column name in `obs` defining conditions
    all_pairs
        Test all pairs of conditions.
    lineages
        Test all lineages independently.

    Returns
    -------
    AnnData object is modified in-place with differential progression test results added to `.uns`
    """
    from scmorph import get_logger

    log = get_logger()
    log.info("Testing differential progression...")
    _test_condiments(
        "test_differential_progression",
        adata,
        conditions,
        all_pairs=all_pairs,
        lineages=lineages,
    )


def test_differential_differentiation(
    adata: AnnData,
    conditions: pd.Series | np.ndarray,
    all_pairs: bool = True,
    pairwise: bool = True,
) -> None:
    """
    Test for differential differentiation using condiments' `differentiationTest`

    Parameters
    ----------
    adata
        AnnData object
    conditions
        Column name in `obs` defining conditions
    all_pairs
        Test all pairs of conditions.
    pairwise
        Test all pairs independently.

    Returns
    -------
    AnnData object is modified in-place with differential differentiation test results added to `.uns`
    """
    from scmorph import get_logger

    log = get_logger()
    log.info("Testing differential differentiation...")
    _test_condiments(
        "test_differential_differentiation",
        adata,
        conditions,
        all_pairs=all_pairs,
        pairwise=pairwise,
    )
