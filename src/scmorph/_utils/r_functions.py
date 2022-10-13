from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from anndata import AnnData


def _initialize_r_functions() -> None:
    """
    Create R functions from raw text to avoid the complexity of shipping R files
    in a Python package. The R code will be provided in a separate repository
    for better readability.
    """
    from rpy2 import robjects

    r_functions = """
run_slingshot <- function(X_pca, clusterLabels, start.clus = NULL, end.clus = NULL, ...){
    suppressPackageStartupMessages(library(slingshot))
    colnames(X_pca) = paste0("PC", 1:ncol(X_pca))
    slingshot_object = slingshot(X_pca,
                    clusterLabels=clusterLabels,
                    start.clus=start.clus,
                    end.clus = end.clus,
                    ...)
   curve_coords = slingCurves(slingshot_object, as.df = TRUE)
   pseudotime = slingPseudotime(slingshot_object)
   cell_assignments = slingCurveWeights(slingshot_object, as.probs=TRUE)
   return(list(slingshot_object = slingshot_object,
               curve_coords = curve_coords,
               pseudotime = pseudotime,
               cell_assingments=cell_assignments))
}

test_common_trajectory <- function(slingshot_object,
                                   conditions,
                                   parallel=TRUE,
                                   ...){
    suppressPackageStartupMessages(library(condiments))
    conditions = as.character(conditions)
    topologyTest(sds = slingshot_object, conditions, parallel=parallel, ...)
}

test_differential_progression <- function(slingshot_object, conditions, global=TRUE, lineages=TRUE, ...){
    suppressPackageStartupMessages(library(condiments))
    progressionTest(slingshot_object,
                    conditions=df$conditions,
                    global=global,
                    lineages=lineages,
                    ...)
}

test_differential_differentiation <- function(slingshot_object,
                                              conditions,
                                              global=TRUE,
                                              pairwise=TRUE, ...){
    suppressPackageStartupMessages(library(condiments))
    differentiationTest(slingshot_object,
                        conditions=df$conditions,
                        global=global,
                        pairwise=pairwise,
                        ...)
}
"""
    robjects.r(r_functions)


def _clean_R_env(except_obj: Optional[Union[str, List[str]]] = None) -> None:
    from rpy2 import robjects
    from rpy2.robjects import default_converter
    from rpy2.robjects.conversion import localconverter

    none_cv = _None_converter()

    if isinstance(except_obj, list):
        except_obj = "|".join(except_obj)

    rm_fun = """
    rm_all <- function(except_obj=NULL){
    obj = ls(all.names=TRUE)
    if (!is.null(except_obj)){
        obj = obj[!grepl(except_obj, obj)]
    }
    rm(list=obj)
    invisible(gc())
    }
    """
    rm_fun = robjects.r(rm_fun)
    with localconverter(default_converter + none_cv):
        rm_fun(except_obj)  # type: ignore


def _load_R_functions(
    function: str = "all",
) -> Union[
    ro.functions.SignatureTranslatedFunction,
    Dict[str, ro.functions.SignatureTranslatedFunction],
]:
    from rpy2 import robjects

    _initialize_r_functions()  # intialize functions in R global environment

    if function != "all":
        return robjects.r[function]

    slingshot = robjects.r["run_slingshot"]
    test_common_trajectory = robjects.r["test_common_trajectory"]
    test_differential_progression = robjects.r["test_differential_progression"]
    test_differential_differentiation = robjects.r["test_differential_differentiation"]

    all_funcs = {
        "slingshot": slingshot,
        "test_common_trajectory": test_common_trajectory,
        "test_differential_progression": test_differential_progression,
        "test_differential_differentiation": test_differential_differentiation,
    }

    return all_funcs


def _None_converter() -> ro.conversion.Converter:
    """Adapted from https://stackoverflow.com/a/65810724"""
    import rpy2.robjects as ro

    def _none2null() -> ro.NULL:
        return ro.r("NULL")

    none_converter = ro.conversion.Converter("None converter")
    none_converter.py2rpy.register(type(None), _none2null)
    return none_converter


def slingshot(
    adata: AnnData,
    cluster_labels: str = "leiden",
    start_clus: Optional[int] = None,
    end_clus: Optional[int] = None,
    n_comps: int = 10,
) -> None:
    """
    Trajectory inference using Slingshot

    Parameters
    ----------
    adata : AnnData
        AnnData object
    cluster_labels : str
        Column name in `obs` defining clusters
    start_clus : Optional[int], optional
        Start cluster label
    end_clus : Optional[int], optional
        End cluster label[s]
    n_comps : int, optional
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
    r_slingshot = _load_R_functions("slingshot")
    none_cv = _None_converter()

    clusters = adata.obs[cluster_labels]

    log.info("Running Slingshot...")
    with localconverter(default_converter + pandas2ri.converter):
        X_pca = adata.obsm["X_pca"][:, :n_comps]
        X_pca = ro.conversion.py2rpy(pd.DataFrame(X_pca))
        cluster = ro.conversion.py2rpy(clusters)

    with localconverter(default_converter + none_cv):
        r_res = r_slingshot(X_pca, cluster, start_clus, end_clus)  # type: ignore

    adata.uns["slingshot_object"] = r_res[0]
    adata.uns["slingshot_curve_coords"] = np.array(r_res[1])
    adata.obsm["slingshot_pseudotime"] = np.array(r_res[2])
    adata.obsm["slingshot_cell_assignments"] = np.array(r_res[3])
    _clean_R_env(except_obj="slingshot_object")


def test_common_trajectory(
    adata: AnnData, conditions: Union[pd.Series, np.array], parallel: bool = True
) -> None:
    """
    Test for common trajectory

    Parameters
    ----------
    adata : AnnData
        AnnData object
    conditions : Union[pd.Series, np.array]
        Column name in `obs` defining conditions
    parallel : bool, optional
        Use parallel processing. Default: True

    Returns
    -------
    AnnData object is modified in-place with common trajectory test results added to `.uns`
    """
    import rpy2.robjects as ro
    from rpy2.robjects import default_converter, pandas2ri
    from rpy2.robjects.conversion import localconverter

    from scmorph.logging import get_logger

    log = get_logger()

    if "slingshot_object" not in adata.uns.keys():
        log.error("Slingshot object not found in `.uns`. Please run `slingshot` first.")
        raise KeyError("adata.uns['slingshot_object'] not found.")
    r_test_common_trajectory = _load_R_functions("test_common_trajectory")
    conditions = adata.obs[conditions]

    with localconverter(default_converter + pandas2ri.converter):
        conditions = ro.conversion.py2rpy(conditions)

    # compute p values
    log.info("Testing for common trajectory...")
    r_res = r_test_common_trajectory(  # type: ignore
        adata.uns["slingshot_object"], conditions, parallel=parallel
    )

    # convert to pandas
    with localconverter(default_converter + pandas2ri.converter):
        res = ro.conversion.rpy2py(r_res)

    print(res)
    adata.uns["slingshot_common_trajectory_test"] = res
