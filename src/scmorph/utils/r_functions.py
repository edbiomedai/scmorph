from collections.abc import Callable


def _initialize_r_functions() -> None:
    """
    Trajectory inference R functions

    Create R functions from raw text to avoid the complexity of shipping R files
    in a Python package.
    """
    import rpy2.robjects as ro

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
   pseudotime = slingPseudotime(slingshot_object, na=FALSE)
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

test_differential_progression <- function(cellWeights, pseudotime, conditions, global=TRUE, lineages=TRUE, ...){
    suppressPackageStartupMessages(library(condiments))
    progressionTest(cellWeights=as.matrix(cellWeights),
                    pseudotime=as.matrix(pseudotime),
                    conditions=conditions,
                    global=global,
                    lineages=lineages,
                    ...)
}

test_differential_differentiation <- function(cellWeights,
                                              conditions,
                                              global=TRUE,
                                              pairwise=TRUE, ...){
    suppressPackageStartupMessages(library(condiments))
    differentiationTest(cellWeights=as.matrix(cellWeights),
                        conditions=conditions,
                        global=global,
                        pairwise=pairwise,
                        ...)
}
"""
    ro.r(r_functions)


def _clean_R_env(except_obj: str | list[str] | None = None) -> None:
    import rpy2.robjects as ro
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
    rm_fun = ro.r(rm_fun)
    with localconverter(default_converter + none_cv):
        rm_fun(except_obj)  # type: ignore


def _load_R_functions(  # type: ignore
    function: str = "all",
) -> Callable | dict:
    import rpy2.robjects as ro

    _initialize_r_functions()  # initialize functions in R global environment

    if function != "all":
        return ro.r[function]

    run_slingshot = ro.r["run_slingshot"]
    test_common_trajectory = ro.r["test_common_trajectory"]
    test_differential_progression = ro.r["test_differential_progression"]
    test_differential_differentiation = ro.r["test_differential_differentiation"]

    all_funcs = {
        "run_slingshot": run_slingshot,
        "test_common_trajectory": test_common_trajectory,
        "test_differential_progression": test_differential_progression,
        "test_differential_differentiation": test_differential_differentiation,
    }

    return all_funcs


def _None_converter():  # type: ignore
    """Adapted from https://stackoverflow.com/a/65810724 by Mike Krassowski"""
    import rpy2.robjects as ro

    def _none2null(none_obj) -> ro.NULL:  # type: ignore
        return ro.r("NULL")

    none_converter = ro.conversion.Converter("None converter")
    none_converter.py2rpy.register(type(None), _none2null)
    return none_converter
