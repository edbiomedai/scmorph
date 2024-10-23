

"""Support functions to run slingshot and condiments

Note: these functions are not actually being called from this file!
They are here for easier reference, but are actually defined in and sourced from
`utils/r_functions.py`. Reason: we do not want to ship R files in a Python package,
as that makes distribution harder. Therefore, take these functions only as reference
that should be in sync with those used in scmorph.
"""

#' @title run_slingshot
#'
#' @description This function conveniently wraps slingshot and returns trajectories and pseudotime
#' @param X_pca Matrix of PC coordinates. Caution: rows should be cells, columns features.
#' @param clusterLabels Vector of assigned cluster labels. Must have as many dimensions as columns in `X_pca`.
#' @param start.clus Integer of assigned start cluster, if any.
#' @param end.clus Vector of integers of assigned terminal cluster(s), if any.
#' @param ... Other arguments passed to `slingshot`
#' @return list with slingshot object, curve coordinates per lineage, pseudotime and lineage probability per cell
#'
#' @import slingshot
#' @author Jesko Wagner
#' @export
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


#' @title test_common_trajectory
#'
#' @description This function wraps `topologyTest` from the condiments package
#' @param slingshot_object slingshot object as returned by `run_slingshot` in slot `slingshot_object`
#' @param conditions Vector with per-cell condition label
#' @param parallel Whether to use multi-core processing
#' @param ... Other parameters passed to `topologyTest`
#' @return data.frame with p-values indicating whether a single trajectory is sufficient
#'
#' @import condiments
#' @author Jesko Wagner
#' @export
#' @note If the p-value reported by this function is < 0.05, consider fitting multiple trajectories.
#' This will currently require manual R implementations following the paper.
test_common_trajectory <- function(slingshot_object,
                                   conditions,
                                   parallel=TRUE,
                                   ...){
    suppressPackageStartupMessages(library(condiments))
    conditions = as.character(conditions)
    topologyTest(sds = slingshot_object, conditions, parallel=parallel, ...)
}


#' @title test_differential_progression
#'
#' @description This function wraps `progressionTest` from the condiments package
#' @param cellWeights the cell weights for each lineage
#' @param pseudotime a matrix of pseudotime values, each row represents a cell and each column represents a lineage.
#' @param conditions Vector with per-cell condition label
#' @param global If TRUE, test for all pairs simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @param ... Other parameters passed to `progressionTest`
#' @return data.frame with p.value per lineage, indicating whether a lineage
#' has differential progression across trajectories
#'
#' @import condiments
#' @author Jesko Wagner
#' @export
test_differential_progression <- function(cellWeights, pseudotime, conditions, global=TRUE, lineages=TRUE, ...){
    suppressPackageStartupMessages(library(condiments))
    progressionTest(as.matrix(cellWeights),
                    pseudotime=as.matrix(pseudotime),
                    conditions=conditions,
                    global=global,
                    lineages=lineages,
                    ...)
}


#' @title test_differential_differentiation
#'
#' @description This function wraps `differentiationTest` from the condiments package
#' @param cellWeights a matrix of cell weights defining the probability that a cell belongs to a particular lineage
#' @param conditions Vector with per-cell condition label
#' @param global If TRUE, test for all pairs simultaneously.
#' @param pairwise If TRUE, test for all pairs independently.
#' @param ... Other parameters passed to `differentiationTest`
#' @return A data frame with 3 columns:
#' \itemize{
#'   \item *pair* for individual pairs, the lineages numbers. For global,
#'   \code{"All"}.
#'   \item *p.value* the pvalue for the test at the global or pair level
#'   \item *statistic* The classifier accuracy
#' }
#' @import condiments
#' @author Jesko Wagner
#' @export
test_differential_differentiation <- function(cellWeights,
                                              conditions,
                                              global=TRUE,
                                              pairwise=TRUE, ...){
    suppressPackageStartupMessages(library(condiments))
    differentiationTest(as.matrix(cellWeights),
                        conditions=conditions,
                        global=global,
                        pairwise=pairwise,
                        ...)
}
