
#' @title run_slingshot
#'
#' @description This function wraps conveniently wraps Slingshot and returns trajectories
#' @param X_pca Matrix of PC coordinates. Caution: rows should be cells, columns features.
#' @param clusterLabels Vector of assigned cluster labels. Must have as many dimensions as columns in `X_pca`.
#' @param start.clus Integer of assigned start cluster, if any.
#' @param end.clus Vector of integers of assigned terminal cluster(s), if any.
#' @param ... Other arguments passed to `slingshot`
#' @return list with curve coordinates per lineage and pseudotime per cell
#'
#' @import slingshot
#' @import SingleCellExperiment
#' @author Jesko Wagner
#' @export
run_slingshot <- function(X_pca, clusterLabels, start.clus = NULL, end.clus = NULL, ...){
    library(slingshot)
    library(SingleCellExperiment)
    sce = SingleCellExperiment(assays=list(matrix(ncol=nrow(X_pca))), # spoof count data
                               reducedDims=list("PCA"=X_pca))
    sce = slingshot(sce,
                    clusterLabels=clusterLabels,
                    reducedDim = "PCA",
                    start.clus=start.clus,
                    end.clus = end.clus,
                    ...)
   traj_results = SlingshotDataSet(sce)
   curve_coords = lapply(traj_results@curves, function(x) x[["s"]])
   curve_coords = lapply(curve_coords, function(x) {colnames(x) <- paste0("PC", 1:ncol(x)); x})
   return(list(curve_coords = curve_coords, pseudotime = sce$slingPseudotime_1))
}
