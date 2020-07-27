#' Compute neighborhood purity
#'
#' Use a hypersphere-based approach to compute the \dQuote{purity} of each cluster based on the number of contaminating observations from different clusters in its neighborhood.
#'
#' @param x A numeric matrix-like object containing observations in rows and variables in columns.
#' @param clusters Vector of length equal to \code{ncol(x)}, specifying the cluster assigned to each observation.
#' @param k Integer scalar specifying the number of nearest neighbors to use to determine the radius of the hyperspheres.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' This should be an algorithm supported by \code{\link{findNeighbors}}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating whether and how parallelization should be performed across genes.
#' @param weighted A logical scalar indicating whether to weight each observation in inverse proportion to the size of its cluster.
#' Alternatively, a numeric vector of length equal to \code{clusters} containing the weight to use for each observation.
#'
#' @return
#' A \linkS4class{DataFrame} with one row per observation in \code{x} and the columns:
#' \itemize{
#' \item \code{purity}, a numeric field containing the purity value for the current observation.
#' \item \code{maximum}, the cluster with the highest proportion of observations neighboring the current observation.
#' }
#' Row names are defined as the row names of \code{x}.
#' 
#' @details
#' The purity of a cluster is quantified by creating a hypersphere around each observation in the cluster
#' and computing the proportion of observations in that hypersphere from the same cluster.
#' If all observations in a cluster have proportions close to 1, this indicates that the cluster is highly pure,
#' i.e., there are few observations from other clusters in its region of the coordinate space.
#' The distribution of purities for each cluster can be used as a measure of separation from other clusters.
#'
#' In most cases, the majority of observations of a cluster will have high purities, corresponding to observations close to the cluster center.
#' A fraction of observations will have low values as these lie at the boundaries of two adjacent clusters.
#' A high degree of over-clustering will manifest as a majority of observations with purities close to zero.
#' The \code{maximum} field in the output can be used to determine the identity of the cluster 
#' with the greatest presence in a observation's neighborhood, usually an adjacent cluster for observations lying on the boundary.
#' 
#' The choice of \code{k} is used only to determine an appropriate value for the hypersphere radius.
#' We use hyperspheres as this is robust to changes in density throughout the coordinate space,
#' in contrast to computing purity based on the proportion of k-nearest neighbors in the same cluster.
#' For example, the latter will fail most obviously when the size of the cluster is less than \code{k}.
#'
#' @section Weighting by frequency:
#' By default, purity values are computed after weighting each observation by the reciprocal of the number of observations in the same cluster.
#' Otherwise, clusters with more observations will have higher purities as any contamination is offset by the bulk of observations,
#' which would compromise comparisons of purities between clusters.
#' One can interpret the weighted purities as the expected value after downsampling all clusters to the same size.
#'
#' Advanced users can achieve greater control by manually supplying a numeric vector of weights to \code{weighted}.
#' For example, we may wish to check the purity of batches after batch correction in single-cell RNA-seq.
#' In this application, \code{clusters} should be set to the \emph{batch blocking factor} (not the cluster identities!)
#' and \code{weighted} should be set to 1 over the frequency of each combination of cell type and batch.
#' This accounts for differences in cell type composition between batches when computing purities.
#' 
#' If \code{weighted=FALSE}, no weighting is performed.
#'
#' @author Aaron Lun
#' @examples
#' m <- matrix(runif(1000), ncol=10)
#' clusters <- clusterRows(m, BLUSPARAM=NNGraphParam())
#' out <- neighborPurity(m, clusters)
#' boxplot(split(out$purity, clusters))
#'
#' # Mocking up a stronger example:
#' centers <- matrix(rnorm(30), nrow=3)
#' clusters <- sample(1:3, 1000, replace=TRUE)
#' y <- centers[clusters,,drop=FALSE]
#' y <- y + rnorm(length(y))
#' 
#' out2 <- neighborPurity(y, clusters)
#' boxplot(split(out2$purity, clusters))
#'
#' @export
#' @importFrom BiocNeighbors KmknnParam buildIndex findKNN findNeighbors
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup
#' @importFrom stats median
#' @importFrom S4Vectors DataFrame
neighborPurity <- function(x, clusters, k=50, weighted=TRUE, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) {
    x <- as.matrix(x)
    idx <- buildIndex(x, BNPARAM=BNPARAM)

    if (!bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    dist <- median(findKNN(BNINDEX=idx, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, last=1, get.index=FALSE)$distance)
    nout <- findNeighbors(BNINDEX=idx, threshold=dist, BNPARAM=BNPARAM, get.distance=FALSE)$index

    # Constructing weights.
    if (isFALSE(weighted)) {
        w <- rep(1, nrow(x))
    } else if (isTRUE(weighted)) {
        w <- 1/table(clusters)
        w <- as.numeric(w[as.character(clusters)])
    } else {
        w <- as.numeric(weighted)
    }

    uclust <- sort(unique(clusters)) # do NOT use as.factor(), we want to preserve type.
    m <- match(clusters, uclust)
    aggregated <- sum_neighbor_weights(length(uclust), nout, m - 1L, w)
    targets <- t(aggregated[[1]])
    totals <- aggregated[[2]]

    DataFrame( 
        purity=targets[cbind(seq_along(m), m)]/totals,
        maximum=uclust[max.col(targets, ties.method="first")],
        row.names=rownames(x) 
    )
}
