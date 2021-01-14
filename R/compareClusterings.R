#' Compare pairs of clusterings
#'
#' Compute the adjusted Rand index between all pairs of clusterings,
#' where larger values indicate a greater similarity between clusterings.
#'
#' @inheritParams linkClusters
#' @inheritParams pairwiseRand
#' 
#' @details
#' The aim of this function is to allow us to easily determine the relationships between clusterings.
#' For example, we might use this to determine which parameter settings have the greatest effect in a sweep by \code{\link{clusterSweep}}.
#' Alternatively, we could use this to obtain an \dQuote{ordering} of clusterings for visualization, e.g., with \pkg{clustree}.
#'
#' This function does not provide any insight into the relationships between individual clusters.
#' A large Rand index only means that two clusterings are similar but does not specify the corresponding set of clusters across clusterings.
#' For that task, we suggest using the \code{\link{linkClusters}} function instead.
#'
#' @return A symmetric square matrix of pairwise (adjusted) Rand indices between all pairs of clusters.
#'
#' @return Aaron Lun
#'
#' @examples
#' clusters <- list(
#'     nngraph = clusterRows(iris[,1:4], NNGraphParam()),
#'     hclust = clusterRows(iris[,1:4], HclustParam(cut.dynamic=TRUE)),
#'     kmeans = clusterRows(iris[,1:4], KmeansParam(20))
#' )
#'
#' aris <- compareClusterings(clusters)
#'
#' # Visualizing the relationships between clusterings.
#' # Here, k-means is forced to be least similar to the others.
#' ari.as.graph <- igraph::graph.adjacency(aris, mode="undirected", weighted=TRUE)
#' plot(ari.as.graph)
#'
#' # Obtain an ordering of clusterings, using the eigenvector 
#' # as a 1-dimensional summary of the matrix:
#' ev1 <- eigen(aris)$vectors[,1]
#' o <- order(ev1)
#' rownames(aris)[o]
#'
#' @seealso
#' \code{\link{linkClusters}}, which identifies relationships between individual clusters across clusterings.
#'
#' \code{\link{pairwiseRand}}, for calculation of the pairwise Rand index.
#' 
#' @export
compareClusterings <- function(clusters, adjusted=TRUE) {
    if (!length(unique(lengths(clusters)))) {
        stop("'clusters' must have elements of the same length")
    }

    xnames <- names(clusters)
    N <- length(xnames)
    output <- matrix(0, N, N, dimnames=list(xnames, xnames))

    for (x in seq_len(N)) {
        for (y in seq_len(x-1)) {
            output[x,y] <- pairwiseRand(clusters[[x]], clusters[[y]], mode="index", adjusted=adjusted)
        }
    }

    # Making it symmetric.
    output + t(output)
}
