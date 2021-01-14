#' Create a graph between different clusterings
#'
#' Create a graph that links together clusters from different clusterings,
#' e.g., generated using different parameter settings or algorithms.
#' This is useful for identifying corresponding clusters between clusterings
#' and to create meta-clusters from multiple clusterings.
#'
#' @param clusters A list of factors or vectors where each entry corresponds to a clustering.
#' All vectors should be of the same length.
#' The list itself should usually be named with a suitable label for each clustering.
#' @param prefix Logical scalar indicating whether the cluster levels should be prefixed with its clustering.
#' If \code{clusters} is not named, numeric prefixes are used instead.
#' @param denominator String specifying how the links between clusters should be computed.
#'
#' @details
#' Links are only formed between clusters from different clusterings.
#' The strength of the link between clusters \eqn{X} in clustering 1 and \eqn{Y} in clustering 2
#' is defined from the number of cells with those two labels in their respective clusterings.
#' A larger number of cells indicates that \eqn{X} and \eqn{Y} are corresponding clusters.
#'
#' Of course, the number of cells also depends on the total number of cells in each cluster.
#' To account for this, we normalize the strength by a function of the total number of cells in the two clusters.
#' The choice of function is determined by \code{denominator} and determines how the strength is adjusted for dissimilar cluster sizes.
#' \itemize{
#' \item For \code{"min"}, the number of shared cells is divided by the smaller of the totals between the two clusters.
#' \item For \code{"max"}, the number of shared cells is divided by the larger of the totals.
#' \item For \code{"union"}, the number of shared cells is divided by the size of the union.
#' The result is equivalent to the Jaccard index.
#' }
#'
#' In situations where \eqn{X} splits into multiple smaller clusters \eqn{Y1}, \eqn{Y2}, etc. in another clustering,
#' \code{denominator="min"} will report strong links between \eqn{X} and its constituent subclusters.
#' \code{denominator="max"} will report weak links and \code{denominator="union"} will be somewhere in between.
#' Conversely, \code{denominator="max"} forms the strongest links when there is a 1:1 mapping between clusters in different clusterings.
#' This usually yields simpler correspondences between clusterings at the cost of orphaning some of the smaller subclusters.
#'
#' The general idea is to use the graph returned by this function in visualization routines or for community-based clustering,
#' to identify \dQuote{clusters of clusters} that can inform about the relationships between clusterings.
#'
#' @return 
#' A \link{graph} object where each node is a cluster level in one of the clusterings in \code{clusters}.
#'
#' @author Aaron Lun
#'
#' @examples
#' clusters <- list(
#'     nngraph = clusterRows(iris[,1:4], NNGraphParam()),
#'     hclust = clusterRows(iris[,1:4], HclustParam(cut.dynamic=TRUE)),
#'     kmeans = clusterRows(iris[,1:4], KmeansParam(5))
#' )
#'
#' g <- linkClusters(clusters)
#' plot(g)
#'
#' igraph::cluster_walktrap(g)
#' @seealso
#' The \pkg{clustree} package, which provides another method for visualizing relationships between clusterings.
#'
#' \code{\link{compareClusterings}}, which computes similarities between the clusterings themselves.
#' 
#' @export
#' @importFrom igraph graph.adjacency
linkClusters <- function(clusters, prefix=TRUE, denominator=c("min", "union", "max")) {
    if (!length(unique(lengths(clusters)))) {
        stop("'clusters' must have elements of the same length")
    }

    totals <- lapply(clusters, table)
    N <- sum(lengths(totals))
    output <- matrix(0, N, N)

    newnames <- unlist(lapply(totals, names))
    if (prefix) {
        clust.names <- names(totals)
        if (is.null(clust.names)) {
            clust.names <- seq_along(totals)
        }
        newnames <- sprintf("%s.%s", rep(clust.names, lengths(totals)), newnames)
    }
    dimnames(output) <- list(newnames, newnames)

    denominator <- match.arg(denominator)
    lastx <- 0L

    for (x in seq_along(clusters)) {
        nx <- as.integer(totals[[x]])
        lasty <- 0L

        for (y in seq_len(x - 1L)) {
            tab <- table(clusters[[x]], clusters[[y]])
            ny <- as.integer(totals[[y]])

            if (denominator=="min") {
                denom <- outer(nx, ny, pmin)
            } else if (denominator=="union") {
                denom <- outer(nx, ny, `+`)
                denom <- denom - tab
            } else {
                denom <- outer(nx, ny, pmax)
            }

            ratio <- tab/denom
            output[lastx + seq_along(nx), lasty + seq_along(ny)] <- ratio
            lasty <- lasty + length(ny)
        }

        lastx <- lastx + length(nx)
    }

    graph.adjacency(output, mode="lower", weighted=TRUE)
}
