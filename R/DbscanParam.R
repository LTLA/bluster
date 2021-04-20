#' Density-based clustering with DBSCAN
#'
#' Perform density-based clustering with a fast re-implementation of the DBSCAN algorithm.
#'
#' @inheritParams clusterRows
#' @param eps Numeric scalar specifying the distance to use to define neighborhoods.
#' If \code{NULL}, this is determined from \code{min.pts} and \code{core.prop}.
#' @param min.pts Integer scalar specifying the minimum number of neighboring observations required for an observation to be a core point.
#' @param core.prop Numeric scalar specifying the proportion of observations to treat as core points.
#' This is only used when \code{eps=NULL}, see Details.
#' @param chunk.size Integer scalar specifying the number of points to process per chunk.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for the neighbor searches.
#' This should be able to support both nearest-neighbor and range queries.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the parallelization to use for the neighbor searches.
#' @param full Logical scalar indicating whether additional statistics should be returned.
#'
#' @details
#' DBSCAN operates by identifying core points, i.e., observations with at least \code{min.pts} neighbors within a distance of \code{eps}.
#' It then identifies which core points are neighbors of each other, forming components of connected core points.
#' All non-core points are then connected to the closest core point within \code{eps}.
#' All groups of points that are connected in this manner are considered to be part of the same cluster.
#' Any unconnected non-core points are treated as noise and reported as \code{NA}.
#'
#' As a suitable value of \code{eps} may not be known beforehand, we can automatically determine it from the data.
#' For all observations, we compute the distance to the \eqn{k}th neighbor where \eqn{k} is defined as \code{round(min.pts * core.prop)}.
#' We then define \code{eps} as the \code{core.prop} quantile of the distances across all observations.
#' The default of \code{core.prop=0.5} means that around half of the observations will be treated as core points.
#'
#' Larger values of \code{eps} will generally result in fewer observations classified as noise, as they are more likely to connect to a core point.
#' It may also promote agglomeration of existing clusters into larger entities if they are connected by regions of (relatively) low density.
#' Conversely, larger values of \code{min.pts} will generally increase the number of noise points and may fragment larger clusters into subclusters.
#' 
#' To modify an existing DbscanParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#' 
#' @return
#' The \code{DbscanParam} constructor will return a \linkS4class{DbscanParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' Note that this may contain \code{NA} values corresponding to noise points.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects} 
#' (a list containing the \code{eps} and \code{min.pts} used in the analysis).
#' 
#' @examples
#' clusterRows(iris[,1:4], DbscanParam())
#' clusterRows(iris[,1:4], DbscanParam(core.prop=0.8))
#'
#' @references
#' Ester M et al. (1996).
#' A density-based algorithm for discovering clusters in large spatial databases with noise.
#' \emph{Proceedings of the Second International Conference on Knowledge Discovery and Data Mining}, 226-231.
#' 
#' @author Aaron Lun
#' @name DbscanParam-class
#' @docType class
#' @aliases
#' show,DbscanParam-method
NULL

#' @export
setClass("DbscanParam", contains="BlusterParam", slots=c(eps="numeric_OR_NULL", min.pts="integer", core.prop="numeric", 
    chunk.size="integer", BNPARAM="BiocNeighborParam", BPPARAM="BiocParallelParam"))

#' @export
#' @rdname DbscanParam-class
DbscanParam <- function(eps=NULL, min.pts=5, core.prop=0.5, chunk.size=1000, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) {
    new("DbscanParam", eps=eps, min.pts=as.integer(min.pts), 
        core.prop=core.prop, chunk.size=as.integer(chunk.size),
        BNPARAM=BNPARAM, BPPARAM=BPPARAM)
}

#' @export
setMethod("show", "DbscanParam", function(object) {
    callNextMethod()
    cat(sprintf("eps: %s\n", if (is.null(object[["eps"]])) "default" else object[["eps"]]))
    for (i in c("min.pts", "core.prop", "chunk.size")) {
        cat(sprintf("%s: %s\n", i, object[[i]]))
    }
    for (i in c("BNPARAM", "BPPARAM")) {
        cat(sprintf("%s: %s\n", i, class(object[[i]])[1]))
    }
})

setValidity2("DbscanParam", function(object) {
    msg <- character(0)

    if (!is.null(object[["eps"]])) {
        msg <- c(msg, .check_positive_slots(object, "eps"))
    }

    msg <- c(msg, .check_positive_slots(object, c("min.pts", "core.prop", "chunk.size")))

    if (object[["core.prop"]] > 1) {
        msg <- c(msg, "'core.prop' should not be greater than 1")
    }

    if (length(msg)) {
        return(msg)
    }
    TRUE
})

#' @export
#' @rdname DbscanParam-class
setMethod("clusterRows", c("ANY", "DbscanParam"), function(x, BLUSPARAM, full=FALSE) {
    out <- .DBSCAN(x, 
        eps=BLUSPARAM[["eps"]],
        min.pts=BLUSPARAM[["min.pts"]],
        core.prop=BLUSPARAM[["core.prop"]],
        chunk.size=BLUSPARAM[["chunk.size"]],
        BNPARAM=BLUSPARAM[["BNPARAM"]],
        BPPARAM=BLUSPARAM[["BPPARAM"]]
    )

    clusters <- out$clusters
    clusters[clusters==0] <- NA_integer_
    clusters <- factor(clusters)

    if (full) {
        out$clusters <- NULL
        list(clusters=clusters, objects=out)
    } else {
        clusters
    }
})

#' @importFrom BiocNeighbors buildIndex findNeighbors findKNN KmknnParam queryKNN
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup
#' @importFrom stats quantile
#' @importFrom igraph make_graph components
#' @importFrom utils head
.DBSCAN <- function(x, eps=NULL, min.pts=5, core.prop=0.5, chunk.size=1000, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) {
    if (!bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam")) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Finding all core points, using a quick NN algorithm.
    all.dist <- findKNN(X=x, k=min.pts, last=1, get.index=FALSE, BNPARAM=BNPARAM, BPPARAM=BPPARAM)$distance
    if (is.null(eps)) {
        eps <- quantile(all.dist, core.prop)
    } 
    is.core <- all.dist <= eps

    # Looping through and finding core-core neighbors in chunks to keep memory usage low.
    n.core <- sum(is.core)
    core.idx <- buildIndex(x[is.core,,drop=FALSE], BNPARAM=BNPARAM)
    remaining <- seq_len(n.core)
    to <- vector("list", n.core)

    while (length(remaining)) {
        chosen <- head(remaining, chunk.size)
        neighbors <- findNeighbors(BNINDEX=core.idx, subset=chosen, threshold=eps, get.distance=FALSE, BPPARAM=BPPARAM)$index
        to[chosen] <- neighbors
        remaining <- setdiff(remaining, unlist(neighbors)) 
    }

    # Identifying core-core components.
    from <- rep(seq_along(to), lengths(to))
    to <- unlist(to)
    g <- make_graph(rbind(from, to), n=n.core, directed=FALSE)
    core.clusters <- components(g)$membership

    # Finding closest core point (within range) for all non-core points.
    clusters <- integer(nrow(x))
    clusters[is.core] <- core.clusters

    if (n.core) {
        closest.core <- queryKNN(BNINDEX=core.idx, query=x[!is.core,,drop=FALSE], BPPARAM=BPPARAM, k=1)
        in.range <- closest.core$distance <= eps
        clusters[which(!is.core)[in.range]] <- core.clusters[closest.core$index[in.range]]
    }

    list(clusters=clusters, eps=eps, min.pts=min.pts)
}
