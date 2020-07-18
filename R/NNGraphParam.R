#' Graph-based clustering
#'
#' Run community detection algorithms on a nearest-neighbor (NN) graph within \code{\link{clusterRows}}.
#'
#' @param shared Logical scalar indicating whether a shared NN graph should be constructed.
#' @param ... Further arguments to pass to \code{\link{makeSNNGraph}} (if \code{shared=TRUE}) or \code{\link{makeKNNGraph}}.
#' @param cluster.fun Function specifying the method to use to detect communities in the NN graph.
#' The first argument of this function should be the NN graph and the return value should be a \link{communities} object.
#' Defaults to \code{\link{cluster_walktrap}}.
#'
#' Alternatively, this may be a string containing the suffix of any \pkg{igraph} community detection algorithm.
#' For example, \code{cluster.fun="louvain"} will instruct \code{\link{clusterRows}} to use \code{\link{cluster_louvain}}.
#' @param cluster.args Further arguments to pass to the chosen \code{cluster.fun}.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{NNGraphParam} object.
#' @param full Logical scalar indicating whether the graph-based clustering objects should be returned.
#'
#' @author Aaron Lun
#'
#' @return 
#' The \code{NNGraphParam} constructor will return a \linkS4class{NNGraphParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects};
#' the latter is a list with \code{graph} (the graph) and \code{communities} (the output of \code{cluster.fun}).
#'
#' @examples
#' clusterRows(iris[,1:4], NNGraphParam())
#' clusterRows(iris[,1:4], NNGraphParam(k=5))
#' clusterRows(iris[,1:4], NNGraphParam(cluster.fun="louvain"))
#'
#' @seealso
#' \code{\link{dist}}, \code{\link{hclust}} and \code{\link{cutree}}, which actually do all the heavy lifting.
#'
#' \code{cutreeDynamic} from the \pkg{dynamicTreeCut} package, for an alternative tree cutting method to use in \code{cut.fun}.
#' @name NNGraphParam-class
NULL

#' @export
#' @rdname NNGraphParam-class
setClass("NNGraphParam", contains="BlusterParam", 
    slots=c(shared="logical", graph.args="list", cluster.params="list"))

#' @export
#' @rdname NNGraphParam-class
NNGraphParam <- function(shared=TRUE, ..., cluster.fun=NULL, cluster.args=list()) {
    new("NNGraphParam", shared=shared, graph.args=list(...), cluster.params=list(fun=cluster.fun, args=cluster.args))
}

#' @export
#' @rdname NNGraphParam-class
#' @importFrom igraph membership cluster_walktrap
setMethod("clusterRows", c("ANY", "NNGraphParam"), function(x, BLUSPARAM, full=FALSE) {
    if (BLUSPARAM@shared) {
        FUN <- makeSNNGraph
    } else {
        FUN <- makeKNNGraph
    }
    g <- do.call(FUN, c(list(x), BLUSPARAM@graph.args))

    cluster.params <- BLUSPARAM@cluster.params
    cFUN <- cluster.params$fun
    if (is.null(cFUN)) {
        cFUN <- cluster_walktrap
    } else if (is.character(cFUN)) {
        cFUN <- getFromNamespace(paste0("cluster_", cFUN), ns="igraph")
    } else if (!is.function(cFUN)) {
        stop("'cluster.fun' should be 'NULL', a string or a function")
    }
    comm <- do.call(cFUN, c(list(g), cluster.params$args))

    clusters <- factor(membership(comm))

    if (full) {
        list(clusters=clusters, objects=list(graph=g, communities=comm))
    } else {
        clusters
    }
})
