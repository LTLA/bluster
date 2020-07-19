#' Graph-based clustering
#'
#' Run community detection algorithms on a nearest-neighbor (NN) graph within \code{\link{clusterRows}}.
#'
#' @param shared Logical scalar indicating whether a shared NN graph should be constructed.
#' @param ... Further arguments to pass to \code{\link{makeSNNGraph}} (if \code{shared=TRUE}) or \code{\link{makeKNNGraph}}.
#' @param cluster.fun Function specifying the method to use to detect communities in the NN graph.
#' The first argument of this function should be the NN graph and the return value should be a \link{communities} object.
#'
#' Alternatively, this may be a string containing the suffix of any \pkg{igraph} community detection algorithm.
#' For example, \code{cluster.fun="louvain"} will instruct \code{\link{clusterRows}} to use \code{\link{cluster_louvain}}.
#' Defaults to \code{\link{cluster_walktrap}}.
#' @param cluster.args Further arguments to pass to the chosen \code{cluster.fun}.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{NNGraphParam} object.
#' @param full Logical scalar indicating whether the graph-based clustering objects should be returned.
#'
#' @author Aaron Lun
#'
#' @details
#' To modify an existing NNGraphParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
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
#' \code{\link{makeSNNGraph}} and related functions, to build the graph.
#'
#' \code{\link{cluster_walktrap}} and related functions, to perform community detection.
#'
#' @name NNGraphParam-class
#' @docType class
#' @aliases
#' show,NNGraphParam-method
NULL

#' @export
#' @rdname NNGraphParam-class
setClass("NNGraphParam", contains="BlusterParam", 
    slots=c(shared="logical", graph.args="list", cluster.fun="character_OR_function", cluster.args="list"))

#' @export
#' @rdname NNGraphParam-class
NNGraphParam <- function(shared=TRUE, ..., cluster.fun="walktrap", cluster.args=list()) {
    new("NNGraphParam", shared=shared, graph.args=list(...), 
        cluster.fun=cluster.fun, cluster.args=cluster.args)
}

setMethod(".extras", "NNGraphParam", function(x) "graph.args")

#' @importFrom S4Vectors setValidity2
setValidity2("NNGraphParam", function(object) {
    msg <- character(0)

    if (!.non_na_scalar(object@shared)) {
        msg <- c(msg, "'shared' must be a non-missing logical scalar")
    }

    cluster.fun <- object@cluster.fun
    if (!is.function(cluster.fun) && !.non_na_scalar(cluster.fun)) {
        msg <- c(msg, "'cluster.fun' must be a non-missing string")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
#' @importFrom S4Vectors coolcat
setMethod("show", "NNGraphParam", function(object) {
    callNextMethod()
    cat(sprintf("shared: %s\n", object@shared))
    coolcat("graph.args(%i): %s", names(object@graph.args))
    cat(sprintf("cluster.fun: %s\n", if (is.function(object@cluster.fun)) "custom" else object@cluster.fun))
    coolcat("cluster.args(%i): %s", names(object@cluster.args))
})

#' @export
#' @rdname NNGraphParam-class
#' @importFrom utils getFromNamespace
#' @importFrom igraph membership cluster_walktrap
setMethod("clusterRows", c("ANY", "NNGraphParam"), function(x, BLUSPARAM, full=FALSE) {
    if (BLUSPARAM@shared) {
        FUN <- makeSNNGraph
    } else {
        FUN <- makeKNNGraph
    }
    g <- do.call(FUN, c(list(x), BLUSPARAM@graph.args))

    cFUN <- BLUSPARAM@cluster.fun
    if (is.character(cFUN)) {
        cFUN <- getFromNamespace(paste0("cluster_", cFUN), ns="igraph")
    }

    comm <- do.call(cFUN, c(list(g), BLUSPARAM@cluster.args))
    clusters <- factor(membership(comm))

    if (full) {
        list(clusters=clusters, objects=list(graph=g, communities=comm))
    } else {
        clusters
    }
})
