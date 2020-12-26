#' Graph-based clustering
#'
#' Run community detection algorithms on a nearest-neighbor (NN) graph within \code{\link{clusterRows}}.
#'
#' @param shared Logical scalar indicating whether a shared NN graph should be constructed.
#' @param ... Further arguments to pass to \code{\link{SNNGraphParam}} (if \code{shared=TRUE}) or \code{\link{KNNGraphParam}}.
#' @inheritParams makeSNNGraph
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
#' The SNNGraphParam and KNNGraphParam classes are both derived from the NNGraphParam virtual class.
#' This former will perform clustering with a shared nearest-neighbor (SNN) graph 
#' while the latter will use a simpler k-nearest neighbor (KNN) graph - see \code{?\link{makeSNNGraph}} for details.
#'
#' To modify an existing NNGraphParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#' 
#' @return 
#' The constructors will return a \linkS4class{NNGraphParam} object with the specified parameters.
#' If \code{shared=TRUE}, this is a SNNGraphParam object; otherwise it is a KNNGraphParam object.
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
setClass("NNGraphParam", contains=c("BlusterParam", "VIRTUAL"),
    slots=c(k="integer", cluster.fun="character_OR_function", cluster.args="list",
        BNPARAM="BiocNeighborParam", BPPARAM="BiocParallelParam"))

#' @export
#' @rdname NNGraphParam-class
setClass("SNNGraphParam", contains="NNGraphParam", slot=c(type="character"))

#' @export
#' @rdname NNGraphParam-class
setClass("KNNGraphParam", contains="NNGraphParam", slot=c(directed="logical"))

#' @export
#' @rdname NNGraphParam-class
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam
NNGraphParam <- function(shared=TRUE, k=10, ..., BNPARAM=KmknnParam(), BPPARAM=SerialParam(), cluster.fun="walktrap", cluster.args=list()) {
    if (shared) {
        SNNGraphParam(k=k, ..., BNPARAM=BNPARAM, BPPARAM=BPPARAM, cluster.fun=cluster.fun, cluster.args=cluster.args)
    } else {
        KNNGraphParam(k=k, ..., BNPARAM=BNPARAM, BPPARAM=BPPARAM, cluster.fun=cluster.fun, cluster.args=cluster.args)
    }
}

#' @export
#' @rdname NNGraphParam-class
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam
SNNGraphParam <- function(k=10, type="rank", BNPARAM=KmknnParam(), BPPARAM=SerialParam(), cluster.fun="walktrap", cluster.args=list()) {
    new("SNNGraphParam", k=as.integer(k), type=type, BNPARAM=BNPARAM, BPPARAM=BPPARAM, cluster.fun=cluster.fun, cluster.args=cluster.args)
}

#' @export
#' @rdname NNGraphParam-class
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam
KNNGraphParam <- function(k=10, directed=FALSE, BNPARAM=KmknnParam(), BPPARAM=SerialParam(), cluster.fun="walktrap", cluster.args=list()) {
    new("KNNGraphParam", k=as.integer(k), directed=directed, BNPARAM=BNPARAM, BPPARAM=BPPARAM, cluster.fun=cluster.fun, cluster.args=cluster.args)
}

#' @importFrom S4Vectors setValidity2
setValidity2("NNGraphParam", function(object) {
    msg <- character(0)

    if (!.positive_number(object@k)) {
        msg <- c(msg, "'k' must be a positive integer scalar")
    }

    cluster.fun <- object@cluster.fun
    if (!is.function(cluster.fun) && !.non_na_scalar(cluster.fun)) {
        msg <- c(msg, "'cluster.fun' must be a non-missing string")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @importFrom S4Vectors setValidity2
setValidity2("SNNGraphParam", function(object) {
    msg <- character(0)

    if (!.non_na_scalar(object@type)) {
        msg <- c(msg, "'type' should be a non-NA string")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @importFrom S4Vectors setValidity2
setValidity2("KNNGraphParam", function(object) {
    msg <- character(0)

    if (!.non_na_scalar(object@directed)) {
        msg <- c(msg, "'directed' should be a non-NA logical scalar")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
#' @importFrom S4Vectors coolcat
setMethod("show", "NNGraphParam", function(object) {
    callNextMethod()
    cat(sprintf("k: %s\n", object@k))
    sub_graph_show(object)
    cat(sprintf("BNPARAM: %s\n", class(object@BNPARAM)[1]))
    cat(sprintf("BPPARAM: %s\n", class(object@BPPARAM)[1]))
    cat(sprintf("cluster.fun: %s\n", if (is.function(object@cluster.fun)) "custom" else object@cluster.fun))
    coolcat("cluster.args(%i): %s", names(object@cluster.args))
})

# Can't quite be bothered to capture the output and re-print it, so I'm just
# going to define more generics to handle the right ordering of show outputs.

setGeneric("sub_graph_show", function(object) standardGeneric("sub_graph_show"))

setMethod("sub_graph_show", "SNNGraphParam", function(object) {
    cat(sprintf("type: %s\n", object@type))
})

setMethod("sub_graph_show", "KNNGraphParam", function(object) {
    cat(sprintf("directed: %s\n", object@directed))
})

#' @export
#' @rdname NNGraphParam-class
setMethod("clusterRows", c("ANY", "SNNGraphParam"), function(x, BLUSPARAM, full=FALSE) {
    g <- makeSNNGraph(x, k=BLUSPARAM[["k"]], type=BLUSPARAM[["type"]],
        BNPARAM=BLUSPARAM[["BNPARAM"]], BPPARAM=BLUSPARAM[["BPPARAM"]]) 
    .cluster_igraph(g, BLUSPARAM, full=full)
})

#' @export
#' @rdname NNGraphParam-class
setMethod("clusterRows", c("ANY", "KNNGraphParam"), function(x, BLUSPARAM, full=FALSE) {
    g <- makeKNNGraph(x, k=BLUSPARAM[["k"]], directed=BLUSPARAM[["directed"]],
        BNPARAM=BLUSPARAM[["BNPARAM"]], BPPARAM=BLUSPARAM[["BPPARAM"]]) 
    .cluster_igraph(g, BLUSPARAM, full=full)
})

#' @importFrom utils getFromNamespace
#' @importFrom igraph membership cluster_walktrap
.cluster_igraph <- function(g, BLUSPARAM, full=TRUE) {
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
}
