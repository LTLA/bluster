#' Hierarchical clustering
#'
#' Run the base \code{\link{hclust}} function on a distance matrix within \code{\link{clusterRows}}.
#'
#' @param clust.fun Function specifying the function to use to do the clustering.
#' The function should apply a hierarchical clustering algorithm and take a data matrix as input.
#' If \code{NULL}, the \code{stats::\link{hclust}} function is used by default.
#' @param metric String specifying the distance metric to use in \code{dist.fun}.
#' If \code{NULL}, the default method of \code{dist.fun} is used.
#' @param dist.fun Function specifying the function to use to compute the distance matrix. 
#' The function should accept a data matrix and a \code{method=} string (used to accept \code{metric}) and return a dissimilarity matrix of type \link{dist}.
#' If \code{NULL}, the \code{stats::\link{dist}} function is used by default.
#' @param method String specifying the agglomeration method to use in \code{\link{hclust}}.
#' @param cut.fun Function specifying the method to use to cut the dendrogram.
#' The first argument of this function should be the output of \code{\link{hclust}},
#' and the return value should be an atomic vector specifying the cluster assignment for each observation.
#' Defaults to \code{\link{cutree}} if \code{cut.dynamic=FALSE} and \code{\link[dynamicTreeCut]{cutreeDynamic}} otherwise.
#' @param cut.dynamic Logical scalar indicating whether a dynamic tree cut should be performed using the \pkg{dynamicTreeCut} package.
#' @param cut.height,cut.number Deprecated, use \code{h} and \code{k} in \code{cut.params} instead.
#' @param cut.params Further arguments to pass to \code{cut.fun}.
#' @param ... Deprecated, more arguments to add to \code{cut.params}.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{HclustParam} object.
#' @param full Logical scalar indicating whether the hierarchical clustering statistics should be returned.
#'
#' @author Aaron Lun
#'
#' @details
#' To modify an existing HclustParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' If \code{cut.fun=NULL}, \code{cut.dynamic=FALSE} and \code{cut.params} does not have \code{h} or \code{k},
#' \code{\link{clusterRows}} will automatically set \code{h} to half the tree height when calling \code{\link{cutree}}. 
#'
#' @return
#' The \code{HclustParam} constructor will return a \linkS4class{HclustParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects} 
#' (a list containing \code{dist}, the distance matrix; and \code{hclust}, the output of \code{\link{hclust}}).
#'
#' @examples
#' clusterRows(iris[,1:4], HclustParam())
#' clusterRows(iris[,1:4], HclustParam(method="ward.D2"))
#' clusterRows(iris[,1:4], HclustParam(metric = "canberra", dist.fun = vegan::vegdist))
#' clusterRows(iris[,1:4], HclustParam(clust.fun=fastcluster::hclust))
#' 
#' @seealso
#' \code{\link{dist}}, \code{\link{hclust}} and \code{\link{cutree}}, which actually do all the heavy lifting.
#'
#' \code{\link[dynamicTreeCut]{cutreeDynamic}}, for an alternative tree cutting method to use in \code{cut.fun}.
#' @name HclustParam-class
#' @docType class
#' @aliases 
#' .defaultScalarArguments,HclustParam-method
#' updateObject,HclustParam-method
#' show,HclustParam-method
#' [[,HclustParam-method
NULL

#' @export
#' @rdname HclustParam-class
setClass("HclustParam", contains="HierarchicalParam", slots=c(method="ANY"))

#' @export
setMethod(".defaultScalarArguments", "HclustParam", function(x) c(callNextMethod(), method="character")) 

#' @export
#' @rdname HclustParam-class
HclustParam <- function(clust.fun=NULL, metric=NULL, dist.fun=NULL, method=NULL, cut.fun=NULL, cut.dynamic=FALSE, cut.height=NULL, cut.number=NULL, cut.params=list(), ...) {
    if (!is.null(cut.number)) {
        .Deprecated(old="cut.number=", new="cut.params=list(k=...)")
        cut.params$k <- cut.number
    } else if (!is.null(cut.height)) {
        .Deprecated(old="cut.height=", new="cut.params=list(h=...)")
        cut.params$h <- cut.height
    }

    extra.args <- list(...)
    if (length(extra.args)) {
        .Deprecated(old="...", new="cut.params=")
        cut.params <- c(cut.params, extra.args)
    }

    new("HclustParam", clust.fun=clust.fun, metric=metric, dist.fun=dist.fun, method=method, cut.fun=cut.fun, cut.dynamic=cut.dynamic, cut.params=cut.params)
}

#' @export
setMethod("[[", "HclustParam", function(x, i) {
    x <- updateObject(x)
    callNextMethod()
})

#' @export
setMethod("show", "HclustParam", function(object) {
    object <- updateObject(object)
    callNextMethod()
    dist.fun <- object@dist.fun
    if (!is.null(dist.fun)) {
        cat("dist.fun: custom\n")
    } else {
        cat("dist.fun: stats::dist\n")
    }

    clust.fun <- object@clust.fun
    if (!is.null(clust.fun)) {
        cat("clust.fun: custom\n")
    } else {
        cat("clust.fun: stats::hclust\n")
    }
})

#' @export
#' @rdname HclustParam-class
#' @importFrom stats dist hclust cutree
setMethod("clusterRows", c("ANY", "HclustParam"), function(x, BLUSPARAM, full=FALSE) {
    dargs <- list(quote(as.matrix(x)))
    if (!is.null(BLUSPARAM@metric)) {
        dargs$method <- BLUSPARAM@metric
    }
    
    if (!is.null(BLUSPARAM@dist.fun)) {
        dst <- do.call(BLUSPARAM@dist.fun, dargs)
    } else {
        dst <- do.call(dist, dargs)
    }

    if (!is.null(BLUSPARAM@clust.fun)) {
        clust.fun <- BLUSPARAM@clust.fun
    } else {
        clust.fun <- hclust
    }

    hargs <- list(quote(dst))
    if (!is.null(BLUSPARAM@method)) {
        hargs$method <- BLUSPARAM@method
    }
    hcl <- do.call(clust.fun, hargs)

    # Moving arguments to their new homes.
    BLUSPARAM <- updateObject(BLUSPARAM)

    clusters <- .cut_hierarchical(hcl, dst, BLUSPARAM)
    if (full) {
        list(clusters=clusters, objects=list(dist=dst, hclust=hcl))
    } else {
        clusters
    }
})

#' @export
#' @importFrom S4Vectors updateObject
setMethod("updateObject", "HclustParam", function(object, ..., verbose=FALSE) {
    has.k <- .hasSlot(object, "cut.number")
    has.h <- .hasSlot(object, "cut.height")
    
    if (has.k || has.h) {
        if (has.k && !is.null(object@cut.number)) {
            object@cut.params$k <- object@cut.number
            if (verbose) {
                message("[updateObject] moving 'cut.number' to 'cut.params$k'")
            }
        } else if (!is.null(object@cut.height)) {
            object@cut.params$h <- object@cut.height
            if (verbose) {
                message("[updateObject] moving 'cut.height' to 'cut.params$h'")
            }
        }
    }

    object
})
