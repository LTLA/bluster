#' Hierarchical clustering
#'
#' Run the base \code{\link{hclust}} function on a distance matrix within \code{\link{clusterRows}}.
#'
#' @param metric String specifying the distance metric to use in \code{\link{dist}}.
#' @param method String specifying the agglomeration method to use in \code{\link{hclust}}.
#' @param cut.fun Function specifying the method to use to cut the dendrogram.
#' The first argument of this function should be the output of \code{\link{hclust}},
#' and the return value should be an atomic vector specifying the cluster assignment for each observation.
#' Defaults to \code{\link{cutree}} if \code{cut.dynamic=FALSE} and \code{\link[dynamicTreeCut]{cutreeDynamic}} otherwise.
#' @param cut.dynamic Logical scalar indicating whether a dynamic tree cut should be performed using the \pkg{dynamicTreeCut} package.
#' @param cut.height Numeric scalar specifying the cut height to use for the tree cut when \code{cut.fun=NULL}.
#' If \code{NULL}, defaults to half the tree height.
#' Ignored if \code{cut.number} is set.
#' @param cut.number Integer scalar specifying the number of clusters to generate from the tree cut when \code{cut.fun=NULL}.
#' @param cut.params Further arguments to pass to \code{cut.fun}, when \code{cut.dynamic=TRUE} or \code{cut.fun} is non-\code{NULL}.
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
#'
#' @seealso
#' \code{\link{dist}}, \code{\link{hclust}} and \code{\link{cutree}}, which actually do all the heavy lifting.
#'
#' \code{\link[dynamicTreeCut]{cutreeDynamic}}, for an alternative tree cutting method to use in \code{cut.fun}.
#' @name HclustParam-class
#' @docType class
#' @aliases 
#' show,HclustParam-method
NULL

#' @export
#' @rdname HclustParam-class
setClass("HclustParam", contains="BlusterParam", 
    slots=c(metric="character", method="character", cut.fun="function_OR_NULL", cut.dynamic="logical",
        cut.height="numeric_OR_NULL", cut.number="integer_OR_NULL", cut.params="list"))

#' @export
#' @rdname HclustParam-class
HclustParam <- function(metric="euclidean", method="complete", 
    cut.fun=NULL, cut.dynamic=FALSE, cut.height=NULL, cut.number=NULL, 
    cut.params=list(), ...)
{
    if (!is.null(cut.number)) {
        cut.number <- as.integer(cut.number)
    }

    extra.args <- list(...)
    if (length(extra.args)) {
        .Deprecated(old="...", new="cut.params=")
        cut.params <- c(cut.params, extra.args)
    }

    new("HclustParam", metric=metric, method=method, cut.fun=cut.fun, cut.dynamic=cut.dynamic,
        cut.height=cut.height, cut.number=cut.number, cut.params=cut.params)
}

#' @importFrom S4Vectors setValidity2
setValidity2("HclustParam", function(object) {
    msg <- character(0)

    for (i in c("metric", "method")) {
        if (!.non_na_scalar(slot(object, i))) {
            msg <- c(msg, sprintf("'%s' must be a non-missing string", i))
        }
    }

    if (!.non_na_scalar(slot(object, "cut.dynamic"))) {
        msg <- c(msg, sprintf("'%s' must be a non-missing logical scalar", i))
    }

    h <- object@cut.height
    if (!is.null(h) && !.positive_number(h)) {
        msg <- c(msg, "'cut.height' must be NULL or a positive number")
    }

    k <- object@cut.number
    if (!is.null(k) && !.positive_number(k)) {
        msg <- c(msg, "'cut.number' must be NULL or a positive number")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
#' @importFrom S4Vectors coolcat
setMethod("show", "HclustParam", function(object) {
    callNextMethod()
    cat(sprintf("metric: %s\n", object@metric))
    cat(sprintf("method: %s\n", object@method))

    fun <- object@cut.fun
    if (is.null(fun)) {
        if (object@cut.dynamic) {
            cat("cut.fun: cutreeDynamic\n")
            coolcat("cut.params(%i): %s", names(object@cut.params))
        } else {
            cat("cut.fun: cutree\n")
            k <- object@cut.number

            if (is.null(k)) {
                h <- object@cut.height
                if (is.null(h)) h <- "default"
                cat(sprintf("cut.height: %s\n", h))
            } else {
                cat(sprintf("cut.number: %s\n", k))
            }
        }

    } else {
        cat("cut.fun: custom\n")
        coolcat("cut.params(%i): %s", names(object@cut.params))
    }
})

#' @export
#' @rdname HclustParam-class
#' @importFrom stats dist hclust cutree
setMethod("clusterRows", c("ANY", "HclustParam"), function(x, BLUSPARAM, full=FALSE) {
    dst <- dist(as.matrix(x), method=BLUSPARAM@metric)
    hcl <- hclust(dst, method=BLUSPARAM@method)

    fun <- BLUSPARAM@cut.fun
    if (is.null(fun)) {
        if (!BLUSPARAM@cut.dynamic) {
            fun <- cutree
            k <- BLUSPARAM@cut.number

            if (is.null(k)) {
                h <- BLUSPARAM@cut.height
                if (is.null(h)) {
                    h <- max(hcl$height)/2
                }
                args <- list(h=h)

            } else {
                args <- list(k=k)
            }
        } else {
            fun <- function(...) unname(dynamicTreeCut::cutreeDynamic(..., verbose=0))
            args <- c(list(dist=as.matrix(dst)), BLUSPARAM@cut.params)
        }
    } else {
        args <- BLUSPARAM@cut.params
    }

    clusters <- do.call(fun, c(list(hcl), args))
    clusters <- factor(clusters)

    if (full) {
        list(clusters=clusters, objects=list(dist=dst, hclust=hcl))
    } else {
        clusters
    }
})
