#' Hierarchical clustering
#'
#' Run the base \code{\link{hclust}} function on a distance matrix within \code{\link{clusterRows}}.
#'
#' @param metric String specifying the distance metric to use in \code{\link{dist}}.
#' @param method String specifying the agglomeration method to use in \code{\link{hclust}}.
#' @param cut.fun Function specifying the method to use to cut the dendrogram.
#' The first argument of this function should be the output of \code{\link{hclust}},
#' and the return value should be an atomic vector specifying the cluster assignment for each observation.
#' Defaults to \code{\link{cutree}}.
#' @param cut.height Numeric scalar specifying the cut height to use for the tree cut when \code{cut.fun=NULL}.
#' If \code{NULL}, defaults to half the tree height.
#' Ignored if \code{cut.number} is set.
#' @param cut.number Integer scalar specifying the number of clusters to generate from the tree cut when \code{cut.fun=NULL}.
#' @param ... Further arguments to pass to a non-\code{NULL} value for \code{cut.fun}.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{HclustParam} object.
#' @param full Logical scalar indicating whether the hierarchical clustering statistics should be returned.
#'
#' @author Aaron Lun
#'
#' @return 
#' The \code{HclustParam} constructor will return a \linkS4class{HclustParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects};
#' the latter is a list with the \code{hclust} entry containing the output of \code{\link{hclust}}.
#'
#' @examples
#' clusterRows(iris[,1:4], HclustParam())
#' clusterRows(iris[,1:4], HclustParam(method="ward.D2"))
#'
#' @seealso
#' \code{\link{dist}}, \code{\link{hclust}} and \code{\link{cutree}}, which actually do all the heavy lifting.
#'
#' \code{cutreeDynamic} from the \pkg{dynamicTreeCut} package, for an alternative tree cutting method to use in \code{cut.fun}.
#' @name HclustParam-class
NULL

#' @export
#' @rdname HclustParam-class
setClass("HclustParam", contains="BlusterParam", 
    slots=c(metric="character", method="character", cut.params="list"))

#' @export
#' @rdname HclustParam-class
HclustParam <- function(metric="euclidean", method="complete", cut.fun=NULL, cut.height=NULL, cut.number=NULL, ...) {
    new("HclustParam", metric=metric, method=method, 
        cut.params=list(fun=cut.fun, height=cut.height, number=cut.number, other=list(...)))
}

#' @export
#' @rdname HclustParam-class
#' @importFrom stats dist hclust cutree
setMethod("clusterRows", c("ANY", "HclustParam"), function(x, BLUSPARAM, full=FALSE) {
    dst <- dist(as.matrix(x), method=BLUSPARAM@metric)
    hcl <- hclust(dst, method=BLUSPARAM@method)

    cut.params <- BLUSPARAM@cut.params
    fun <- cut.params$fun

    if (is.null(fun)) {
        fun <- cutree
        k <- cut.params$number

        if (is.null(k)) {
            h <- cut.params$height
            if (is.null(h)) {
                h <- max(hcl$height)/2
            }
            args <- list(h=h)

        } else {
            args <- list(k=k)
        }
    } else {
        args <- cut.param$other
    }

    clusters <- do.call(fun, c(list(hcl), args))
    clusters <- factor(clusters)

    if (full) {
        list(clusters=clusters, objects=list(hclust=hcl))
    } else {
        clusters
    }
})
