#' Divisive analysis clustering
#'
#' Use the \code{\link{diana}} function to perform divisive analysis clustering.
#'
#' @inheritParams HclustParam
#' @param stand Further arguments to pass to \code{\link{diana}}.
#'
#' @author Aaron Lun
#'
#' @details
#' To modify an existing DianaParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' If \code{cut.fun=NULL}, \code{cut.dynamic=FALSE} and \code{cut.params} does not have \code{h} or \code{k},
#' \code{\link{clusterRows}} will automatically set \code{h} to half the tree height when calling \code{\link{cutree}}. 
#'
#' @return 
#' The \code{DianaParam} constructor will return a \linkS4class{DianaParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects}
#' (a list containing \code{diana}, the function output; \code{dist}, the dissimilarity matrix; and \code{hclust}, a \link{hclust} object created from \code{diana}).
#'
#' @examples
#' clusterRows(iris[,1:4], DianaParam())
#' clusterRows(iris[,1:4], DianaParam(metric="manhattan"))
#'
#' @seealso
#' \code{\link{diana}}, which actually does all the heavy lifting.
#'
#' \linkS4class{HclustParam}, for the more commonly used implementation of hierarchical clustering.
#'
#' @name DianaParam-class
#' @docType class
#' @aliases 
#' show,DianaParam-method
#' .defaultScalarArguments,DianaParam-method
NULL

#' @export
#' @rdname DianaParam-class
setClass("DianaParam", contains="HierarchicalParam", slots=c(stand="ANY"))

#' @export
setMethod(".defaultScalarArguments", "DianaParam", function(x) c(callNextMethod(), stand="logical"))

#' @export
#' @rdname DianaParam-class
DianaParam <- function(metric=NULL, stand=NULL, cut.fun=NULL, cut.dynamic=FALSE, cut.params=list()) {
    new("DianaParam", metric=metric, stand=stand, cut.fun=cut.fun, cut.dynamic=cut.dynamic, cut.params=cut.params)
}

#' @export
#' @rdname DianaParam-class
#' @importFrom cluster diana 
#' @importFrom stats as.hclust
setMethod("clusterRows", c("ANY", "DianaParam"), function(x, BLUSPARAM, full=FALSE) {
    args <- c(
        list(
            quote(as.matrix(x)),
            keep.diss = TRUE,
            keep.data = FALSE
        ),
        .extractScalarArguments(BLUSPARAM)
    )

    out <- do.call(diana, args)
    dst <- out$diss
    hcl <- as.hclust(out)

    clusters <- .cut_hierarchical(hcl, dst, BLUSPARAM)
    if (full) {
        list(clusters=clusters, objects=list(diana=out, dist=dst, hclust=hcl))
    } else {
        clusters
    }
})
