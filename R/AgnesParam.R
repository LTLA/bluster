#' Agglomerative nesting 
#'
#' Run the \code{\link{agnes}} function on a distance matrix within \code{\link{clusterRows}}.
#'
#' @inheritParams HclustParam
#' @param metric,stand,method,par.method Further arguments to pass to \code{\link{agnes}}.
#'
#' @author Aaron Lun
#'
#' @details
#' To modify an existing AgnesParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' If \code{cut.fun=NULL}, \code{cut.dynamic=FALSE} and \code{cut.params} does not have \code{h} or \code{k},
#' \code{\link{clusterRows}} will automatically set \code{h} to half the tree height when calling \code{\link{cutree}}. 
#'
#' @return 
#' The \code{AgnesParam} constructor will return a \linkS4class{AgnesParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects}
#' (a list containing \code{agnes}, the function output; \code{dist}, the dissimilarity matrix; and \code{hclust}, a \link{hclust} object created from \code{agnes}).
#'
#' @examples
#' clusterRows(iris[,1:4], AgnesParam())
#' clusterRows(iris[,1:4], AgnesParam(method="ward"))
#'
#' @seealso
#' \code{\link{agnes}}, which actually does all the heavy lifting.
#'
#' \linkS4class{HclustParam}, for the more commonly used implementation of hierarchical clustering.
#'
#' @name AgnesParam-class
#' @docType class
#' @aliases 
#' show,AgnesParam-method
#' .defaultScalarArguments,AgnesParam-method
NULL

#' @export
#' @rdname AgnesParam-class
setClass("AgnesParam", contains="HierarchicalParam", slots=c(stand="ANY", method="ANY", par.method="ANY"))

#' @export
setMethod(".defaultScalarArguments", "AgnesParam", function(x) c(callNextMethod(), stand="logical", method="character", par.method="numeric")) 

#' @export
#' @rdname AgnesParam-class
AgnesParam <- function(metric=NULL, stand=NULL, method=NULL, par.method=NULL, cut.fun=NULL, cut.dynamic=FALSE, cut.params=list()) {
    new("AgnesParam", metric=metric, stand=stand, method=method, par.method=par.method, cut.fun=cut.fun, cut.dynamic=cut.dynamic, cut.params=cut.params)
}

#' @export
#' @rdname AgnesParam-class
#' @importFrom cluster agnes
#' @importFrom stats as.hclust
setMethod("clusterRows", c("ANY", "AgnesParam"), function(x, BLUSPARAM, full=FALSE) {
    args <- c(
        list(
            quote(as.matrix(x)),
            keep.diss = TRUE,
            keep.data = FALSE
        ),
        .extractScalarArguments(BLUSPARAM)
    )

    out <- do.call(agnes, args)
    dst <- out$diss
    hcl <- as.hclust(out)

    clusters <- .cut_hierarchical(hcl, dst, BLUSPARAM)
    if (full) {
        list(clusters=clusters, objects=list(agnes=out, dist=dst, hclust=hcl))
    } else {
        clusters
    }
})
