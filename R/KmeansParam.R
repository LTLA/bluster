#' K-means clustering
#'
#' Run the base \code{\link{kmeans}} function with the specified number of centers within \code{\link{clusterRows}}.
#'
#' @param centers Integer scalar specifying the number of centers.
#' @param ... Further arguments to pass to \code{\link{kmeans}}.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{KmeansParam} object.
#' @param full Logical scalar indicating whether the full k-means statistics should be returned.
#'
#' @author Aaron Lun
#'
#' @return 
#' The \code{KmeansParam} constructor will return a \linkS4class{KmeansParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects};
#' the latter will contain the direct output of \code{\link{kmeans}}.
#'
#' @examples
#' clusterRows(iris[,1:4], KmeansParam(centers=4))
#' clusterRows(iris[,1:4], KmeansParam(centers=4, algorithm="Lloyd"))
#'
#' @seealso
#' \code{\link{kmeans}}, which actually does all the heavy lifting.
#' @name KmeansParam-class
NULL

#' @export
#' @rdname KmeansParam-class
setClass("KmeansParam", contains="BlusterParam", slots=c(centers="integer", extra.args="list"))

#' @export
#' @rdname KmeansParam-class
KmeansParam <- function(centers, ...) {
    new("KmeansParam", centers=as.integer(centers), extra.args=list(...))
}

#' @export
#' @rdname KmeansParam-class
#' @importFrom stats kmeans
setMethod("clusterRows", c("ANY", "KmeansParam"), function(x, BLUSPARAM, full=FALSE) {
    args <- c(list(x=as.matrix(x), centers=BLUSPARAM@centers), BLUSPARAM@extra.args)
    stats <- do.call(kmeans, args)
    clusters <- factor(stats$cluster)

    if (full) {
        list(clusters=clusters, objects=stats)
    } else {
        clusters
    }
})
