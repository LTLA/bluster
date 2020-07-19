#' K-means clustering
#'
#' Run the base \code{\link{kmeans}} function with the specified number of centers within \code{\link{clusterRows}}.
#'
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers.
#' @param ... Further arguments to pass to \code{\link{kmeans}}.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{KmeansParam} object.
#' @param full Logical scalar indicating whether the full k-means statistics should be returned.
#'
#' @author Aaron Lun
#'
#' @details
#' The standard \linkS4class{KmeansParam} class requires the user to specify the number of clusters beforehand.
#' However, we can also allow the number of clusters to vary as a function of the number of observations.
#' The latter is occasionally useful, e.g., to allow the clustering to automatically become more granular for large datasets.
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
#' clusterRows(iris[,1:4], KmeansParam(centers=sqrt))
#' @seealso
#' \code{\link{kmeans}}, which actually does all the heavy lifting.
#' @name KmeansParam-class
#' @docType class
#' @aliases 
#' show,KmeansParam-method
NULL

#' @export
#' @rdname KmeansParam-class
setClass("KmeansParam", contains="BlusterParam", slots=c(centers="integer_OR_function", extra.args="list"))

#' @export
#' @rdname KmeansParam-class
KmeansParam <- function(centers, ...) {
    if (!is.function(centers)) {
        centers <- as.integer(centers)
    }
    new("KmeansParam", centers=centers, extra.args=list(...))
}

setMethod(".extras", "KmeansParam", function(x) "extra.args")

#' @export
setMethod("show", "KmeansParam", function(object) {
    callNextMethod()
    cat(sprintf("centers: %s\n", if (is.function(object@centers)) "variable" else object@centers))
    coolcat("extra.args(%i): %s", names(object@extra.args))
})

#' @importFrom S4Vectors setValidity2
setValidity2("KmeansParam", function(object) {
    msg <- character(0)

    val <- object@centers
    if (!is.function(val) && !.positive_number(val)) {
        msg <- c(msg, "'centers' must be a positive number")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
#' @rdname KmeansParam-class
#' @importFrom stats kmeans
setMethod("clusterRows", c("ANY", "KmeansParam"), function(x, BLUSPARAM, full=FALSE) {
    centers <- BLUSPARAM@centers
    if (is.function(centers)) {
        centers <- centers(nrow(x))
    }

    args <- c(list(x=as.matrix(x), centers=centers), BLUSPARAM@extra.args)
    stats <- do.call(kmeans, args)
    clusters <- factor(stats$cluster)

    if (full) {
        list(clusters=clusters, objects=stats)
    } else {
        clusters
    }
})
