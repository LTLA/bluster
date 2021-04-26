#' Clustering Large Applications
#'
#' Run the CLARA algorithm, an extension of the PAM method for large datasets.
#'
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers.
#' @param metric,stand,samples,sampsize Further arguments to pass to \code{\link{clara}}.
#' Set to the function defaults if not supplied.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{ClaraParam} object.
#' @param full Logical scalar indicating whether the full PAM statistics should be returned.
#'
#' @author Aaron Lun
#'
#' @details
#' This class usually requires the user to specify the number of clusters beforehand.
#' However, we can also allow the number of clusters to vary as a function of the number of observations.
#' The latter is occasionally useful, e.g., to allow the clustering to automatically become more granular for large datasets.
#'
#' To modify an existing ClaraParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' Note that \code{clusterRows} will always use \code{rngR=TRUE}, for greater consistency with other algorithms of the \linkS4class{FixedNumberParam} class;
#' and \code{pamLike=TRUE}, for consistency with the PAM implementation from which it was derived.
#'
#' @return
#' The \code{ClaraParam} constructor will return a \linkS4class{ClaraParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects}
#' (a list containing \code{clara}, the direct output of \code{\link{clara}}).
#'
#' @examples
#' clusterRows(iris[,1:4], ClaraParam(centers=4))
#' clusterRows(iris[,1:4], ClaraParam(centers=4, sampsize=50))
#' clusterRows(iris[,1:4], ClaraParam(centers=sqrt))
#' @seealso
#' \code{\link{clara}}, which actually does all the heavy lifting.
#'
#' \linkS4class{PamParam}, for the original PAM algorithm.
#'
#' @name ClaraParam-class
#' @docType class
#' @aliases
#' show,ClaraParam-method
#' .defaultScalarArguments,ClaraParam-method
NULL

#' @export
#' @rdname ClaraParam-class
setClass("ClaraParam", contains="FixedNumberParam", slots=c(metric="ANY", stand="ANY", samples="ANY", sampsize="ANY"))

#' @export
#' @rdname ClaraParam-class
ClaraParam <- function(centers, metric=NULL, stand=NULL, samples=NULL, sampsize=NULL) {
    if (!is.function(centers)) {
        centers <- as.integer(centers)
    }
    new("ClaraParam", centers=centers, metric=metric, stand=stand, samples=samples, sampsize=sampsize)
}

#' @export
setMethod("show", "ClaraParam", function(object) {
    callNextMethod()
    .showScalarArguments(object)
})

#' @export
setMethod(".defaultScalarArguments", "ClaraParam", function(x) 
    c(metric="character", stand="logical", samples="numeric", sampsize="numeric"))

#' @export
#' @rdname ClaraParam-class
#' @importFrom cluster clara
setMethod("clusterRows", c("ANY", "ClaraParam"), function(x, BLUSPARAM, full=FALSE) {
    centerx <- centers(BLUSPARAM, n=nrow(x))

    args <- c(
        list(
            quote(as.matrix(x)),
            k=centerx,
            keep.data=FALSE,
            rngR=TRUE,
            pamLike=TRUE
        ),
        .extractScalarArguments(BLUSPARAM)
    )

    stats <- do.call(clara, args)
    clusters <- factor(stats$clustering)

    if (full) {
        list(clusters=clusters, objects=list(clara=stats))
    } else {
        clusters
    }
})
