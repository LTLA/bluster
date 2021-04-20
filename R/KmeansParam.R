#' K-means clustering
#'
#' Run the base \code{\link{kmeans}} function with the specified number of centers within \code{\link{clusterRows}}.
#'
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers.
#' @param iter.max,nstart,algorithm Further arguments to pass to \code{\link{kmeans}}.
#' Set to the \code{\link{kmeans}} defaults if not supplied.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{KmeansParam} object.
#' @param full Logical scalar indicating whether the full k-means statistics should be returned.
#'
#' @author Aaron Lun
#'
#' @details
#' This class usually requires the user to specify the number of clusters beforehand.
#' However, we can also allow the number of clusters to vary as a function of the number of observations.
#' The latter is occasionally useful, e.g., to allow the clustering to automatically become more granular for large datasets.
#'
#' To modify an existing KmeansParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' @return 
#' The \code{KmeansParam} constructor will return a \linkS4class{KmeansParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects}
#' (a list containing \code{kmeans}, the direct output of \code{\link{kmeans}}).
#'
#' @examples
#' clusterRows(iris[,1:4], KmeansParam(centers=4))
#' clusterRows(iris[,1:4], KmeansParam(centers=4, algorithm="Lloyd"))
#' clusterRows(iris[,1:4], KmeansParam(centers=sqrt))
#' @seealso
#' \code{\link{kmeans}}, which actually does all the heavy lifting.
#'
#' \linkS4class{MbkmeansParam}, for a faster but more approximate version of the k-means algorithm.
#' @name KmeansParam-class
#' @docType class
#' @aliases 
#' show,KmeansParam-method
#' updateObject,KmeansParam-method
NULL

#' @export
#' @rdname KmeansParam-class
setClass("KmeansParam", contains="FixedNumberParam", slots=c(iter.max="integer", nstart="integer", algorithm="character"))

#' @export
#' @rdname KmeansParam-class
KmeansParam <- function(centers, iter.max = NULL, nstart = NULL, algorithm = NULL) {
    if (!is.function(centers)) {
        centers <- as.integer(centers)
    }

    # Filling in missing values with the defaults.
    current <- list(iter.max=iter.max, nstart=nstart, algorithm=algorithm)
    notpresent <- vapply(current, is.null, FALSE)
    if (any(notpresent)) {
        defaults <- .get_kmeans_defaults()
        current[notpresent] <- defaults[notpresent]
    }

    new("KmeansParam", centers=centers, iter.max=as.integer(current$iter.max), 
        nstart=as.integer(current$nstart), algorithm=current$algorithm)
}

.get_kmeans_defaults <- function() {
    args <- formals(kmeans)
    out <- args[c("iter.max", "nstart", "algorithm")]
    out$algorithm <- eval(out$algorithm)[1]
    out
}

#' @export
setMethod("show", "KmeansParam", function(object) {
    callNextMethod()
    cat(sprintf("iter.max: %s\n", object@iter.max))
    cat(sprintf("nstart: %s\n", object@nstart))
    cat(sprintf("algorithm: %s\n", object@algorithm))
})

#' @importFrom S4Vectors setValidity2
setValidity2("KmeansParam", function(object) {
    msg <- character(0)

    if (length(object@iter.max) > 1) {
        msg <- c(msg, "'iter.max' must be an integer scalar")
    }
    if (length(object@nstart) > 1) {
        msg <- c(msg, "'nstart' must be an integer scalar")
    }
    if (length(object@algorithm) > 1) {
        msg <- c(msg, "'algorithm' must be a string")
    }

    if (length(msg)) return(msg)
    TRUE
})

#' @export
#' @rdname KmeansParam-class
#' @importFrom stats kmeans
setMethod("clusterRows", c("ANY", "KmeansParam"), function(x, BLUSPARAM, full=FALSE) {
    centerx <- centers(BLUSPARAM, n=nrow(x))

    stats <- kmeans(as.matrix(x), centers=centerx, iter.max=BLUSPARAM@iter.max,
        nstart=BLUSPARAM@nstart, algorithm=BLUSPARAM@algorithm)
    clusters <- factor(stats$cluster)

    if (full) {
        list(clusters=clusters, objects=list(kmeans=stats))
    } else {
        clusters
    }
})
