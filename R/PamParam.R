#' Partitioning around medoids
#'
#' Partition observations into k-medoids as a more robust version of k-means.
#'
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers.
#' @param metric,medoids,nstart,stand,do.swap,variant Further arguments to pass to \code{\link[cluster]{pam}}.
#' Set to the function defaults if not supplied.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{PamParam} object.
#' @param full Logical scalar indicating whether the full PAM statistics should be returned.
#'
#' @author Aaron Lun
#'
#' @details
#' This class usually requires the user to specify the number of clusters beforehand.
#' However, we can also allow the number of clusters to vary as a function of the number of observations.
#' The latter is occasionally useful, e.g., to allow the clustering to automatically become more granular for large datasets.
#'
#' To modify an existing PamParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' @return
#' The \code{PamParam} constructor will return a \linkS4class{PamParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects}
#' (a list containing \code{pam}, the direct output of \code{\link[cluster]{pam}}).
#'
#' @examples
#' clusterRows(iris[,1:4], PamParam(centers=4))
#' clusterRows(iris[,1:4], PamParam(centers=4, variant="faster", do.swap=FALSE))
#' clusterRows(iris[,1:4], PamParam(centers=sqrt))
#' @seealso
#' \code{\link[pam]{pam}}, which actually does all the heavy lifting.
#'
#' \linkS4class{KmeansParam}, for the more commonly used k-means algorithm.
#' @name PamParam-class
#' @docType class
#' @aliases
#' show,PamParam-method
NULL

#' @export
#' @rdname PamParam-class
setClass("PamParam", contains="FixedNumberParam", slots=c(metric="ANY", medoids="ANY", nstart="ANY", stand="ANY", do.swap="ANY", variant="ANY"))

#' @export
#' @rdname PamParam-class
PamParam <- function(centers, metric=NULL, medoids=NULL, nstart=NULL, stand=NULL, do.swap=NULL, variant=NULL) {
    if (!is.function(centers)) {
        centers <- as.integer(centers)
    }
    new("PamParam", centers=centers, metric=metric, medoids=medoids, nstart=nstart, stand=stand, do.swap=do.swap, variant=variant)
}

pam.args <- c(metric="character", medoids="character", nstart="numeric", stand="logical", do.swap="logical", variant="character")

#' @export
setMethod("show", "PamParam", function(object) {
    callNextMethod()
    for (x in names(pam.args)) {
        val <- slot(object, x)
        if (is.null(val)) {
            val <- "[default]"
        }
        cat(sprintf("%s: %s\n", x, val))
    }
})

#' @importFrom S4Vectors setValidity2
setValidity2("PamParam", function(object) {
    for (x in names(pam.args)) {
        val <- slot(object, x)
        if (!is.null(val)) {
            if (length(val)!=1 || !is(val, pam.args[[x]])) {
                return(sprintf("'%s' should be NULL or a %s scalar", x, pam.args[[x]]))
            }
        }
    }
    TRUE
})

#' @export
#' @rdname PamParam-class
#' @importFrom stats kmeans
setMethod("clusterRows", c("ANY", "PamParam"), function(x, BLUSPARAM, full=FALSE) {
    centerx <- centers(BLUSPARAM, n=nrow(x))

    args <- list(
        as.matrix(x), 
        k=centerx,
        keep.data=FALSE
    )

    for (i in names(pam.args)) {
        val <- slot(BLUSPARAM, i)
        if (!is.null(val)) {
            args[[i]] <- val
        }
    }

    stats <- do.call(cluster::pam, args)
    clusters <- factor(stats$clustering)

    if (full) {
        list(clusters=clusters, objects=list(pam=stats))
    } else {
        clusters
    }
})
