#' Mini-batch k-means clustering
#'
#' Run the mini-batch k-means \code{\link[mbkmeans]{mbkmeans}} function with the specified number of centers within \code{\link{clusterRows}}.
#' This sacrifices some accuracy for speed compared to the standard k-means algorithm.
#' Note that this requires installation of the \pkg{mbkmeans} package.
#'
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers.
#' @inheritParams clusterRows
#' @param batch_size,max_iters,num_init,init_fraction,initializer,calc_wcss,early_stop_iter,tol,BPPARAM
#' Further arguments to pass to \code{\link[mbkmeans]{mbkmeans}}.
#' @param BLUSPARAM A \linkS4class{MbkmeansParam} object.
#' @param full Logical scalar indicating whether the full mini-batch k-means statistics should be returned.
#'
#' @author Stephanie Hicks
#'
#' @details
#' This class usually requires the user to specify the number of clusters beforehand.
#' However, we can also allow the number of clusters to vary as a function of the number of observations.
#' The latter is occasionally useful, e.g., to allow the clustering to automatically become more granular for large datasets.
#'
#' To modify an existing MbkmeansParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' For \code{batch_size} and \code{init_fraction}, a value of \code{NULL} means that the default arguments in the \code{\link[mbkmeans]{mbkmeans}} function signature are used.
#' These defaults are data-dependent and so cannot be specified during construction of the MbkmeansParam object, but instead are defined within the \code{clusterRows} method.
#'
#' @return 
#' The \code{MbkmeansParam} constructor will return a \linkS4class{MbkmeansParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects}
#' (a list containing \code{mbkmeans}, the direct output of \code{\link[mbkmeans]{mbkmeans}}).
#'
#' @examples
#' clusterRows(iris[,1:4], MbkmeansParam(centers=3))
#' clusterRows(iris[,1:4], MbkmeansParam(centers=3, batch_size=10))
#' clusterRows(iris[,1:4], MbkmeansParam(centers=3, init_fraction=0.5))
#' @seealso
#' \code{\link[mbkmeans]{mbkmeans}} from the \pkg{mbkmeans} package, which actually does all the heavy lifting.
#'
#' \linkS4class{KmeansParam}, for dispatch to the standard k-means algorithm.
#'
#' @name MbkmeansParam-class
#' @docType class
#' @aliases 
#' show,MbkmeansParam-method
NULL

#' @export
#' @rdname MbkmeansParam-class
setClass("MbkmeansParam", contains="FixedNumberParam", 
    slots=c(batch_size="integer_OR_NULL", 
        max_iters="integer", 
        num_init="integer", 
        init_fraction="numeric_OR_NULL", 
        initializer="character",
        calc_wcss="logical", 
        early_stop_iter="integer", 
        tol="numeric", 
        BPPARAM="BiocParallelParam"))

setValidity2("MbkmeansParam", function(object) {
    msg <- character(0)

    msg <- c(msg, .check_positive_slots(object, c("batch_size", "max_iters", "init_fraction", "num_init", "early_stop_iter", "tol")))

    if (!is.null(val <- object@init_fraction)) {
        if (val > 1) {
            msg <- c(msg, "'init_fraction' cannot be greater than 1")
        }
    }

    msg <- c(msg, .check_nonna_slots(object, c("initializer", "calc_wcss")))

    if (length(msg)) {
        return(msg)
    }
    TRUE
})

#' @export
setMethod("show", "MbkmeansParam", function(object) {
    callNextMethod()
    cat(sprintf("batch_size: %s\n", if (is.null(object@batch_size)) "default" else object@batch_size)) 
    cat(sprintf("max_iters: %s\n", object@max_iters))
    cat(sprintf("num_init: %s\n", object@num_init))
    cat(sprintf("init_fraction: %s\n", if (is.null(object@init_fraction)) "default" else object@init_fraction)) 
    cat(sprintf("initializer: %s\n", object@initializer))
    cat(sprintf("calc_wcss: %s\n", object@calc_wcss))
    cat(sprintf("early_stop_iter: %s\n", object@early_stop_iter))
    cat(sprintf("tol: %s\n", object@tol))
    cat(sprintf("BPPARAM: %s\n", class(object@BPPARAM)[1]))
})

#' @export
#' @rdname MbkmeansParam-class
#' @importFrom BiocParallel SerialParam
MbkmeansParam <- function(centers, 
    batch_size = NULL,
    max_iters = 100,
    num_init = 1,
    init_fraction = NULL,
    initializer = "kmeans++",
    calc_wcss = FALSE,
    early_stop_iter = 10,
    tol = 1e-04,
    BPPARAM = SerialParam())
{
    if (!is.function(centers)) {
        centers <- as.integer(centers)
    }
    if (!is.null(batch_size)) {
        batch_size <- as.integer(batch_size)
    }

    new("MbkmeansParam", centers=centers, 
        batch_size=batch_size, 
        max_iters=as.integer(max_iters), 
        num_init=as.integer(num_init),
        init_fraction=init_fraction, 
        initializer=initializer, 
        calc_wcss=calc_wcss,
        early_stop_iter=as.integer(early_stop_iter), 
        tol=tol, 
        BPPARAM=BPPARAM)
}

#' @export
#' @rdname MbkmeansParam-class
setMethod("clusterRows", c("ANY", "MbkmeansParam"), function(x, BLUSPARAM, full=FALSE) {
    centerx <- centers(BLUSPARAM, n=nrow(x))

    # Setting the data-dependent defaults.
    batch_size <- BLUSPARAM[["batch_size"]]
    if (is.null(batch_size)) {
        batch_size <- min(500, nrow(x))
    }

    init_fraction <- BLUSPARAM[["init_fraction"]]
    if (is.null(init_fraction)) {
        init_fraction <- batch_size/nrow(x)
    }

    if (is.matrix(x) && is.object(x)) {
        # their C++ code can't handle S3 classes that inherit from base matrix objects.
        x <- unclass(x)
    }

    stats <- mbkmeans::mbkmeans(t(x), clusters=centerx,
        batch_size=batch_size,
        max_iters=BLUSPARAM[["max_iters"]], 
        num_init=BLUSPARAM[["num_init"]],
        init_fraction=init_fraction,
        initializer=BLUSPARAM[["initializer"]], 
        calc_wcss=BLUSPARAM[["calc_wcss"]],
        early_stop_iter=BLUSPARAM[["early_stop_iter"]],
        tol=BLUSPARAM[["tol"]], 
        BPPARAM=BLUSPARAM[["BPPARAM"]])

    vec_clusters <- factor(stats$Clusters)
    if (full) {
        list(clusters=vec_clusters, objects=list(mbkmeans=stats))
    } else {
        vec_clusters
    }
})
