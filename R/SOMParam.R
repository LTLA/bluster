#' SOM-based clustering
#'
#' Use the self-organizing map implementation in the \pkg{kohonen} package to cluster observations into the specified number of nodes.
#' Note that this requires the installation of the \pkg{kohonen} package.
#'
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers.
#' @inheritParams clusterRows
#' @param dim.ratio A positive numeric scalar in specifying how \code{centers} should be distributed between the \code{x} and \code{y} dimensions.
#' Defaults to equal distribution, i.e., both dimensions will be of length equal to the square root of \code{centers}.
#' Values above 1 will distribute more nodes to \code{x} while values below 1 will distribute mode nodes to \code{y}.
#' @param topo,neighbourhood.fct,toroidal
#' Further arguments to pass to the \code{\link[kohonen]{somgrid}} function in the \pkg{kohonen} package.
#' @param rlen,alpha,radius,dist.fct
#' Further arguments to pass to the \code{\link[kohonen]{som}} function in the \pkg{kohonen} package.
#' @param BLUSPARAM A \linkS4class{SOMParam} object.
#' @param full Logical scalar indicating whether the full mini-batch k-means statistics should be returned.
#' 
#' @details
#' This class usually requires the user to specify the number of clusters beforehand.
#' However, we can also allow the number of clusters to vary as a function of the number of observations.
#' The latter is occasionally useful, e.g., to allow the clustering to automatically become more granular for large datasets.
#'
#' Note that the final number of clusters may not be exactly equal to \code{centers}, depending on how \code{dim.ratio} is specified.
#' For example, if \code{centers} is a perfect square and \code{dim.ratio=1}, we will get exactly the requested number of points.
#'
#' To modify an existing SOMParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' For \code{radius}, a value of \code{NULL} means that the default argument in the \code{\link[kohonen]{som}} function signature is used.
#' This is are data-dependent and so cannot be specified during construction of the SOMParam object.
#'
#' For \code{dist.fct}, users can specify any string that can be used in the \code{dist.fcts} arguments in \code{\link[kohonen]{som}}.
#' In practice, the only real alternative is \code{"manhattan"}.
#'
#' @return
#' The \code{SOMParam} constructor will return a \linkS4class{SOMParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects};
#' the latter will contain the direct output of \code{\link[kohonen]{som}}.
#'
#' @author Aaron Lun
#' @examples
#' clusterRows(iris[,1:4], SOMParam(centers=16))
#' clusterRows(iris[,1:4], SOMParam(centers=12, dim.ratio=3/4))
#'
#' @seealso
#' \code{\link[kohonen]{som}} from the \pkg{kohonen} package, which does all of the heavy lifting.
#'
#' \linkS4class{FixedNumberParam}, the parent of the SOMParam class.
#' @name SOMParam-class
#' @docType class
#' @aliases 
#' show,SOMParam-method
NULL

#' @export
#' @rdname SOMParam-class
setClass("SOMParam", contains="FixedNumberParam", slots=c(
    dim.ratio="numeric", 
    topo="character",
    neighbourhood.fct="character",
    toroidal="logical",
    rlen="integer", 
    alpha="numeric",
    radius="numeric_OR_NULL",
    dist.fct="character"
))

setValidity2("SOMParam", function(object) {
    msg <- character(0)

    msg <- c(msg, .check_positive_slots(object, c("dim.ratio", "rlen"))) 

    val <- object@alpha
    if (length(val)!=2 || !all(is.finite(val)) || !all(val > 0)) {
        msg <- c(msg, "'alpha' must be a positive numeric vector of length 2")
    }

    val <- object@radius
    if (!is.null(val)) {
        if (length(val)!=2 || !all(is.finite(val)) || !all(val > 0)) {
            msg <- c(msg, "'radius' must be a positive numeric vector of length 2")
        }
    }

    msg <- c(msg, .check_nonna_slots(object, c("dist.fct", "topo", "neighbourhood.fct", "toroidal")))

    if (length(msg)) {
        return(msg)
    }
    TRUE
})

#' @export
setMethod("show", "SOMParam", function(object) {
    callNextMethod()

    for (i in c("dim.ratio", "topo", "neighbourhood.fct", "toroidal", "rlen")) {
        cat(sprintf("%s: %s\n", i, slot(object, i)))
    }

    cat(sprintf("alpha: %s\n", paste(object@alpha, collapse=" ")))
    cat(sprintf("radius: %s\n", if (is.null(object@radius)) "default" else paste(object@radius, collapse=" ")))
    cat(sprintf("dist.fct: %s\n", object@dist.fct))
})

#' @export
#' @rdname SOMParam-class
SOMParam <- function(centers, 
    dim.ratio = 1,
    topo = "rectangular",
    neighbourhood.fct = "bubble",
    toroidal = FALSE,
    rlen = 100,
    alpha = c(0.05, 0.01),
    radius = NULL,
    dist.fct = "sumofsquares")
{
    if (!is.function(centers)) {
        centers <- as.integer(centers)
    }
    new("SOMParam", 
        centers=centers, 
        dim.ratio=dim.ratio,
        topo=topo,
        neighbourhood.fct=neighbourhood.fct,
        toroidal=toroidal,
        rlen=as.integer(rlen),
        alpha=alpha,
        radius=radius,
        dist.fct=dist.fct
    )
}

#' @export
#' @rdname SOMParam-class
setMethod("clusterRows", c("ANY", "SOMParam"), function(x, BLUSPARAM, full=FALSE) {
    centerx <- centers(BLUSPARAM, n=nrow(x))
    dim.ratio <- BLUSPARAM[["dim.ratio"]]
    xdim <- max(1, round(sqrt(centerx * dim.ratio)))
    ydim <- max(1, round(sqrt(centerx / dim.ratio)))

    args <- list(
        X = as.matrix(x),
        grid = kohonen::somgrid(
            xdim = xdim,
            ydim = ydim,
            topo = BLUSPARAM@topo,
            neighbourhood.fct = BLUSPARAM@neighbourhood.fct,
            toroidal = BLUSPARAM@toroidal
        ),
        rlen = BLUSPARAM@rlen,
        alpha = BLUSPARAM@alpha,
        dist.fcts = list(BLUSPARAM@dist.fct)
    )

    if (!is.null(BLUSPARAM[["radius"]])) {
        args$radius <- BLUSPARAM[["radius"]]
    }

    stats <- do.call(kohonen::som, args)
    clusters <- factor(stats$unit.classif)

    if (full) {
        list(clusters=clusters, objects=stats)
    } else {
        clusters 
    }
})
