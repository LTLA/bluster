#' FlowSOM-based clustering
#'
#' Use the self-organizing map implementation in the \pkg{FlowSOM} package to cluster observations into the specified number of nodes.
#' Note that this requires the installation of the \pkg{FlowSOM} package.
#'
#' @param centers An integer scalar specifying the number of centers.
#' Alternatively, a function that takes the number of observations and returns the number of centers.
#' @inheritParams clusterRows
#' @param dim.ratio A positive numeric scalar in specifying how \code{centers} should be distributed between the \code{x} and \code{y} dimensions.
#' Defaults to equal distribution, i.e., both dimensions will be of length equal to the square root of \code{centers}.
#' Values above 1 will distribute more nodes to \code{x} while values below 1 will distribute mode nodes to \code{y}.
#' @param rlen,mst,alpha,radius,init,initf,distf
#' Further arguments to pass to the \code{\link[FlowSOM]{SOM}} function in the \pkg{FlowSOM} package.
#' @param BLUSPARAM A \linkS4class{FlowSOMParam} object.
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
#' To modify an existing FlowSOMParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' For \code{radius}, a value of \code{NULL} means that the default argument in the \code{\link[FlowSOM]{SOM}} function signature is used.
#' This is are data-dependent and so cannot be specified during construction of the FlowSOMParam object.
#' \code{initf} has the same behavior.
#'
#' @return
#' The \code{FlowSOMParam} constructor will return a \linkS4class{FlowSOMParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects};
#' the latter will contain the direct output of \code{\link[FlowSOM]{SOM}}.
#'
#' @author Aaron Lun
#' @examples
#' clusterRows(iris[,1:4], FlowSOMParam(centers=16))
#' clusterRows(iris[,1:4], FlowSOMParam(centers=12, dim.ratio=3/4))
#'
#' @seealso
#' \code{\link[FlowSOM]{SOM}} from the \pkg{FlowSOM} package, which does all of the heavy lifting.
#'
#' \linkS4class{FixedNumberParam}, the parent of the FlowSOMParam class.
#' @name FlowSOMParam-class
#' @docType class
#' @aliases 
#' show,FlowSOMParam-method
NULL

#' @export
#' @rdname FlowSOMParam-class
setClass("FlowSOMParam", contains="FixedNumberParam", slots=c(
    dim.ratio="numeric", 
    rlen="integer", 
    mst="integer",
    alpha="numeric",
    radius="numeric_OR_NULL",
    init="logical",
    initf="function_OR_NULL",
    distf="integer"
))

#' @export
setMethod("show", "FlowSOMParam", function(object) {
    callNextMethod()
    cat(sprintf("dim.ratio: %s\n", object@dim.ratio))
    cat(sprintf("rlen: %s\n", object@rlen))
    cat(sprintf("mst: %s\n", object@mst))
    cat(sprintf("alpha: %s\n", paste(object@alpha, collapse=" ")))
    cat(sprintf("radius: %s\n", if (is.null(object@radius)) "default" else paste(object@radius, collapse=" ")))
    cat(sprintf("init: %s\n", object@init))
    cat(sprintf("initf: %s\n", if (is.null(object@initf)) "default" else "custom"))
    cat(sprintf("distf: %s\n", object@distf))
})

#' @export
#' @rdname FlowSOMParam-class
FlowSOMParam <- function(centers, 
    dim.ratio = 1,
    rlen = 10,
    mst = 1,
    alpha = c(0.05, 0.01),
    radius = NULL,
    init = FALSE,
    initf = NULL,
    distf = 2)
{
    if (!is.function(centers)) {
        centers <- as.integer(centers)
    }
    new("FlowSOMParam", 
        centers=centers, 
        dim.ratio=dim.ratio,
        rlen=as.integer(rlen),
        mst=as.integer(mst),
        alpha=alpha,
        radius=radius,
        init=init,
        initf=initf,
        distf=as.integer(distf)
    )
}

#' @export
#' @rdname FlowSOMParam-class
setMethod("clusterRows", c("ANY", "FlowSOMParam"), function(x, BLUSPARAM, full=FALSE) {
    centerx <- centers(BLUSPARAM, n=nrow(x))
    dim.ratio <- BLUSPARAM[["dim.ratio"]]
    xdim <- max(1, round(sqrt(centerx * dim.ratio)))
    ydim <- max(1, round(sqrt(centerx / dim.ratio)))

    args <- list(
        data = as.matrix(x),
        xdim = xdim,
        ydim = ydim,
        rlen = BLUSPARAM[["rlen"]],
        mst = BLUSPARAM[["mst"]],
        alpha = BLUSPARAM[["alpha"]],
        init = BLUSPARAM[["init"]],
        initf = BLUSPARAM[["initf"]],
        radius = BLUSPARAM[["radius"]],
        distf = BLUSPARAM[["distf"]],
        codes = NULL
    )

#    if (!is.null(BLUSPARAM[["initf"]])) {
#        args$initf <- BLUSPARAM[["initf"]]
#    }
#    if (!is.null(BLUSPARAM[["radius"]])) {
#        args$radius <- BLUSPARAM[["radius"]]
#    }

    stats <- do.call(.SOM, args)
    clusters <- factor(stats$mapping)

    if (full) {
        list(clusters=clusters, objects=stats)
    } else {
        clusters 
    }
})

#' @importFrom stats quantile dist
.SOM <- function(data, xdim, ydim, rlen, mst, alpha, radius, codes, init, initf, distf) {
    # Initialize the grid
    grid <- expand.grid(seq_len(xdim),seq_len(ydim))
    nCodes <- nrow(grid)
    if (is.null(codes)){
        if (init) {
            stopifnot(is.function(initf))
            codes <- initf(data, xdim, ydim)
        } else {
            codes <- data[sample(nrow(data), nCodes),, drop = FALSE]
        }
    }

    # Initialize the neighbourhood
    nhbrdist <- dist(grid, method = "maximum")

    if (is.null(radius)) {
        radius <- c(quantile(nhbrdist, 0.67), 0)
    }

    # Initialize the radius
    if (mst==1) {
        radius <- list(radius)
        alpha <- list(alpha)
    } else {
        radius <- seq(radius[1], radius[2], length.out = mst+1)
        radius <- lapply(seq_len(mst), function(i) c(radius[i], radius[i+1]))
        alpha <- seq(alpha[1], alpha[2], length.out = mst+1)
        alpha <- lapply(seq_len(mst), function(i) c(alpha[i], alpha[i+1]))
    }

    # Compute the SOM
    data <- t(as.matrix(data))
    codes <- t(as.matrix(codes))

    for (i in seq_len(mst)) {
        res <- flow_som(data = data,
            original_codes = codes,
            nhbrdist = as.matrix(nhbrdist),
            alpha = alpha[[i]],
            radii = radius[[i]],
            rlen = rlen,
            dist = distf)

        codes <- res[[1]]
        nhbrdist <- distance.along.mst(codes)
    }

    names(res) <- c("codes", "mapping")
    rownames(res$codes) <- rownames(data)
    res$codes <- t(res$codes)
    res$mapping <- res$mapping + 1L

    res
}

#' @importFrom stats dist
#' @importFrom igraph graph.adjacency minimum.spanning.tree shortest.paths V
distance.along.mst <- function(X) {
    adjacency <- dist(t(X), method = "euclidean")
    fullGraph <- graph.adjacency(as.matrix(adjacency), mode = "undirected", weighted = TRUE)
    mst <- minimum.spanning.tree(fullGraph)
    shortest.paths(mst, v=V(mst), to=V(mst), weights=NA)
}
