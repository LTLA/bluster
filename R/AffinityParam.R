#' Affinity propogation 
#'
#' Use affinity propagation from the \pkg{apcluster} package to cluster observations.
#' Note that this requires the installation of the \pkg{apcluster} package.
#'
#' @inheritParams clusterRows
#' @param s A function that accepts a matrix of observations by dimensions and returns a similarity matrix.
#' If \code{NULL}, defaults to the output of \code{\link[apcluster]{negDistMat}} with \code{r=2}.
#' @param p,q Numeric scalars controlling the input preference, i.e., the resolution of the clustering.
#' These are passed to the \code{\link[apcluster]{apcluster}} function, where values of \code{NA} are the default.
#' @param maxits,convits,lam,nonoise
#' Further arguments to pass to the \code{\link[apcluster]{apcluster}} function.
#' @param BLUSPARAM A \linkS4class{AffinityParam} object.
#' @param full Logical scalar indicating whether the full affinity propagation statistics should be returned.
#' 
#' @details
#' To modify an existing AffinityParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' Setting \code{q} (and less typically, \code{p}) allows us to tune the resolution of the clustering.
#' In particular, when \code{p=NA}, it is computed based on the setting of \code{q}:
#' \itemize{
#' \item If the specified \code{q} lies in [0, 1], \code{p} is defined as the \code{q}-quantile of the finite similarities across all pairs of observations.
#' When \code{q=NA}, it defaults to 0.5.
#' \item If \code{q} is negative, \code{p} is defined as the \code{M + abs(M) * q} where \code{M} is the smallest finite similarity across all pairs.
#' This yields smaller \code{p} values while still responding to the scale of the similarities. 
#' }
#' The resulting value is used as the self-preference, i.e., the diagonal of the availability matrix.
#' Larger values yield more clusters as each data point is more inclined to form its own cluster.
#' 
#' @return
#' The \code{AffinityParam} constructor will return a \linkS4class{AffinityParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects}
#' (a list containing \code{similarity}, the similarity matrix; and \code{apcluster}, the direct output of \code{\link[apcluster]{apcluster}}).
#'
#' @author Aaron Lun
#' @examples
#' clusterRows(iris[,1:4], AffinityParam())
#' clusterRows(iris[,1:4], AffinityParam(q=0.9))
#' clusterRows(iris[,1:4], AffinityParam(s=apcluster::expSimMat()))
#'
#' @seealso
#' \code{\link[apcluster]{apcluster}} from the \pkg{apcluster} package, which does all of the heavy lifting.
#'
#' @name AffinityParam-class
#' @docType class
#' @aliases 
#' show,AffinityParam-method
NULL

#' @export
setClass("AffinityParam", contains="BlusterParam", slots=c(s="function_OR_NULL",
    p="numeric", q="numeric", maxits="integer", convits="integer", lam="numeric", nonoise="logical"))

#' @export
#' @rdname AffinityParam-class
AffinityParam <- function(s=NULL, p=NA, q=NA, maxits=1000, convits=100, lam=0.9, nonoise=FALSE) {
    new("AffinityParam", s=s, p=as.numeric(p), q=as.numeric(q), 
        maxits=as.integer(maxits), convits=as.integer(convits), 
        lam=as.numeric(lam), nonoise=as.logical(nonoise))
}

#' @export
setMethod("show", "AffinityParam", function(object) {
    callNextMethod()

    cat(sprintf("s: %s\n", if (is.null(object[["s"]])) "default" else "custom"))
    for (i in c("p", "q", "maxits", "convits", "lam", "nonoise")) {
        cat(sprintf("%s: %s\n", i, slot(object, i)))
    }
})

setValidity2("AffinityParam", function(object) {
    msg <- .check_positive_slots(object, c("maxits", "convits"))
    
    msg <- .check_nonna_slots(object, c("nonoise", "lam"))

    if (object@lam >= 1 || object@lam < 0.5) {
        msg <- c(msg, "'lam' should be in [0.5, 1)")
    }

    if (length(object[["p"]])!=1) {
        msg <- c(msg, "'p' should be of length 1")
    }

    q <- object[["q"]]
    if (length(q)!=1) {
        msg <- c(msg, "'q' should be of length 1")
    } else {
        if (!is.na(q) && q > 1) {
            msg <- c("'q' should lie in (-Inf, 1]")
        }
    }

    if (object@lam >= 1 || object@lam < 0.5) {
        msg <- c(msg, "'lam' should be in [0.5, 1)")
    }

    if (length(msg)) {
        return(msg)
    }
    TRUE
})

#' @export
#' @rdname AffinityParam-class
setMethod("clusterRows", c("ANY", "AffinityParam"), function(x, BLUSPARAM, full=FALSE) {
    s <- BLUSPARAM[["s"]]
    if (is.null(s)) {
        s <- apcluster::negDistMat(r=2)
    }
    mat <- s(as.matrix(x))

    p <- BLUSPARAM[["p"]]
    q <- BLUSPARAM[["q"]]
    if (is.na(p)) {
        if (!is.na(q) && q < 0) {
            keep <- mat > -Inf
            diag(keep) <- FALSE
            p <- min(mat[keep])
            p <- p + abs(p) * q
            q <- NA
        }
    }

    res <- apcluster::apcluster(s=mat, 
        p=p,
        q=q,
        maxits=BLUSPARAM[["maxits"]],
        convits=BLUSPARAM[["convits"]],
        lam=BLUSPARAM[["lam"]],
        nonoise=BLUSPARAM[["nonoise"]]
    )

    out <- apcluster::labels(res, type="enum")
    out <- factor(out)

    if (full) {
        list(clusters=out, objects=list(similarity=mat, apcluster=res))
    } else {
        out
    }
})
