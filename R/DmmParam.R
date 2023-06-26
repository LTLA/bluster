#' Dirichlet Multinomial Mixtures clustering
#' 
#' Use the DMM algorithm from the \code{DirichletMultinomial} package on a 
#' distance matrix within \code{\link{clusterRows}}.
#' 
#' @param k An integer vector indicating the number of clusters to choose from (or the possible number of clusters if the vector has multiple values)
#' @param type A string specifying the fit function to use to find the optimal number of clusters. 
#'   Use in combination with \code{k} being a list of integers. 
#'   Must be equal to \code{"laplace"}, \code{"AIC"} or \code{"BIC"}.
#' @param transposed Logical scalar, is x transposed with samples in rows?
#' @param seed An integer specifying the seed to use when doing the DMM 
#'   algorithm. Will be random if not set.
#' 
#' @inheritParams clusterRows
#' 
#' @author Basil Courbayre
#'
#' @references
#' Holmes I, Harris K and Quince C (2012).
#' Dirichlet multinomial mixtures: generative models for microbial metagenomics.
#' \emph{PLoS ONE, vol. 7(2)}, 1-15
#' 
#' @details
#' The DMM algorithm (see Holmes et al. 2012) is commonly used in microbial 
#' ecology along with metagenomic and 16S rRNA count data.
#' 
#' For example, this algorithm could be used on an assay coming from a 
#' \code{SummarizedExperiment}.
#' 
#' To modify an existing DmmParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} 
#' is any argument used in the constructor.
#'
#' The default for \code{k} is \code{1:3}.
#' 
#' The default for \code{type} is \code{"laplace"}.
#' 
#' @return 
#' The \code{DmmParam} constructor will return a \linkS4class{DmmParam} object 
#' with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to 
#' \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} 
#' (the factor, as above) and \code{objects} 
#' (a list containing \code{dmm}, the output of 
#' \code{\link[DirichletMultinomial]{dmn}} for each value contained in 
#' \code{k}; \code{k}, the number of clusters; \code{prob}, the array of 
#' probabilities; and \code{seed} the seed for the dmm).
#'
#' @examples
#' clusterRows(t(iris[, 1:4]), DmmParam())
#' clusterRows(t(iris[, 1:4]), DmmParam(k = 2))
#' 
#' fl <- system.file(package="DirichletMultinomial", "extdata", "Twins.csv")
#'         counts <- t(as.matrix(read.csv(fl, row.names=1)))
#'         clusterRows(counts, DmmParam(k=1:3, type="laplace"))
#'
#' @name DmmParam-class
#' @docType class
#' @aliases 
#' show,DmmParam-method
NULL

#' @export
#' @rdname DmmParam-class
setClass("DmmParam", contains="BlusterParam", slots=c(k="integer_OR_NULL", 
                                                      type="character",
                                                      transposed="logical",
                                                      seed = "integer_OR_NULL"))

#' @export
#' @rdname DmmParam-class
DmmParam <- function(k=1:3, type="laplace", transposed=FALSE, seed=NULL) {
    # Filling in missing values with the defaults.
    if (is.null(seed)) {
        seed=runif(1, 0, .Machine$integer.max)
    }
    
    # Formatting data
    k <- as.integer(k)
    seed <- as.integer(seed)

    new("DmmParam", k=k, type=type, 
        transposed=transposed, seed=seed)
}

#' @export
setMethod("show", "DmmParam", function(object) {
    callNextMethod()
    if (length(object@k) > 1) {
        cat(sprintf("k: %s\n", paste(object@k, collapse = ", ")))
    } else {
        cat(sprintf("k: %s\n", object@k))
    }
    cat(sprintf("type: %s\n", object@type))
    cat(sprintf("transposed: %s\n", object@transposed))
    cat(sprintf("seed: %s\n", object@seed))
})

setValidity2("DmmParam", function(object) {
    msg <- character(0)
    if (!is.integer(object@k) ||
        length(object@k) < 1 ||
        any(object@k <= 0)) {
        msg <- c(msg, "'k' must be a strictly positive integer vector")
    }
    if (length(object@type) > 1) {
        msg <- c(msg, "'type' must be a string")
    } else {
        match.arg(object@type, c("laplace","AIC","BIC"))
    }
    
    if (length(msg) > 0) {
        msg
    } else {
        TRUE
    }
})

#' @export
#' @rdname DmmParam-class
setMethod("clusterRows", c("ANY", "DmmParam"), function(x, 
                                                        BLUSPARAM, 
                                                        full=FALSE) {
    k <- BLUSPARAM[["k"]]
    type <- BLUSPARAM[["type"]]
    transposed <- BLUSPARAM[["transposed"]]
    seed <- BLUSPARAM[["seed"]]

    if (!transposed) {
        x <- t(x)
    }
    
    dmm <- .get_dmm(x, k=k, seed = seed)
    
    # Finding optimal number of clusters if necessary
    if (length(k) > 1) {
        k <- .best_dmm_fit(dmm, type)
    } 
    
    # Get the index corresponding to k in dmm list
    i <- which(
        vapply(dmm, function(x) ncol(DirichletMultinomial::mixture(x)) == k,
               logical(1)))

    prob <- DirichletMultinomial::mixture(dmm[[i]])
    
    # Formatting the result
    colnames(prob) <- 1:k
    clusters <- colnames(prob)[max.col(prob, ties.method = "first")]
    clusters <- factor(clusters)
    names(clusters) <- rownames(x)
    
    if (full) {
        list(clusters=clusters, 
             objects=list(dmm=dmm, k=k, prob=prob, seed=seed))
    } else {
        clusters
    }
})

#' @importFrom stats runif
#' @importFrom BiocParallel bplapply
.get_dmm <- function(x, k = 1, BPPARAM = SerialParam(),
                     seed = runif(1, 0, .Machine$integer.max), ...){
    old <- DelayedArray::getAutoBPPARAM()
    DelayedArray::setAutoBPPARAM(BPPARAM)
    on.exit(DelayedArray::setAutoBPPARAM(old))
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }
    
    dmm <- BiocParallel::bplapply(k, DirichletMultinomial::dmn, count = x,
                                  seed = seed, ...,
                                  BPPARAM = BPPARAM)
    dmm
}

.get_dmm_fit_FUN <- function(type){
    fit_FUN <- switch(type,
                      laplace = DirichletMultinomial::laplace,
                      AIC = DirichletMultinomial::AIC,
                      BIC = DirichletMultinomial::BIC)
    fit_FUN
}

.best_dmm_fit <- function(dmm, type){
    fit_FUN <- .get_dmm_fit_FUN(type)
    fit <- vapply(dmm, fit_FUN, numeric(1))
    which.min(fit)
}