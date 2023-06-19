#' DMM clustering
#' 
#' Run the \code{DMM} function on a distance matrix 
#' within \code{\link{clusterRows}}.
#' 
#' @param k An integer indicating the number of clusters, or a list of integers 
#'   defining the possible number of clusters to choose from.
#' @param type A string specifying the fit function to use to find the optimal 
#'   number of clusters. Use in combination with k being a list of integers. 
#'   Must be equal to \code{laplace}, \code{AIC} or \code{BIC}.
#' @param transposed Logical scalar, is x transposed with samples in rows?
#' @param seed An integer specifying the seed to use when doing the DMM 
#'   algorithm. Will be random if not set.
#' 
#' @inheritParams clusterRows
#' 
#' @author Basil Courbayre
#'
#' @references
#' Ian Holmes, Keith Harris, Christopher Quince (2012)
#' Dirichlet multinomial mixtures: generative models for microbial metagenomics.
#' \url{https://pubmed.ncbi.nlm.nih.gov/22319561/}
#' 
#' @details
#' The DMM algorithm (see Holmes et al. 2012) is commonly used in microbial ecology along with metagenomic & 16S rRNA count data.
#' 
#' If we have microbiota data (for example the Twins data set in the DMM package), the suitable algorithm can be DMM.
#' To modify an existing DMMParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' If \code{k=NULL}, the cluster search range will be 1:3.
#' The default for the fit function is \code{laplace}.
#' 
#' @return 
#' The \code{DMMParam} constructor will return a \linkS4class{DMMParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects} 
#' (a list containing \code{dmm}, the dmm object; \code{k}, the number of clusters; 
#' \code{prob}, the array of probabilities; and \code{seed} the seed for the dmm).
#'
#' @examples
#' clusterRows(iris[,1:4], DMMParam())
#' clusterRows(iris[,1:4], DMMParam(k=2))
#' \dontrun{
#' fl <- system.file(package="DirichletMultinomial", "extdata", "Twins.csv")
#'         counts <- t(as.matrix(read.csv(fl, row.names=1)))
#'         clusterRows(counts, DMMParam(k=1:3, type="laplace"))
#' }
#' # Example using a SummarizedExperiment
#' data("GlobalPatterns", package="mia")
#' tse <-GlobalPatterns[1:1000, ]
#' x <- assay(tsebis, "counts")
#' dmm <- clusterRows(x, DMMParam(k=1:3))
#'
#' @name DMMParam-class
#' @docType class
#' @aliases 
#' .defaultScalarArguments,DMMParam-method
#' show,DMMParam-method
NULL

#' @export
#' @rdname DMMParam-class
setClass("DMMParam", contains="BlusterParam", slots=c(k="integer_OR_NULL", 
                                                      type="character",
                                                      transposed="logical",
                                                      seed = "integer_OR_NULL"))

#' @export
#' @rdname DMMParam-class
DMMParam <- function(k=NULL, type=NULL, transposed=FALSE, seed=NULL) {
    # Filling in missing values with the defaults.
    current <- list(k=k, type=type, seed=seed)
    notpresent <- vapply(current, is.null, FALSE)
    if (any(notpresent)) {
        defaults <- list(k=1:3, 
                         type="laplace", 
                         seed=runif(1, 0, .Machine$integer.max))
        current[notpresent] <- defaults[notpresent]
    }
    
    # Formatting data
    if (!is.null(current$k)) {
        current$k <- as.integer(current$k)
    }
    if (!is.null(current$seed)) {
        current$seed <- as.integer(current$seed)
    }

    new("DMMParam", k=current$k, type=current$type, 
        transposed=transposed, seed=current$seed)
}

#' @export
setMethod("show", "DMMParam", function(object) {
    callNextMethod()
    if (length(object@k) > 1) {
        cat(sprintf("k: %s\n", paste(range(object@k), collapse = ":")))
    } else {
        cat(sprintf("k: %s\n", object@k))
    }
    cat(sprintf("type: %s\n", object@type))
    cat(sprintf("transposed: %s\n", object@transposed))
    cat(sprintf("seed: %s\n", object@seed))
})

setValidity2("DMMParam", function(object) {
    msg <- character(0)
    if (!is.integer(object@k) ||
        length(object@k) < 1 ||
        any(object@k <= 0)) {
        msg <- c(msg, "'k' must be a strictly positive integer vector.")
    }
    if (length(object@type) > 1) {
        msg <- c(msg, "'type' must be a string")
    } else if (!(object@type %in% c("laplace","AIC","BIC"))) {
        msg <- c(msg, "'type' must be equal to 'laplace', 'AIC', or 'BIC'")
    }
    
    if (length(msg) > 0) {
        msg
    } else {
        TRUE
    }
})

#' @export
#' @rdname DMMParam-class
#' @importFrom package function
setMethod("clusterRows", c("ANY", "DMMParam"), function(x, 
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
        k <- bestDMMFit(dmm, type)
    } 
    
    # Get the index corresponding to k in dmm list
    i <- which(
        sapply(dmm, function(x) ncol(DirichletMultinomial::mixture(x)) == k))

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
    type <- match.arg(type, c("laplace","AIC","BIC"))
    fit_FUN <- switch(type,
                      laplace = DirichletMultinomial::laplace,
                      AIC = DirichletMultinomial::AIC,
                      BIC = DirichletMultinomial::BIC)
    fit_FUN
}

#' @export
#' @rdname DMMParam-class
bestDMMFit <- function(dmm, type){
    fit_FUN <- .get_dmm_fit_FUN(type)
    fit <- vapply(dmm, fit_FUN, numeric(1))
    which.min(fit)
}