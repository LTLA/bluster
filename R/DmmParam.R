#' Dirichlet multinomial mixture clustering
#' 
#' Apply the Dirichlet multinomial mixture (DMM) algorithm from the \pkg{DirichletMultinomial} package.
#' This is commonly used in microbial ecology and in analyses of metagenomic and 16S rRNA count data.
#'
#' @param x A numeric matrix-like object where rows represent observations and columns represent variables.
#' Values are expected to be counts.
#' @inheritParams clusterRows
#' @param k An integer vector indicating the number of clusters to create with the DMM algorithm.
#' A vector containing two or more values will instruct \code{\link{clusterRows}} to perform clustering on each number,
#' and choose the optimal number of clusters based on \code{type}.
#' @param type A string specifying the method to use to find the optimal number of clusters. 
#' Must be equal to \code{"laplace"}, \code{"AIC"} or \code{"BIC"}.
#' Only used when \code{k} contains multiple values.
#' @param seed Integer scalar specifying the seed to use.
#' If \code{NULL}, a random value is used on each invocation of \code{\link{clusterRows}}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how multiple clusterings should be parallelized.
#' Only relevant if \code{k} contains multiple values.
#' 
#' @author Basil Courbayre
#'
#' @references
#' Holmes I, Harris K and Quince C (2012).
#' Dirichlet multinomial mixtures: generative models for microbial metagenomics.
#' \emph{PLoS ONE}, 7(2), 1-15
#' 
#' @details
#' To modify an existing DmmParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#' 
#' @return 
#' The \code{DmmParam} constructor will return a \linkS4class{DmmParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects};
#' the latter is a list containing:
#' \itemize{
#' \item \code{dmm}, a list containing the output of \code{\link[DirichletMultinomial]{dmn}} for each value of \code{k}.
#' \item \code{best}, an integer scalar specifying the best choice of \code{k} according to the method of \code{type}.
#' \item \code{prob}, a matrix array of probabilities where each row is an observation and each column is a cluster.
#' The number of columns is set to the best number of clusters in \code{best}.
#' \item \code{seed}, an integer scalar specifying the seed used for clustering.
#' }
#'
#' @examples
#' # Mocking up a small example.
#' nfeatures <- 50
#' out1 <- matrix(rpois(20 * nfeatures, lambda = rgamma(nfeatures, 5)), ncol=nfeatures, byrow=TRUE)
#' out2 <- matrix(rpois(20 * nfeatures, lambda = rgamma(nfeatures, 5)), ncol=nfeatures, byrow=TRUE)
#' out <- rbind(out1, out2)
#' clusterRows(out, DmmParam())
#'
#' @name DmmParam-class
#' @docType class
#' @aliases 
#' show,DmmParam-method
NULL

#' @export
#' @rdname DmmParam-class
setClass("DmmParam", contains="BlusterParam", slots=c(
    k = "integer_OR_NULL", 
    type = "character", 
    seed = "integer_OR_NULL",
    BPPARAM = "BiocParallelParam"
))

.allowed.dmm.type <- c("laplace", "AIC", "BIC")

#' @export
#' @rdname DmmParam-class
#' @importFrom BiocParallel SerialParam
DmmParam <- function(k = 1:3, type = "laplace", seed = NULL, BPPARAM = SerialParam()) {
    new("DmmParam", k = as.integer(k), type=match.arg(type, .allowed.dmm.type), seed=seed, BPPARAM = BPPARAM)
}

#' @export
#' @importFrom S4Vectors coolcat
setMethod("show", "DmmParam", function(object) {
    callNextMethod()
    coolcat("k(%i): %s", as.character(object@k))
    cat(sprintf("type: %s\n", object@type))
    cat(sprintf("seed: %s\n", object@seed))
    cat(sprintf("BPPARAM: %s\n", class(object@seed)[1]))
})

setValidity2("DmmParam", function(object) {
    msg <- character(0)

    if (!is.integer(object@k) || length(object@k) < 1 || any(object@k <= 0)) {
        msg <- c(msg, "'k' must be a strictly positive integer vector")
    }

    if (!(object@type %in% .allowed.dmm.type)) {
        msg <- c(msg, "'type' must be one of \"AIC\", \"BIC\" or \"laplace\"")
    }

    if (length(msg) > 0) {
        msg
    } else {
        TRUE
    }
})

#' @export
#' @rdname DmmParam-class
#' @importFrom stats runif
#' @importFrom BiocParallel bplapply bpisup bpstart bpstop
setMethod("clusterRows", c("ANY", "DmmParam"), function(x, BLUSPARAM, full=FALSE) {
    k <- BLUSPARAM[["k"]]
    type <- BLUSPARAM[["type"]]

    # Assigning the seed here to ensure reproducibility
    # when multiple threads are used.
    seed <- BLUSPARAM[["seed"]]
    if (is.null(seed)) {
        seed <- sample(.Machine$integer.max, 1)
    }

    BPPARAM <- BLUSPARAM[["BPPARAM"]]
    if (!(bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam"))) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add = TRUE)
    }

    if (is.null(dimnames(x))) { # otherwise dmn() complains.
        dimnames(x) <- list(NULL, NULL)
    }
    dmm <- bplapply(k, DirichletMultinomial::dmn, count = x, seed = seed, BPPARAM = BPPARAM)

    # Finding optimal number of clusters if necessary
    best <- 1L
    if (length(k) > 1) {
        fit_FUN <- switch(type,
            laplace = DirichletMultinomial::laplace,
            AIC = DirichletMultinomial::AIC,
            BIC = DirichletMultinomial::BIC)
        fit <- vapply(dmm, fit_FUN, numeric(1))
        best <- which.min(fit)
    }

    prob <- DirichletMultinomial::mixture(dmm[[best]])
    colnames(prob) <- seq_len(ncol(prob))
    clusters <- colnames(prob)[max.col(prob, ties.method = "first")]
    clusters <- factor(clusters)
    names(clusters) <- rownames(x)

    if (full) {
        list(clusters=clusters, objects=list(dmm=dmm, best=k[[best]], prob=prob, seed=seed))
    } else {
        clusters
    }
})
