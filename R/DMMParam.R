#' DMM clustering
#' 


#' @export
setClass("DMMParam", contains="BlusterParam", slots=c(variable="factor", 
                                                      k="integer_OR_NULL", 
                                                      type="character",
                                                      transposed="logical",
                                                      seed = "integer_OR_NULL")) # Put type of type character_OR_NULL if possible

#' @export
#' @rdname DMMParam-class
DMMParam <- function(variable, k=NULL, type=NULL, transposed=FALSE, seed=NULL) {
    if (!is.null(k)) {
        k <- as.integer(k)
    }
    if (!is.null(seed)) {
        seed <- as.integer(seed)
    }
    # Filling in missing values with the defaults.
    current <- list(k=k, type=type, seed=seed)
    notpresent <- vapply(current, is.null, FALSE)
    if (any(notpresent)) {
        defaults <- .get_dmm_defaults()
        current[notpresent] <- defaults[notpresent]
    }
    if (!is.factor(variable) && is.character(variable)){
        variable <- factor(variable, unique(variable))
    }
    new("DMMParam", variable=variable, k=current$k,
        type=current$type, transposed=transposed, seed=current$seed)
}

.get_dmm_defaults <- function() {
    out <-list(k=as.integer(1:2), 
               type="laplace", 
               seed=as.integer(runif(1, 0, .Machine$integer.max)))
    out
}

#' @export
setMethod("show", "DMMParam", function(object) {
    callNextMethod()
    cat(sprintf("variable: %s\n", object@variable))
    cat(sprintf("k: %s\n", object@k))
    cat(sprintf("type: %s\n", object@type))
    cat(sprintf("seed: %s\n", object@seed))
})

setValidity2("DMMParam", function(object) {
    msg <- character(0)
    
    if (!is.factor(object@variable)) {
        msg <- c(msg, "'variable' must be a factor or a character value.")
    }
    if (!is.integer(object@k) ||
        length(object@k) < 1) {
        msg <- c(msg, "'variable' must be an integer.")
    }
    if (length(object@type) > 1) {
        msg <- c(msg("'type' must be a string"))
    }
    
    if (length(msg)) {
        msg
    }
    TRUE
})

#' @export
#' @rdname DMMParam-class
#' @importFrom package function
setMethod("clusterRows", c("ANY", "DMMParam"), function(x, 
                                                        BLUSPARAM, 
                                                        full=FALSE) {
    seed <- BLUSPARAM[["seed"]]
    transposed <- BLUSPARAM[["transposed"]]
    k <- BLUSPARAM[["k"]]
    type <- BLUSPARAM[["type"]]
    variable <- BLUSPARAM[["variable"]]
    if (!transposed) {
        x <- t(x)
    }
    
    dmm <- .get_dmm(x, k=k, seed = seed)
    if (length(k) > 1) {
        fit_FUN <- .get_dmm_fit_FUN(type)
        k <- .get_best_nb_clusters(dmm, fit_FUN)
    }
    dmm_group <- .calculate_dmm_group(x, 
                                      variable = variable,
                                      k = k,
                                      seed = .Machine$integer.max)
    
    # Get the index corresponding to k in dmm list
    i <- which(sapply(dmm, 
                      function(x, k) ncol(DirichletMultinomial::mixture(x)) == k, 
                      k=k))[1]
    prob <- DirichletMultinomial::mixture(dmm[[i]])
    show(prob)
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

#' @importFrom DirichletMultinomial dmn
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

.get_best_nb_clusters <- function(dmm, fit_FUN){
    fit <- vapply(dmm, fit_FUN, numeric(1))
    which.min(fit)
}

#' @importFrom DirichletMultinomial dmngroup
#' @importFrom stats runif
.calculate_dmm_group <- function(x, variable, k = 1,
                                seed = runif(1, 0, .Machine$integer.max), 
                                ...) {
    variable <- droplevels(variable)
    DirichletMultinomial::dmngroup(x, variable, k = k, seed = seed, ...)
}

