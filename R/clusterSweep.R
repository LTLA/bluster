#' Clustering parameter sweeps
#'
#' Perform a sweep across combinations of parameters to obtain different clusterings from the same algorithm.
#'
#' @inheritParams clusterRows
#' @param ... Named vectors or lists specifying the parameters to sweep over.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the sweep should be parallelized.
#' @param args A named list of additional arguments to use with \code{...}.
#' This is provided in case there is a name conflict with the existing arguments in this function signature.
#' 
#' @return
#' A \linkS4class{List} containing:
#' \itemize{
#' \item \code{clusters}, a \linkS4class{DataFrame} with number of rows equal to that of \code{x},
#' where each column corresponds to (and is named after) a specific combination of clustering parameters.
#' \item \code{parameters}, another DataFrame with number of rows equal to the number of columns in the previous \code{clusters} DataFrame.
#' Each row contains the specific parameter combination for each column of \code{clusters}.
#' \item If \code{full=TRUE}, \code{objects} is an additional list of length equal to the number of rows in \code{clusters}.
#' This contains the \code{objects} produced by each run.
#' }
#'
#' @details
#' This function allows users to conveniently test out a range of clustering parameters in a single call.
#' The name of each argument in \code{...} should be a legitimate argument to \code{x[[i]]},
#' and will be used to modify any existing values in \code{BLUSPARAM} to obtain a new set of parameters.
#' (For all other parameters, the existing values in \code{BLUSPARAM} are used.)
#' If multiple arguments are provided, all combinations are tested.
#'
#' We attempt to create a unique name for each column based on its parameter combination.
#' This has the format of \code{<NAME1>.<VALUE1>_<NAME2>.<VALUE2>_...} based on the parameter names and values.
#' Note that any non-atomic values are simply represented by the name of their class;
#' no attempt is made to convert these into a compact string.
#' 
#' If an entry of \code{...} is a \emph{named} list of vectors, we expand those to generate all possible combinations of values.
#' For example, if we passed:
#' \preformatted{    blah.args = list(a = 1:5, b = LETTERS[1:3])}
#' This would be equivalent to manually specifying:
#' \preformatted{    blah.args = list(list(a = 1, b = "A"), list(a = 1, b = "B"), ...)}
#' The auto-expansion mechanism allows us to conveniently test parameter combinations when those parameters are stored inside \code{x} as a list.
#' The algorithm is recursive so any internal named lists containing vectors are similarly expanded.
#' Expansion can be disabled by wrapping vectors in \code{\link{I}}, in which case they are passed verbatim.
#' No expansion is performed for non-vector arguments.
#'
#' @author Aaron Lun
#'
#' @examples
#' out <- clusterSweep(iris[,1:4], KmeansParam(10), 
#'     centers=4:10, algorithm=c("Lloyd", "Hartigan-Wong"))
#' out$clusters[,1:5]
#' out$parameters
#'
#' out <- clusterSweep(iris[,1:4], NNGraphParam(), k=c(5L, 10L, 15L, 20L),
#'     cluster.fun=c("louvain", "walktrap"))
#' out$clusters[,1:5]
#' out$parameters
#' 
#' # Combinations are automatically expanded inside named lists:
#' out <- clusterSweep(iris[,1:4], NNGraphParam(), k=c(5L, 10L, 15L, 20L), 
#'     cluster.args=list(steps=3:4))
#' out$clusters[,1:5]
#' out$parameters
#' 
#' @seealso
#' \code{\link{clusterRows}}, which manages the dispatch to specific methods based on \code{BLUSPARAM}.
#'
#' \linkS4class{BlusterParam}, which determines which algorithm is actually used.
#'
#' @export
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom S4Vectors I DataFrame List
clusterSweep <- function(x, BLUSPARAM, ..., full=FALSE, BPPARAM=SerialParam(), args=list()) {
    args <- c(list(...), args)
    if (is.null(names(args)) || anyDuplicated(names(args))) {
        stop("all arguments in '...' and 'args' must have unique names")
    }

    for (i in names(args)) {
        current <- args[[i]]
        if (is.list(current) && !is.null(names(current))) {
            args[[i]] <- .expand_list(current)
        }
    }

    N <- lengths(args)
    all.choices <- do.call(expand.grid, lapply(N, seq_len))

    out <- bplapply(seq_len(nrow(all.choices)), FUN=.run_combination, BPPARAM=BPPARAM, 
        x=x, combinations=all.choices, parameters=args, full=full, BLUSPARAM=BLUSPARAM)

    param.list <- lapply(seq_along(all.choices), FUN=function(i) args[[i]][all.choices[,i]])
    names(param.list) <- names(args)
    param.df <- do.call(DataFrame, lapply(param.list, I))
    rownames(param.df) <- names(out) <- .create_combination_names(param.list)

    output <- List()
    if (full) {
        output$clusters <- DataFrame(lapply(out, function(x) x$clusters), check.names=FALSE)
        output$parameters <- param.df
        output$objects <- lapply(out, function(x) x$objects)
    } else {
        output$clusters <- DataFrame(out, check.names=FALSE)
        output$parameters <- param.df
    }

    output
}

.run_combination <- function(i, x, BLUSPARAM, combinations, parameters, full) {
    current <- combinations[i,]
    for (j in seq_along(parameters)) {
        curparams <- parameters[[j]]
        BLUSPARAM[[names(parameters)[j]]] <- curparams[[current[[j]]]]
    }
    clusterRows(x, BLUSPARAM=BLUSPARAM, full=full)
}

.expand_list <- function(values) {
    collected <- list(list())
    for (field in names(values)) {
        current <- values[[field]]
        if (is.list(current) && !is.null(names(current))) {
            current <- .expand_list(current)
        } else if (!is.vector(current)) {
            if (is(current, "AsIs")) {
                classes <- class(current)
                new_class <- classes[classes != "AsIs"]
                attr(new_class, "package") <- attr(classes, "package")
                class(current) <- new_class
            }
            current <- list(current)
        }

        # Expanding by these values.
        for (i in seq_along(collected)) {
            added <- vector("list", length(current))
            for (j in seq_along(current)) {
                stub <- collected[[i]]
                stub[[field]] <- current[[j]]
                added[[j]] <- stub
            }
            collected[[i]] <- added
        }
        collected <- unlist(collected, recursive=FALSE)
    }

    collected 
}

.create_combination_names <- function(parameters) {
    for (i in names(parameters)) {
        values <- parameters[[i]]
        if (!is.atomic(values)) {
            values <- rep(class(values)[1], length(values))
        }
        parameters[[i]] <- sprintf("%s.%s", i, values)
    }
    do.call(paste, c(parameters, list(sep="_")))
}
