#' Cluster rows of a matrix
#'
#' Cluster rows of a matrix-like object with a variety of algorithms.
#'
#' @param x A numeric matrix-like object where rows represent observations and columns represent variables.
#' @param BLUSPARAM A \linkS4class{BlusterParam} object specifying the algorithm to use.
#' @param full Logical scalar indicating whether the full clustering statistics should be returned for each method.
#'
#' @return
#' By default, a factor of length equal to \code{nrow(x)} containing cluster assignments for each row of \code{x}.
#'
#' If \code{full=TRUE}, a list is returned containing \code{clusters}, a factor as described above;
#' and \code{objects}, an arbitrary object containing algorithm-specific statistics or intermediate objects.
#'
#' @details
#' This generic allows users to write agile code that can use a variety of clustering algorithms.
#' By simply changing \code{BLUSPARAM}, we can tune the clustering procedure in analysis workflows and package functions.
#'
#' @seealso
#' \linkS4class{HclustParam}, \linkS4class{KmeansParam} and \linkS4class{NNGraphParam} 
#' for some examples of values for \code{BLUSPARAM}.
#' 
#' @author Aaron Lun
#'
#' @examples
#' m <- matrix(runif(10000), ncol=10)
#'
#' clusterRows(m, KmeansParam(10L))
#' clusterRows(m, HclustParam())
#' clusterRows(m, NNGraphParam())
#' @export
setGeneric("clusterRows", function(x, BLUSPARAM, full=FALSE) standardGeneric("clusterRows"))

#' @export
setGeneric("centers", function(x, n=NULL) standardGeneric("centers"))

#' @export
setGeneric("centers<-", function(x, value) standardGeneric("centers<-"))

#' Define the default arguments
#'
#' Provide a consistent mechanism to handle specification of default arguments to the underlying clustering functions.
#'
#' @param x,object A \linkS4class{BlusterParam} object.
#'
#' @return
#' For \code{.defaultScalarArguments}, a named character vector is returned.
#' Each entry corresponds to an argument to the clustering function - the name is the argument name, and the value is the argument type.
#'
#' For \code{.extractScalarArguments}, a named list of non-default scalar arguments is returned.
#' Any arguments set to their default values are omitted from the list.
#'
#' For \code{.showScalarArguments}, the values of the arguments are printed to screen.
#' Default values are marked with \code{[default]}.
#'
#' @details
#' The idea is to simplify the derivation of new \linkS4class{BlusterParam} objects,
#' by allowing developers to indicate that the underlying function default should be used for particular arguments.
#' This avoids duplication of the default arguments in the object constructor;
#' instead, default arguments can be indicated as such by setting them to \code{NULL},
#' in which case they will not be explicitly passed to the underlying clustering function.
#'
#' @author Aaron Lun
#' @examples
#' .defaultScalarArguments(PamParam(10))
#' .extractScalarArguments(PamParam(10))
#' .extractScalarArguments(PamParam(10, variant="faster"))
#'
#' @export
#' @aliases
#' .defaultScalarArguments,BlusterParam-method
#'
#' @rdname defaultArguments
setGeneric(".defaultScalarArguments", function(x) standardGeneric(".defaultScalarArguments"))
