#' Cluster rows of a matrix
#'
#' Cluster rows of a matrix-like object with a variety of algorithms.
#'
#' @param x A numeric matrix-like object where rows represent observations and columns represent variables.
#' @param BLUSPARAM A \linkS4class{BlusterParam} object specifying the algorithm to use.
#'
#' @return
#' A factor of length equal to \code{nrow(x)} containing cluster assignments for each row of \code{x}.
#'
#' @details
#' This generic allows users to write agile code that can use a variety of clustering algorithms.
#' By simply changing \code{BLUSPARAM}, we can tune the clustering procedure in analysis workflows and package functions.
#'
#' @seealso
#' \linkS4class{HclustParam}, \linkS4class{KmeansParam} and \linkS4class{SnnGraphParam} 
#' for some examples of values for \code{BLUSPARAM}.
#' 
#' @author Aaron Lun
#'
#' @examples
#' m <- matrix(runif(10000), ncol=10)
#'
#' clusterRows(m)
#' @export
setGeneric("clusterRows", function(x, BLUSPARAM) standardGeneric("clusterRows"))
