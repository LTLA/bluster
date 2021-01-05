#' The BlusterParam class
#'
#' The BlusterParam class is a virtual base class controlling S4 dispatch in \code{\link{clusterRows}} and friends.
#' Concrete subclasses specify the choice of clustering algorithm,
#' while the slots of an instance of such a subclass represent the parameters for that algorithm.
#'
#' @section Available methods:
#' In the following code snippets, \code{x} is a \linkS4class{BlusterParam} object or one of its subclasses.
#' \itemize{
#' \item \code{x[[i]]} will return the value of the parameter \code{i}.
#' Refer to the documentation for each concrete subclass for more details on the available parameters.
#' \item \code{x[[i]] <- value} will set the value of the parameter \code{i} to \code{value}.
#' \item \code{show(x)} will print some information about the class instance.
#' }
#' 
#' @author Aaron Lun
#' @seealso
#' \linkS4class{HclustParam}, \linkS4class{KmeansParam} and \linkS4class{NNGraphParam} 
#' for some examples of concrete subclasses.
#'
#' @aliases
#' [[,BlusterParam-method
#' [[<-,BlusterParam-method
#' show,BlusterParam-method
#' @docType class
#' @export
setClass("BlusterParam", contains="VIRTUAL")

setClassUnion("function_OR_NULL", c("function", "NULL"))

setClassUnion("numeric_OR_NULL", c("numeric", "NULL"))

setClassUnion("integer_OR_NULL", c("integer", "NULL"))

setClassUnion("integer_OR_function", c("function", "integer"))

setClassUnion("character_OR_function", c("function", "character"))

#' The FixedNumberParam class
#'
#' The FixedNumberParam is a virtual subclass of the \linkS4class{BlusterParam} class.
#' It causes \code{\link{clusterRows}} to dispatch to clustering algorithms that rely on a pre-specified number of clusters, e.g., \linkS4class{KmeansParam}.
#' 
#' @section Available methods:
#' \code{centers(x, n=NULL)} will return the specified number of centers in a FixedNumberParam \code{x}.
#' This can be an positive integer, or a function that accepts the number of observations and returns a positive number.
#' If a function and \code{n} is supplied, the function is called on \code{n} and the result is rounded to obtain an integer.
#'
#' \code{centers(x) <- value} will replace the specified number of centers in \code{x} with an integer scalar or function \code{value}.
#' The function should accept a single argument and return a positive integer.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \linkS4class{KmeansParam}, for the archetypal example of a concrete subclass.
#'
#' @aliases
#' centers
#' centers,FixedNumberParam-method
#' centers<-
#' centers<-,FixedNumberParam-method
#' show,FixedNumberParam-method
#' 
#' @docType class
#' @export
setClass("FixedNumberParam", contains=c("BlusterParam", "VIRTUAL"), slots=c(centers="integer_OR_function"))
