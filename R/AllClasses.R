#' The BlusterParam class
#'
#' The BlusterParam class is a virtual base class controlling S4 dispatch in \code{\link{clusterRows}} and friends.
#' Concrete subclasses specify the choice of clustering algorithm,
#' while the slots of an instance of such a subclass represent the parameters for that algorithm.
#'
#' @author Aaron Lun
#' @seealso
#' \linkS4class{HclustParam}, \linkS4class{KmeansParam} and \linkS4class{SnnGraphParam} 
#' for some examples of concrete subclasses.
#'
#' @docType class
#' @export
setClass("BlusterParam", contains="VIRTUAL")
