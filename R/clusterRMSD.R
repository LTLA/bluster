#' Compute the RMSD per cluster
#'
#' Compute the root mean-squared deviation (RMSD) for each cluster.
#'
#' @param x Numeric matrix containing observations in rows and variables in columns.
#' @param clusters Vector containing the assigned cluster for each observation.
#' @param sum Logical scalar indicating whether to compute the sum of squares.
#'
#' @return 
#' Numeric vector of RMSD values per cluster.
#' If \code{sum=TRUE}, a numeric vector of the sum of squares per cluster is returned instead.
#'
#' @details
#' The RMSD for each cluster is a measure of its dispersion;
#' clusters with large internal heterogeneity will have high RMSDs and are good candidates for further subclustering.
#'
#' @examples
#' x <- matrix(rnorm(10000), ncol=10)
#' kout <- kmeans(x, 5)
#' clusterRMSD(x, kout$cluster)
#'
#' @author Aaron Lun
#' @export
#' @importFrom stats var
clusterRMSD <- function(x, clusters, sum=FALSE) {
    by.cluster <- split(seq_along(clusters), clusters)
    for (i in names(by.cluster)) {
        chosen <- x[by.cluster[[i]],,drop=FALSE]
        by.cluster[[i]] <- sum(apply(chosen, 2, var))
        if (sum) {
            by.cluster[[i]] <- by.cluster[[i]] * (nrow(chosen) - 1)            
        } else {
            by.cluster[[i]] <- sqrt(by.cluster[[i]])
        }
    }

    if (length(by.cluster)) {
        unlist(by.cluster)    
    } else {
        numeric(0)
    }
}
