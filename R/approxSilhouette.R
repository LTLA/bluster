#' Approximate silhouette width
#'
#' Given a clustering, compute a fast approximate silhouette width for each cell.
#'
#' @inheritParams neighborPurity 
#' 
#' @details
#' The silhouette width is a general-purpose method for evaluating the separation between clusters 
#' but requires calculating the average distance between pairs of observations within or between clusters.
#' This function instead approximates the average distances for faster computation in large datasets.
#'
#' For a given observation, let \eqn{\tilde D} be the approximate average distance to all cells in cluster \eqn{X}.
#' This is defined as the square root of the sum of:
#' \itemize{
#' \item The squared distance from the current observation to the centroid of cluster \eqn{X}.
#' This is most accurate when the observation is distant to \eqn{X} relative to the latter's variation.
#' \item The summed variance of all variables across observations in cluster \eqn{X}.
#' This is most accurate when the observation lies close to the close to the centroid of \eqn{X}.
#' }
#' This is also equivalent to the root-square-mean distance from the current observation to all cells in \eqn{X}.
#'
#' % In case you're wondering, here it is. Let O be the current observation.
#' % - The squared distance to the centroid is (O - Y)^2 where Y = E(X).
#' % - The variance within X is E[(X - Y)^2].
#' % - Their sum is:
#' %   (O - Y)^2 + E[(X - Y)^2] = O^2 - 2OY + Y^2 + E(X^2) - 2E(X)Y + Y^2
#' %                            = O^2 - 2OY + E(X^2)
#' %                            = O^2 - 2OE(X) + E(X^2)
#' %                            = E[(O - X)^2]
#'
#' The approximate silhouette width for each cell can then be calculated with the relevant two values of \eqn{\tilde D},
#' computed by setting \eqn{X} to the cluster of the current cell or the closest other cluster.
#' 
#' @return
#' A \linkS4class{DataFrame} with one row per cell in \code{x} and the columns:
#' \itemize{
#' \item \code{width}, a numeric field containing the approximate silhouette width of the current cell.
#' \item \code{other}, the closest cluster other than the one to which the current cell is assigned.
#' }
#' Row names are defined as the column names of \code{x}.
#' 
#' @author Aaron Lun
#' 
#' @examples
#' m <- matrix(rnorm(10000), ncol=10)
#' clusters <- clusterRows(m, BLUSPARAM=KmeansParam(5))
#' out <- approxSilhouette(m, clusters)
#' boxplot(split(out$width, clusters))
#'
#' # Mocking up a stronger example:
#' centers <- matrix(rnorm(30), nrow=3)
#' clusters <- sample(1:3, 1000, replace=TRUE)
#'
#' y <- centers[clusters,]
#' y <- y + rnorm(length(y), sd=0.1)
#' 
#' out2 <- approxSilhouette(y, clusters)
#' boxplot(split(out2$width, clusters))
#'
#' @seealso
#' \code{silhouette} from the \pkg{cluster} package, for the exact calculation.
#'
#' \code{\link{neighborPurity}}, for another method of evaluating cluster separation.
#'
#' @export
#' @importFrom S4Vectors DataFrame
#' @importFrom DelayedMatrixStats colVars
approxSilhouette <- function(x, clusters) {
    x <- as.matrix(x)
    uclust <- sort(unique(clusters))
    averaged <- list(length(uclust))
    clust.var <- numeric(length(uclust))

    for (i in seq_along(uclust)) {
        current <- uclust[i]==clusters
        xcurrent <- x[current,,drop=FALSE]
        centroid <- colMeans(xcurrent)
        averaged[[i]] <- centroid
        clust.var[i] <- sum(colVars(xcurrent, center=centroid))
    }

    self.dist <- other.dist <- rep(Inf, nrow(x))
    other.clust <- integer(nrow(x))
    tx <- t(x)

    for (i in seq_along(uclust)) {
        D <- sqrt(colSums((tx - averaged[[i]])^2) + clust.var[i])

        is.self <- uclust[i]==clusters
        self.dist[is.self] <- D[is.self]

        is.other <- !is.self
        other.D <- D[is.other]
        better <- other.D < other.dist[is.other]
        other.dist[is.other][better] <- other.D[better]
        other.clust[is.other][better] <- i
    }

    DataFrame(
        width=(other.dist - self.dist)/pmax(other.dist, self.dist),
        other=uclust[other.clust],
        row.names=rownames(x)
    )
}
