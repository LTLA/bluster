#' Approximate silhouette width
#'
#' Given a clustering, quickly compute an approximate silhouette width for each observation.
#'
#' @inheritParams neighborPurity 
#' 
#' @details
#' The silhouette width is a general-purpose method for evaluating the separation between clusters 
#' but requires calculating the average distance between pairs of observations within or between clusters.
#' This function instead approximates the average distance with the root-mean-squared-distance, which can be computed very efficiently for large datasets.
#' The approximated averages are then used to compute the silhouette width using the usual definition.
#'
#' % In case you're wondering, here it is. Let O be the current observation.
#' % - The squared distance to the centroid is (O - Y)^2 where Y = mean(X) for the vector of observations X.
#' % - The mean of squares within X is mean[(X - Y)^2].
#' % - Their sum is:
#' %   (O - Y)^2 + mean[(X - Y)^2] = O^2 - 2OY + Y^2 + mean(X^2) - 2 mean(X) Y + Y^2
#' %                               = O^2 - 2OY + mean(X^2)
#' %                               = O^2 - 2OE(X) + mean(X^2)
#' %                               = mean[(O - X)^2]
#' %   After which we can just take the square root.
#' 
#' @return
#' A \linkS4class{DataFrame} with one row per observation in \code{x} and the columns:
#' \itemize{
#' \item \code{cluster}, the assigned cluster for each observation in \code{x}.
#' \item \code{other}, the closest cluster other than the one to which the current observation is assigned.
#' \item \code{width}, a numeric field containing the approximate silhouette width of the current cell.
#' }
#' Row names are defined as the row names of \code{x}.
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
        clust.var[i] <- sum(colMeans(sweep(xcurrent, 2, centroid)^2))
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
        cluster=clusters,
        other=uclust[other.clust],
        width=(other.dist - self.dist)/pmax(other.dist, self.dist),
        row.names=rownames(x)
    )
}
