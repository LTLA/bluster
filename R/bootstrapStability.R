#' Assess cluster stability by bootstrapping
#'
#' Generate bootstrap replicates and recluster on them to determine the stability of clusters with respect to sampling noise.
#'
#' @param x A numeric matrix-like object containing observations in the rows and variables in the columns.
#' @param FUN A function that takes \code{x} as its first argument and returns a vector or factor of cluster identities.
#' @param clusters A vector or factor of cluster identities equivalent to that obtained by calling \code{FUN(x, ...)}.
#' This is provided as an additional argument in the case that the clusters have already been computed,
#' in which case we can save a single round of computation.
#' @param iterations A positive integer scalar specifying the number of bootstrap iterations.
#' @param average String specifying the method to use to average across bootstrap iterations.
#' @param ... Further arguments to pass to \code{FUN} to control the clustering procedure.
#' @param compare A function that accepts the original clustering and the bootstrapped clustering,
#' and returns a numeric vector or matrix containing some measure of similarity between them - see Details.
#' @param mode,adjusted Further arguments to pass to \code{\link{pairwiseRand}} when \code{compare=NULL}.
#'
#' @return 
#' If \code{compare=NULL} and \code{mode="ratio"}, a numeric matrix is returned with upper triangular entries set to the ratio of the adjusted observation pair counts (see \code{?\link{pairwiseRand}}) for each pair of clusters in \code{clusters}.
#' Each ratio is averaged across bootstrap iterations as specified by \code{average}.
#' 
#' If \code{compare=NULL} and \code{mode="index"}, a numeric scalar containing the average ARI between \code{clusters} and the bootstrap replicates across iterations is returned.
#'
#' If \code{compare} is provided, a numeric array of the same type as the output of \code{compare} is returned, containing the average statistic(s) across bootstrap replicates.
#'
#' @details
#' Bootstrapping is conventionally used to evaluate the precision of an estimator by applying it to an \emph{in silico}-generated replicate dataset.
#' We can (ab)use this framework to determine the stability of the clusters given the original dataset.
#' We sample observations with replacement from \code{x}, perform clustering with \code{FUN} and compare the new clusters to \code{clusters}.
#'
#' For comparing clusters, we compute the ratio matrix from \code{\link{pairwiseRand}} and average its values across bootstrap iterations.
#' High on-diagonal values indicate that the corresponding cluster remains coherent in the bootstrap replicates,
#' while high off-diagonal values indicate that the corresponding pair of clusters are still separated in the replicates.
#' If a single value is necessary, we can instead average the adjusted Rand indices across iterations with \code{mode="index"}.
#' 
#' We use the ratio matrix by default as it is more interpretable than a single value like the ARI or the Jaccard index (see the \pkg{fpc} package).
#' It focuses on the relevant differences between clusters, allowing us to determine which aspects of a clustering are stable.
#' For example, A and B may be well separated but A and C may not be, which is difficult to represent in a single stability measure for A.
#' If our main interest lies in the A/B separation, we do not want to be overly pessimistic about the stability of A, even though it might not be well-separated from all other clusters.
#'
#' @section Using another comparison function:
#' We can use a different method for comparing clusterings by setting \code{compare}.
#' This is expected to be a function that takes two arguments - 
#' the original clustering first, and the bootstrapped clustering second - 
#' and returns some kind of numeric scalar, vector or matrix containing 
#' statistics for the similarity or difference between the original and bootstrapped clustering.
#' These statistics are then averaged across all bootstrap iterations.
#' 
#' Any numeric output of \code{compare} is acceptable as long as the dimensions are only dependent on the \emph{levels} of the original clustering - including levels that have no observations, due to resampling! - and thus do not change across bootstrap iterations.
#'
#' @section Statistical note on bootstrap comparisons:
#' Technically speaking, some mental gymnastics are required to compare the original and bootstrap clusters in this manner.
#' After bootstrapping, the sampled observations represent distinct entities from the original dataset (otherwise it would be difficult to treat them as independent replicates) for which the original clusters do not immediately apply.
#' Instead, we assume that we perform label transfer using a nearest-neighbors approach - which, in this case, is the same as using the original label for each observation, as the nearest neighbor of each resampled observation to the original dataset is itself.
#'
#' Needless to say, bootstrapping will only generate replicates that differ by sampling noise.
#' Real replicates will differ due to composition differences, variability in expression across individuals, etc.
#' Thus, any stability inferences from bootstrapping are likely to be overly optimistic.
#' 
#' @author Aaron Lun
#' @examples
#' m <- matrix(runif(10000), ncol=10)
#' 
#' # BLUSPARAM just gets passed to the default FUN=clusterRows:
#' bootstrapStability(m, BLUSPARAM=KmeansParam(4), iterations=10)
#'
#' # Defining your own clustering function:
#' kFUN <- function(x) kmeans(x, 2)$cluster  
#' bootstrapStability(m, FUN=kFUN)
#'
#' # Using an alternative comparison, in this case the Rand index:
#' bootstrapStability(m, FUN=kFUN, compare=pairwiseRand)
#'
#' @seealso
#' \code{\link{clusterRows}}, for the default clustering function.
#'
#' \code{\link{pairwiseRand}}, for the calculation of the ARI.
#' 
#' @export
#' @importFrom stats median
bootstrapStability <- function(x, FUN=clusterRows, clusters=NULL, iterations=20, 
    average=c("median", "mean"), ..., compare=NULL, mode="ratio", adjusted=TRUE)
{
    if (is.null(clusters)) {
        clusters <- FUN(x, ...)
    }
    clusters <- as.factor(clusters)

    if (iterations <= 0L) {
        stop("'iterations' must be a positive integer")
    }

    collated <- vector("list", iterations)
    if (is.null(compare)) {
        compare <- function(...) pairwiseRand(..., mode=mode, adjusted=adjusted)
    }

    for (i in seq_len(iterations)) {
        chosen <- sample(nrow(x), nrow(x), replace=TRUE)
        resampled <- x[chosen,,drop=FALSE]
        reclusters <- FUN(resampled, ...)
        collated[[i]] <- compare(clusters[chosen], reclusters)
    }

    if (length(unique(lapply(collated, dim))) > 1L) { 
        stop("'compare' output should have constant dimension")
    }

    # A robust way of computing the average that handles NAs.
    as.mat <- do.call(cbind, lapply(collated, as.numeric))

    average <- match.arg(average)
    if (average=="mean") {
        averaged <- rowMeans(as.mat, na.rm=TRUE)
    } else {
        averaged <- apply(as.mat, 1, median, na.rm=TRUE)
    }

    dim(averaged) <- dim(collated[[1]])
    dimnames(averaged) <- dimnames(collated[[1]])
    averaged 
}
