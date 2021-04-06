#' Map nested clusterings 
#'
#' Map an alternative clustering to a reference clustering,
#' where the latter is expected to be nested within the former.
#'
#' @inheritParams pairwiseRand
#'
#' @return
#' A list containing:
#' \itemize{
#' \item \code{proportions}, a matrix where each row corresponds to one of the \code{alt} clusters and each column corresponds to one of the \code{ref} clusters.
#' Each matrix entry represents the proportion of cells in \code{alt} that are assigned to each cluster in \code{ref}.
#' (That is, the proportions across all \code{ref} clusters should sum to unity for each \code{alt} cluster.)
#' \item \code{alt.mapping}, a \linkS4class{DataFrame} with one row per cluster in \code{alt}.
#' This contains the columns \code{max}, a numeric vector specifying the maximum value of \code{statistic} for that \code{alt} cluster;
#' and \code{which}, a character vector specifying the \code{ref} cluster in which the maximum value occurs.
#' \item \code{ref.score}, a numeric vector of length equal to the number of \code{ref} clusters.
#' This represents the degree of nesting of \code{alt} clusters within each \code{ref} cluster, see Details.
#' }
#' 
#' @details
#' This function identifies mappings between two clusterings on the same set of cells where \code{alt} is potentially nested within \code{ref} (e.g., as it is computed at higher resolution).
#' To do so, we take each \code{alt} cluster and compute the the proportion of its cells that are derived from each \code{ref} cluster.
#' The corresponding \code{ref} cluster is identified as that with the highest proportion, as reported by the \code{which} field in the \code{mapping} DataFrame. 
#'
#' The quality of the mapping is determined by \code{max} in the output \code{mapping} DataFrame.
#' A low value indicates that \code{alt} does not have a clear counterpart in \code{ref}, representing loss of heterogeneity.
#' Note that this is not a symmetrical inference; multiple \code{alt} clusters can map to the same \code{ref} cluster without manifesting as a low \code{max}.
#' This implicitly assumes that an increase in resolution in \code{alt} is not problematic.
#'
#' The \code{ref.score} value for each cluster \code{ref} is formally defined as the probability of randomly picking a cell that belongs to \code{ref},
#' conditional on the event that the chosen cell belongs to the same \code{alt} cluster as a randomly chosen cell from \code{ref}.
#' This probability is equal to unity when \code{ref} is an exact superset of all \code{alt} clusters that contain its cells, corresponding to perfect 1:many nesting.
#' In contrast, if the \code{alt} clusters contain a mix of cells from different \code{ref}, this probability will be low and can be used as a diagnostic for imperfect nesting.
#'
#' @examples
#' m <- matrix(runif(10000), ncol=10)
#' clust1 <- kmeans(m,10)$cluster
#' clust2 <- kmeans(m,20)$cluster
#' nestedClusters(clust1, clust2)
#'
#' # The ref.score is 1 in cases of perfect nesting.
#' nestedClusters(clust1, clust1)$ref.score
#'
#' nest.clust <- paste0(clust1, sample(letters, length(clust1), replace=TRUE))
#' nestedClusters(clust1, nest.clust)$ref.score
#'
#' # In contrast, it is much lower when nesting is bad.
#' nestedClusters(clust1, sample(clust1))$ref.score
#' 
#' @seealso
#' \code{\link{linkClusters}}, to do this in a symmetric manner (i.e., without nesting).
#'
#' \code{\link{pairwiseRand}}, for another way of comparing two sets of clusterings.
#'
#' @export
#' @importFrom S4Vectors DataFrame
nestedClusters <- function(ref, alt) {
    tab <- table(alt, ref)
    stats <- tab/rowSums(tab)
    m <- max.col(stats, ties.method="first")
    list(
         proportions=.untable(stats),
         alt.mapping=DataFrame(row.names=rownames(stats),
                               max=stats[cbind(seq_along(m), m)],
                               which=colnames(stats)[m]),

         # To explain: 'stats' represents the probability of picking a cell
         # that belongs to 'ref', conditional on having chosen a particular
         # 'alt' cluster.  The probability of choosing that 'alt' cluster is
         # equal to the distribution of cells across 'alt' clusters for a given
         # 'ref', which is why we multiply by 'tab'. We then divide by the
         # total number of cells in each 'ref' (i.e., 'colSums(tab)') to obtain
         # an actual probability.
         ref.score=colSums(stats * tab) / colSums(tab)
    )
}
