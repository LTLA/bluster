#' Build a nearest-neighbor graph
#'
#' Build a shared or k-nearest-neighbors graph of observations for downstream community detection.
#'
#' @param x A matrix-like object containing expression values for each observation (row) and dimension (column).
#' @param k An integer scalar specifying the number of nearest neighbors to consider during graph construction.
#' @param type A string specifying the type of weighting scheme to use for shared neighbors.
#' @param directed A logical scalar indicating whether the output of \code{buildKNNGraph} should be a directed graph.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' @param BPPARAM Deprecated, use \code{num.threads} instead.
#' @param num.threads Integer scalar specifying the number of threads to use. 
#' @param indices An integer matrix where each row corresponds to an observation
#' and contains the indices of the \code{k} nearest neighbors (by increasing distance and excluding self) from that observation.
#' 
#' @details
#' The \code{makeSNNGraph} function builds a shared nearest-neighbour graph using observations as nodes.
#' For each observation, its \code{k} nearest neighbours are identified using the \code{\link{findKNN}} function,
#' based on distances between their expression profiles (Euclidean by default).
#' An edge is drawn between all pairs of observations that share at least one neighbour,
#' weighted by the characteristics of the shared nearest neighbors - see \dQuote{Weighting Schemes} below.
#' 
#' The aim is to use the SNN graph to perform clustering of observations via community detection algorithms in the \pkg{igraph} package.
#' This is faster and more memory efficient than hierarchical clustering for large numbers of observations.
#' In particular, it avoids the need to construct a distance matrix for all pairs of observations.
#' Only the identities of nearest neighbours are required, which can be obtained quickly with methods in the \pkg{BiocNeighbors} package.
#' 
#' The choice of \code{k} controls the connectivity of the graph and the resolution of community detection algorithms.
#' Smaller values of \code{k} will generally yield smaller, finer clusters, while increasing \code{k} will increase the connectivity of the graph and make it more difficult to resolve different communities.
#' The value of \code{k} can be roughly interpreted as the anticipated size of the smallest subpopulation.
#' If a subpopulation in the data has fewer than \code{k+1} observations, \code{buildSNNGraph} and \code{buildKNNGraph} will forcibly construct edges between observations in that subpopulation and observations in other subpopulations. 
#' This increases the risk that the subpopulation will not form its own cluster as it is more interconnected with the rest of the observations in the dataset.
#' 
#' The \code{makeKNNGraph} method builds a simpler k-nearest neighbour graph.
#' Observations are again nodes, and edges are drawn between each observation and its k-nearest neighbours.
#' No weighting of the edges is performed.
#' In theory, these graphs are directed as nearest neighour relationships may not be reciprocal.
#' However, by default, \code{directed=FALSE} such that an undirected graph is returned.
#'
#' The \code{neighborsToSNNGraph} and \code{neighborsToKNNGraph} functions operate directly on a matrix of nearest neighbor indices,
#' obtained using functions like \code{\link{findKNN}}.
#' This may be useful for constructing a graph from precomputed nearest-neighbor search results.
#' Note that the user is responsible for ensuring that the indices are valid, i.e., \code{range(indices)} is positive and no greater than \code{max(indices)}.
#' 
#' @section Shared neighbor weighting schemes:
#' If \code{type="rank"}, the weighting scheme defined by Xu and Su (2015) is used.
#' The weight between two nodes is \eqn{k - r/2} where \eqn{r} is the smallest sum of ranks for any shared neighboring node.
#' For example, if one node was the closest neighbor of each of two nodes, the weight between the two latter nodes would be \eqn{k - 1}.
#' For the purposes of this ranking, each node has a rank of zero in its own nearest-neighbor set. 
#' More shared neighbors, or shared neighbors that are close to both observations, will generally yield larger weights.
#'
#' If \code{type="number"}, the weight between two nodes is simply the number of shared nearest neighbors between them.
#' The weight can range from zero to \eqn{k + 1}, as the node itself is included in its own nearest-neighbor set.
#' This is a simpler scheme that is also slightly faster but does not account for the ranking of neighbors within each set.
#'
#' If \code{type="jaccard"}, the weight between two nodes is the Jaccard similarity between the two sets of neighbors for those nodes.
#' This weight can range from zero to 1, and is a monotonic transformation of the weight used by \code{type="number"}.
#' It is provided for consistency with other clustering algorithms such as those in \pkg{seurat}.
#' 
#' Tehcnically, edges with zero weights are assigned a nominal small positive weight of the order of \code{1e-6}.
#' This is done only to satisfy the requirements for positive weights in many \pkg{igraph} clustering algorithms.
#' We do not just remove these edges as that might lead to the situation where some observations have no edges at all and thus form single-observation clusters.
#'
#' Note that the behavior of \code{k} for \code{type="rank"} is slightly different from that used in the original SNN-Cliq implementation by Xu and Su.
#' The original implementation considers each observation to be its first nearest neighbor that contributes to \code{k}.
#' Here, the \code{k} nearest neighbours refers to the number of \emph{other} observations.
#' 
#' @return
#' A \link{graph} where nodes are cells and edges represent connections between nearest neighbors.
#' For \code{buildSNNGraph}, these edges are weighted by the number of shared nearest neighbors.
#' For \code{buildKNNGraph}, edges are not weighted but may be directed if \code{directed=TRUE}.
#' 
#' @author 
#' Aaron Lun, with KNN code contributed by Jonathan Griffiths.
#' 
#' @seealso
#' See \code{\link{make_graph}} for details on the graph output object.
#' 
#' See \code{\link{cluster_walktrap}}, \code{\link{cluster_louvain}} and related functions in \pkg{igraph} for clustering based on the produced graph.
#' 
#' Also see \code{\link{findKNN}} for specifics of the nearest-neighbor search.
#' 
#' @references
#' Xu C and Su Z (2015).
#' Identification of cell types from single-cell transcriptomes using a novel clustering method.
#' \emph{Bioinformatics} 31:1974-80
#' 
#' @examples
#' m <- matrix(rnorm(10000), ncol=10)
#'
#' g <- makeSNNGraph(m)
#' clusters <- igraph::cluster_fast_greedy(g)$membership
#' table(clusters)
#'
#' # Any clustering method from igraph can be used:
#' clusters <- igraph::cluster_walktrap(g)$membership
#' table(clusters)
#'
#' # Smaller 'k' usually yields finer clusters:
#' g <- makeSNNGraph(m, k=5)
#' clusters <- igraph::cluster_walktrap(g)$membership
#' table(clusters)
#'
#' @name makeSNNGraph
NULL

#' @export
#' @rdname makeSNNGraph
#' @importFrom BiocNeighbors KmknnParam findKNN
#' @importFrom BiocParallel SerialParam 
makeSNNGraph <- function(x, k=10, type=c("rank", "number", "jaccard"), BNPARAM=KmknnParam(), num.threads=1, BPPARAM=SerialParam()) { 
    x <- as.matrix(x)
    nn.out <- findKNN(x, k=k, BNPARAM=BNPARAM, num.threads=num.threads, BPPARAM=BPPARAM, get.distance=FALSE)
    neighborsToSNNGraph(nn.out$index, type=match.arg(type), num.threads=max(num.threads, BiocParallel::bpnworkers(BPPARAM)))
}

#' @export
#' @rdname makeSNNGraph
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam findKNN
makeKNNGraph <- function(x, k=10, directed=FALSE, BNPARAM=KmknnParam(), num.threads=1, BPPARAM=SerialParam()) { 
    x <- as.matrix(x)
    nn.out <- findKNN(x, k=k, BNPARAM=BNPARAM, num.threads=num.threads, BPPARAM=BPPARAM, get.distance=FALSE)
    neighborsToKNNGraph(nn.out$index, directed=directed)
}

#' @export
#' @rdname makeSNNGraph
#' @importFrom igraph make_graph E "E<-"
neighborsToSNNGraph <- function(indices, type=c("rank", "number", "jaccard"), num.threads=1) {
    type <- match.arg(type)
    g.out <- build_snn_graph(t(indices), type, num_threads=num.threads)
    edges <- g.out[[1]] 
    weights <- g.out[[2]]
    g <- make_graph(edges, directed=FALSE)
    E(g)$weight <- weights
    g
}

#' @export
#' @rdname makeSNNGraph
#' @importFrom igraph make_graph simplify
neighborsToKNNGraph <- function(indices, directed=FALSE) {
    start <- as.vector(row(indices))
    end <- as.vector(indices)
    interleaved <- as.vector(rbind(start, end))
    
    if (directed) { 
        g <- make_graph(interleaved, directed=TRUE)
    } else {
        g <- make_graph(interleaved, directed=FALSE)
        g <- simplify(g, edge.attr.comb = "first")
    }
    g
}
