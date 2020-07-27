#' Two step clustering with vector quantization
#'
#' For large datasets, we can perform vector quantization (e.g., with k-means clustering) to create centroids.
#' These centroids are then subjected to a slower clustering technique such as graph-based community detection.
#' The label for each cell is set to the label of the centroid to which it was assigned.
#'
#' @param first A \linkS4class{BlusterParam} object specifying a fast vector quantization technique.
#' @param second A \linkS4class{BlusterParam} object specifying the second clustering technique on the centroids.
#' @inheritParams clusterRows
#' @param BLUSPARAM A \linkS4class{KmeansParam} object.
#' @param full Logical scalar indicating whether the clustering statistics from both steps should be returned.
#' 
#' @details
#' Here, the idea is to use a fast clustering algorithm to perform vector quantization and reduce the size of the dataset,
#' followed by a slower algorithm that aggregates the centroids for easier interpretation.
#' The exact choice of the number of clusters is less relevant to the first clustering step
#' as long as not too many centroids are generated but the clusters are still sufficiently granular.
#' The second step can take more care (and computational time) summarizing the centroids into meaningful \dQuote{meta-clusters}.
#' 
#' The default choice is to use k-means for the first step, with number of clusters set to the root of the number of observations;
#' and graph-based clustering for the second step, which automatically detects a suitable number of clusters.
#' K-means also eliminates density differences in the data that can introduce variable resolution from graph-based methods.
#'
#' To modify an existing TwoStepParam object \code{x},
#' users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
#'
#' @return 
#' The \code{TwoStepParam} constructor will return a \linkS4class{TwoStepParam} object with the specified parameters.
#'
#' The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
#' If \code{full=TRUE}, a list is returned with a \code{clusters} factor and an \code{objects} list containing:
#' \itemize{
#' \item \code{first}, a list of objects from the first clustering step.
#' This is equal to the \code{objects} list in the output of \code{\link{clusterRows}} with the \code{first} BlusterParam.
#' \item \code{centroids}, a numeric matrix of centroids generated from the first clustering step.
#' \item \code{second}, a list of objects from the second clustering step on the centroids.
#' This is equal to the \code{objects} list in the output of \code{\link{clusterRows}} with the \code{second} BlusterParam.
#' }
#'
#' @author Aaron Lun
#' @examples
#' m <- matrix(runif(100000), ncol=10)
#' stuff <- clusterRows(m, TwoStepParam())
#' table(stuff)
#'
#' @name TwoStepParam-class
#' @aliases
#' show,TwoStepParam-method
NULL

#' @export
#' @rdname TwoStepParam-class
setClass("TwoStepParam", contains="BlusterParam", slots=c(first="BlusterParam", second="BlusterParam"))

#' @export
#' @rdname TwoStepParam-class
TwoStepParam <- function(first=KmeansParam(centers=sqrt), second=NNGraphParam()) {
    new("TwoStepParam", first=first, second=second)
}

#' @export
#' @importFrom utils capture.output
setMethod("show", "TwoStepParam", function(object) {
    callNextMethod()
    cat("first:\n")
    fout <- capture.output(show(object@first))
    cat(paste0("  ", paste(fout, collapse="\n  "), "\n"))
    cat("second:\n")
    sout <- capture.output(show(object@second))
    cat(paste0("  ", paste(sout, collapse="\n  "), "\n"))
})

#' @export
#' @rdname TwoStepParam-class
setMethod("clusterRows", c("ANY", "TwoStepParam"), function(x, BLUSPARAM, full=FALSE) {
    first <- clusterRows(x, BLUSPARAM@first, full=TRUE)

    fclust <- first$clusters
    centroids <- rowsum(x, fclust)
    rn <- rownames(centroids)
    centroids <- centroids/as.integer(table(fclust)[rn])

    second <- clusterRows(centroids, BLUSPARAM@second, full=TRUE)
    clusters <- second$clusters[match(levels(fclust), rn)[fclust]]
    clusters <- factor(clusters)

    if (!full) {
        clusters
    } else {
        list(clusters=clusters, objects=list(centroids=centroids, first=first$objects, second=second$objects))
    }
})

