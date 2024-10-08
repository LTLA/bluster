% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DbscanParam.R
\docType{class}
\name{DbscanParam-class}
\alias{DbscanParam-class}
\alias{show,DbscanParam-method}
\alias{DbscanParam}
\alias{clusterRows,ANY,DbscanParam-method}
\title{Density-based clustering with DBSCAN}
\usage{
DbscanParam(
  eps = NULL,
  min.pts = 5,
  core.prop = 0.5,
  chunk.size = 1000,
  BNPARAM = KmknnParam(),
  num.threads = 1,
  BPPARAM = NULL
)

\S4method{clusterRows}{ANY,DbscanParam}(x, BLUSPARAM, full = FALSE)
}
\arguments{
\item{eps}{Numeric scalar specifying the distance to use to define neighborhoods.
If \code{NULL}, this is determined from \code{min.pts} and \code{core.prop}.}

\item{min.pts}{Integer scalar specifying the minimum number of neighboring observations required for an observation to be a core point.}

\item{core.prop}{Numeric scalar specifying the proportion of observations to treat as core points.
This is only used when \code{eps=NULL}, see Details.}

\item{chunk.size}{Integer scalar specifying the number of points to process per chunk.}

\item{BNPARAM}{A \linkS4class{BiocNeighborParam} object specifying the algorithm to use for the neighbor searches.
This should be able to support both nearest-neighbor and range queries.}

\item{num.threads}{Integer scalar specifying the number of threads to use.}

\item{BPPARAM}{Deprecated and ignored, use \code{num.threads} instead.}

\item{x}{A numeric matrix-like object where rows represent observations and columns represent variables.}

\item{BLUSPARAM}{A \linkS4class{BlusterParam} object specifying the algorithm to use.}

\item{full}{Logical scalar indicating whether additional statistics should be returned.}
}
\value{
The \code{DbscanParam} constructor will return a \linkS4class{DbscanParam} object with the specified parameters.

The \code{clusterRows} method will return a factor of length equal to \code{nrow(x)} containing the cluster assignments.
Note that this may contain \code{NA} values corresponding to noise points.
If \code{full=TRUE}, a list is returned with \code{clusters} (the factor, as above) and \code{objects} 
(a list containing the \code{eps} and \code{min.pts} used in the analysis).
}
\description{
Perform density-based clustering with a fast re-implementation of the DBSCAN algorithm.
}
\details{
DBSCAN operates by identifying core points, i.e., observations with at least \code{min.pts} neighbors within a distance of \code{eps}.
It identifies which core points are neighbors of each other, one \code{chunk.size} at a time, forming components of connected core points.
All non-core points are then connected to the closest core point within \code{eps}.
All groups of points that are connected in this manner are considered to be part of the same cluster.
Any unconnected non-core points are treated as noise and reported as \code{NA}.

As a suitable value of \code{eps} may not be known beforehand, we can automatically determine it from the data.
For all observations, we compute the distance to the \eqn{k}th neighbor where \eqn{k} is defined as \code{round(min.pts * core.prop)}.
We then define \code{eps} as the \code{core.prop} quantile of the distances across all observations.
The default of \code{core.prop=0.5} means that around half of the observations will be treated as core points.

Larger values of \code{eps} will generally result in fewer observations classified as noise, as they are more likely to connect to a core point.
It may also promote agglomeration of existing clusters into larger entities if they are connected by regions of (relatively) low density.
Conversely, larger values of \code{min.pts} will generally increase the number of noise points and may fragment larger clusters into subclusters.

To modify an existing DbscanParam object \code{x},
users can simply call \code{x[[i]]} or \code{x[[i]] <- value} where \code{i} is any argument used in the constructor.
}
\examples{
clusterRows(iris[,1:4], DbscanParam())
clusterRows(iris[,1:4], DbscanParam(core.prop=0.8))

}
\references{
Ester M et al. (1996).
A density-based algorithm for discovering clusters in large spatial databases with noise.
\emph{Proceedings of the Second International Conference on Knowledge Discovery and Data Mining}, 226-231.
}
\author{
Aaron Lun
}
