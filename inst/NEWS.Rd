\name{blusterNews}
\title{bluster News}
\encoding{UTF-8}

\section{Version 1.2.0}{\itemize{
\item
Previously zero-weight edges are now assigned a nominal positive weight in \code{makeSNNGraph()}.

\item 
Added \code{MbkmeansParam()} to wrap mini-batch k-means from \pkg{mbkmeans}.

\item 
Added \code{SOMParam()} to wrap self-organizing map implementation from \pkg{kohonen}.

\item 
Added \code{AffinityParam()} to wrap the affinity propagation code from \pkg{apcluster}.

\item
Added \code{DbscanParam()} to provide a custom DBSCAN implementation with automatic \code{eps} choice.

\item 
Added \code{clusterSweep()} to easily perform parameter sweeps via \code{clusterRows()}.

\item 
Added \code{linkClusters()} to find relationships between clusters in different clusterings.

\item 
Added \code{compareClusterings()} to compute similarities between multiple clusterings.

\item
Added \code{nestedClusters()} to quantify the degree of nesting across two clusterings.
}}

\section{Version 1.0.0}{\itemize{
\item
New package bluster for modular clustering in Bioconductor packages.
}}

