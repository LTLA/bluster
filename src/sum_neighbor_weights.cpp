#include "Rcpp.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix sum_neighbor_weights(int nclusters, Rcpp::List neighbors, Rcpp::IntegerVector clusters, Rcpp::NumericVector weights) {
    const int npts = neighbors.size();
    Rcpp::NumericMatrix output(nclusters, npts);

    for (int i=0; i<npts; ++i) {
        auto COL = output.column(i);
        COL[clusters[i]] = weights[i]; // setting the weight of the point itself.

        Rcpp::IntegerVector NN = neighbors[i];
        for (auto n : NN) {
            auto current = n - 1; // get back to 0-based indexing.
            COL[clusters[current]] += weights[current];
        }
    }

    return output;
}
