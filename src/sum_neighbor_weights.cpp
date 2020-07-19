#include "Rcpp.h"

// [[Rcpp::export(rng=false)]]
Rcpp::List sum_neighbor_weights(int nclusters, Rcpp::List neighbors, Rcpp::IntegerVector clusters, Rcpp::NumericVector weights) {
    const int npts=neighbors.size();
    Rcpp::NumericMatrix output(nclusters, npts);
    Rcpp::NumericVector totals(npts);

    for (int i=0; i<npts; ++i) {
        auto COL=output.column(i);
        double& TOT=totals[i];
        Rcpp::IntegerVector NN=neighbors[i];

        for (auto it=NN.begin(); it!=NN.end(); ++it) {
            auto current=*it - 1;
            TOT += weights[current];
            COL[clusters[current]] += weights[current];
        }
    }

    return Rcpp::List::create(output, totals);
}
