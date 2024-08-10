#include "Rcpp.h"
#include "scran_graph_cluster/build_snn_graph.hpp"
#include "tatami/tatami.hpp"

// [[Rcpp::export(rng=false)]]
Rcpp::List build_snn_graph(Rcpp::IntegerMatrix neighbors, std::string scheme) {
    const int* nptr = neighbors.begin();
    size_t nrow = neighbors.rows();
    scran_graph_cluster::BuildSnnGraphResults<int, double> output;

    scran_graph_cluster::BuildSnnGraphOptions opt;
    if (scheme == "rank") {
        opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::RANKED;
    } else if (scheme == "number") {
        opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::NUMBER;
    } else if (scheme == "jaccard") {
        opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::JACCARD;
    } else {
        throw std::runtime_error("unknown weighting scheme '" + scheme + "'");
    }

    size_t ncells = neighbors.cols();
    scran_graph_cluster::build_snn_graph(
        ncells,
        [&](int i) -> tatami::ArrayView<int> {
            return tatami::ArrayView<int>(nptr + nrow * static_cast<size_t>(i), nrow);
        },
        [](int i) -> int {
            return i - 1;
        },
        opt,
        output
    );

    size_t nedges = output.edges.size();
    Rcpp::IntegerVector edges(nedges);
    int* eptr = edges.begin();
    for (size_t e = 0; e < nedges; ++e) {
        eptr[e] = output.edges[e] + 1; // get to 1-based indexing.
    }

    return Rcpp::List::create(
        Rcpp::Named("edges") = edges,
        Rcpp::Named("weights") = Rcpp::NumericVector(output.weights.begin(), output.weights.end())
    );
}
