// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// build_snn_rank
Rcpp::List build_snn_rank(Rcpp::IntegerMatrix neighbors);
RcppExport SEXP _bluster_build_snn_rank(SEXP neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type neighbors(neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_snn_rank(neighbors));
    return rcpp_result_gen;
END_RCPP
}
// build_snn_number
Rcpp::List build_snn_number(Rcpp::IntegerMatrix neighbors);
RcppExport SEXP _bluster_build_snn_number(SEXP neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type neighbors(neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_snn_number(neighbors));
    return rcpp_result_gen;
END_RCPP
}
// flow_som
Rcpp::List flow_som(Rcpp::NumericMatrix data, Rcpp::NumericMatrix original_codes, Rcpp::NumericMatrix nhbrdist, Rcpp::NumericVector alphas, Rcpp::NumericVector radii, int rlen, int dist);
RcppExport SEXP _bluster_flow_som(SEXP dataSEXP, SEXP original_codesSEXP, SEXP nhbrdistSEXP, SEXP alphasSEXP, SEXP radiiSEXP, SEXP rlenSEXP, SEXP distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type original_codes(original_codesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type nhbrdist(nhbrdistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type radii(radiiSEXP);
    Rcpp::traits::input_parameter< int >::type rlen(rlenSEXP);
    Rcpp::traits::input_parameter< int >::type dist(distSEXP);
    rcpp_result_gen = Rcpp::wrap(flow_som(data, original_codes, nhbrdist, alphas, radii, rlen, dist));
    return rcpp_result_gen;
END_RCPP
}
// sum_neighbor_weights
Rcpp::List sum_neighbor_weights(int nclusters, Rcpp::List neighbors, Rcpp::IntegerVector clusters, Rcpp::NumericVector weights);
RcppExport SEXP _bluster_sum_neighbor_weights(SEXP nclustersSEXP, SEXP neighborsSEXP, SEXP clustersSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type nclusters(nclustersSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type neighbors(neighborsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type clusters(clustersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_neighbor_weights(nclusters, neighbors, clusters, weights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bluster_build_snn_rank", (DL_FUNC) &_bluster_build_snn_rank, 1},
    {"_bluster_build_snn_number", (DL_FUNC) &_bluster_build_snn_number, 1},
    {"_bluster_flow_som", (DL_FUNC) &_bluster_flow_som, 7},
    {"_bluster_sum_neighbor_weights", (DL_FUNC) &_bluster_sum_neighbor_weights, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bluster(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
