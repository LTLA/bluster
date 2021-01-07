/* flow_som.cpp
 *
 * This is based on som.c in the FlowSOM package,
 * which in turn is based on code of Ron Wehrens
 */

#include <cmath>
#include "Rcpp.h"

double eucl2(const double * p1, const double * p2, int n) {
    double xdist = 0.0;
    for (int j = 0; j < n; ++j, ++p1, ++p2) {
        double tmp = *p1 - *p2;
        xdist += tmp * tmp;
    }
    return xdist; // save some time by avoiding a sqrt.
}

double manh(const double * p1, const double * p2, int n) {
    double xdist = 0.0;
    for (int j = 0; j < n; ++j, ++p1, ++p2) {
        double tmp = *p1 - *p2;
        xdist += std::abs(tmp);
    }
    return xdist;
}

double chebyshev(const double * p1, const double * p2, int n) {
    double xdist = 0.0;
    for (int j = 0; j < n; ++j, ++p1, ++p2) {
        double tmp = std::abs(*p1 - *p2);
        if (tmp > xdist) {
            xdist = tmp;
        }
    }
    return xdist;
}

double cosine(const double * p1, const double * p2, int n) {
    double nom = 0;
    double denom1 = 0;
    double denom2 = 0;
    for (int j = 0; j < n; ++j, ++p1, ++p2) {
        nom += (*p1) * (*p2);
        denom1 += (*p1) * (*p1);
        denom2 += (*p2) * (*p2);
    }
    return 1 - nom/(std::sqrt(denom1 * denom2));
}

//[[Rcpp::export]]
Rcpp::List flow_som(Rcpp::NumericMatrix data,
    Rcpp::NumericMatrix original_codes,
    Rcpp::NumericMatrix nhbrdist,
    Rcpp::NumericVector alphas, 
    Rcpp::NumericVector radii,
    int rlen,
    int dist)
{
    Rcpp::NumericMatrix codes = Rcpp::clone(original_codes);

    // Transposed; rows are dimensions, columns are samples/codes.
    // This gives us better cache efficiency.
    const int n = data.ncol();
    const int ncodes = codes.ncol();
    const int px = data.nrow();

    double (*distf)(const double*,const double*,int);
    if(dist == 1){
        distf = &manh;
    } else if (dist == 2){
        distf = &eucl2;
    } else if (dist == 3){
        distf = &chebyshev;
    } else if (dist == 4){
        distf = &cosine;
    } else {
        distf = &eucl2;
    }

    double change=0;
    const int niter = std::round(rlen * n);
    double threshold = radii[0];
    double thresholdStep = (radii[0] - radii[1]) / static_cast<double>(niter);

    for (int k = 0; k < niter; ++k) {
        if(k%n == 0){
            if(change < 1){
                k = niter;
            }
            change = 0.0;
        }    
    
        int i = static_cast<int>(n * R::unif_rand()); /* Select a random sample */
    
        int nearest = 0;
        double neardist = R_PosInf;
        auto curobs = data.column(i);

        /* calculate distances in x and y spaces, and keep track of the
        nearest code */
        for (int cd = 0; cd < ncodes; ++cd) {
            auto curcode = codes.column(cd);
            double curdist = distf(curobs.begin(), curcode.begin(), px);
            if (curdist < neardist) {
                neardist = curdist;
                nearest = cd;
            }
        }

        if (threshold < 1.0) {
            threshold = 0.5;
        }
        double alpha = alphas[0] - (alphas[0] - alphas[1]) * static_cast<double>(k)/static_cast<double>(niter);

        for (int cd = 0; cd < ncodes; ++cd) {
            if (nhbrdist.column(cd)[nearest] > threshold) {
                continue;
            }

            auto coIt = curobs.begin();
            auto curcode = codes.column(cd);
            auto ccIt = curcode.begin();

            for(int j = 0; j < px; ++j, ++coIt, ++ccIt) {
                double tmp = *coIt - *ccIt;
                change += std::abs(tmp);
                codes[cd + j*ncodes] += tmp * alpha; 
            }
        }
    
        threshold -= thresholdStep;
    }

    // Mapping everyone to their closest code.
    Rcpp::IntegerVector assigned(n);
    auto aIt = assigned.begin();

    for (int i = 0; i < n; ++i, ++aIt) {
        int nearest = 0;
        double neardist = R_PosInf;
        auto curobs = data.column(i);

        for (int cd = 0; cd < ncodes; ++cd) {
            auto curcode = codes.column(cd);
            double curdist = distf(curobs.begin(), curcode.begin(), px);
            if (curdist < neardist) {
                neardist = curdist;
                nearest = cd;
            }
        }

        *aIt = nearest;
    }

    return Rcpp::List::create(codes, assigned); 
}
