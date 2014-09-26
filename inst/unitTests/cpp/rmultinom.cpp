// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/rmultinom.h>
#include <RcppArmadilloExtensions/fixprob.h>

using namespace Rcpp ;
// [[Rcpp::export]]
IntegerVector rmultinomC( int n, int size,
    NumericVector prob 
) {
    IntegerMatrix draws(prob.size(), n);
    // FixProb modifies in-place
    NumericVector fixprob = clone(prob);
    RcppArmadillo::FixProb(fixprob, 1, true);
    RNGScope scope;
    for (int ii=0; ii<n; ii++){
        draws(_,ii) = RcppArmadillo::rmultinom( size, fixprob);
    }
    return draws;
}
