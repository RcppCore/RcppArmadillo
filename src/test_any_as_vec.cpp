// [ [ Rcpp::depends(RcppArmadillo) ] ]
#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export(.veccany_as_v_test)]]
arma::vec veccany_as_v_test(arma::vec v) { return(v); }
// [[Rcpp::export(.vecrany_as_v_test)]]
arma::rowvec vecrany_as_v_test(arma::rowvec v) { return(v); }
