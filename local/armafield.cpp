
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
SEXP field2x2() {
    mat A = randn(2,3);
    mat B = randn(3,2);

    field<mat> F(2,2);
    F(0,1) = A;
    F(1,0) = B;

  F.print("F2x2:");
  return Rcpp::wrap(F);
}

// [[Rcpp::export]]
SEXP field2x3() {
    mat A = randn(2,3);
    mat B = randn(3,2);
    mat C = randn(2,4);

    field<mat> F(2,3);
    F(0,1) = A;
    F(1,0) = B;
    F(1,2) = C;

    F.print("F2x3:");
    return Rcpp::wrap(F);
}

// [[Rcpp::export]]
SEXP field2x3x2() {
    mat A = randn(2,3);
    mat B = randn(3,2);
    mat C = randn(2,4);

    field<mat> F(2,3,2);
    F(0,1,0) = A;
    F(1,0,1) = B;
    F(1,2,0) = C;

    F.print("F2x3x2:");
    return Rcpp::wrap(F);
}

/*** R
a1 <- field2x2();
print(a1[[3]])
a2 <- field2x3();
a3 <- field2x3x2();
*/
