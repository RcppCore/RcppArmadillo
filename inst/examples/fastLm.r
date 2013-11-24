#!/usr/bin/r
##
## fastLm.r: Benchmarking lm() via RcppArmadillo and directly
##
## Copyright (C)  2010 - 2013  Dirk Eddelbuettel, Romain Francois and Douglas Bates
##
## This file is part of RcppArmadillo.
##
## RcppArmadillo is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppArmadillo is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

library(RcppArmadillo)
library(rbenchmark)

src <- '
Rcpp::List fLmTwoCasts(Rcpp::NumericMatrix Xr, Rcpp::NumericVector yr) {
    int n = Xr.nrow(), k = Xr.ncol();
    arma::mat X(Xr.begin(), n, k, false);
    arma::colvec y(yr.begin(), yr.size(), false);
    int df = n - k;

    // fit model y ~ X, extract residuals
    arma::colvec coef = arma::solve(X, y);
    arma::colvec res  = y - X*coef;

    double s2 = std::inner_product(res.begin(), res.end(),
                                   res.begin(), 0.0)/df;
    // std.errors of coefficients
    arma::colvec sderr = arma::sqrt(s2 *
       arma::diagvec(arma::pinv(arma::trans(X)*X)));

    return Rcpp::List::create(Rcpp::Named("coefficients")=coef,
                              Rcpp::Named("stderr")      =sderr,
                              Rcpp::Named("df")          =df);
}
'
cppFunction(code=src, depends="RcppArmadillo")

src <- '
Rcpp::List fLmOneCast(arma::mat X, arma::colvec y) {
    int df = X.n_rows - X.n_cols;

    // fit model y ~ X, extract residuals
    arma::colvec coef = arma::solve(X, y);
    arma::colvec res  = y - X*coef;

    double s2 = std::inner_product(res.begin(), res.end(),
                                   res.begin(), 0.0)/df;
    // std.errors of coefficients
    arma::colvec sderr = arma::sqrt(s2 *
       arma::diagvec(arma::pinv(arma::trans(X)*X)));

    return Rcpp::List::create(Rcpp::Named("coefficients")=coef,
                              Rcpp::Named("stderr")      =sderr,
                              Rcpp::Named("df")          =df);
}
'
cppFunction(code=src, depends="RcppArmadillo")

src <- '
Rcpp::List fLmConstRef(const arma::mat & X, const arma::colvec & y) {
    int df = X.n_rows - X.n_cols;

    // fit model y ~ X, extract residuals
    arma::colvec coef = arma::solve(X, y);
    arma::colvec res  = y - X*coef;

    double s2 = std::inner_product(res.begin(), res.end(),
                                   res.begin(), 0.0)/df;
    // std.errors of coefficients
    arma::colvec sderr = arma::sqrt(s2 *
       arma::diagvec(arma::pinv(arma::trans(X)*X)));

    return Rcpp::List::create(Rcpp::Named("coefficients")=coef,
                              Rcpp::Named("stderr")      =sderr,
                              Rcpp::Named("df")          =df);
}
'
cppFunction(code=src, depends="RcppArmadillo")


fastLmPureDotCall <- function(X, y) {
    .Call("RcppArmadillo_fastLm", X, y, PACKAGE = "RcppArmadillo")
}


y <- log(trees$Volume)
X <- cbind(1, log(trees$Girth))
frm <- formula(log(Volume) ~ log(Girth))

res <- benchmark(fLmOneCast(X, y),             	# inline'd above
                 fLmTwoCasts(X, y),            	# inline'd above
                 fLmConstRef(X, y),            	# inline'd above
                 fastLmPure(X, y),              # similar, but with 2 error checks
                 fastLmPureDotCall(X, y),       # now without the 2 error checks
                 fastLm(frm, data=trees),       # using model matrix
                 lm.fit(X, y),                  # R's fast function, no stderr
                 lm(frm, data=trees),           # R's standard function
                 columns = c("test", "replications", "relative",
                             "elapsed", "user.self", "sys.self"),
                 order="relative",
                 replications=5000)

print(res[,1:4])
