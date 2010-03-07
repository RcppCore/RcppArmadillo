// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// fastLm.cpp: Rcpp/Armadillo glue example of a simple lm() alternative
//
// Copyright (C)  2010 Dirk Eddelbuettel and Romain Francois
//
// This file is part of RcppArmadillo.
//
// RcppArmadillo is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppArmadillo is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>

extern "C" SEXP fastLm(SEXP ys, SEXP Xs) {

    Rcpp::NumericVector yr(ys);			// creates Rcpp vector from SEXP
    Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
    int n = Xr.nrow(), k = Xr.ncol();

#if ARMA_VERSION_GE_090
    arma::mat X(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
    arma::colvec y(yr.begin(), yr.size(), false);
#else
    arma::mat X(Xr.begin(), n, k);
    arma::colvec y(yr.begin(), yr.size());
#endif

    arma::colvec coef = solve(X, y);            // fit model y ~ X

    arma::colvec resid = y - X*coef; 		// residuals

#if ARMA_VERSION_GE_090
    double sig2 = arma::as_scalar( trans(resid)*resid/(n-k) );
#else
    double sig2 = ( trans(resid)*resid/(n-k) );
#endif
    						// std.error of estimate 
    arma::colvec stderrest = sqrt( sig2 * diagvec( arma::inv(arma::trans(X)*X)) );

    Rcpp::Pairlist res(Rcpp::Named( "coefficients", coef),
                       Rcpp::Named( "stderr", stderrest));
    return res;
}

