// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// fastLm.cpp: Rcpp/Armadillo glue example of a simple lm() alternative
//
// Copyright (C)  2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

    try {
	arma::colvec y = Rcpp::as<arma::colvec>(ys);	// direct from SEXP to arma::mat
	arma::mat X    = Rcpp::as<arma::mat>(Xs);
	int df = X.n_rows - X.n_cols;

	arma::colvec coef = arma::solve(X, y);      	// fit model y ~ X
	arma::colvec res  = y - X*coef;			// residuals

	double s2 = std::inner_product(res.begin(), res.end(), res.begin(), double())/df;
							// std.errors of coefficients
	arma::colvec std_err = arma::sqrt(s2 * arma::diagvec( arma::pinv(arma::trans(X)*X) ));	

	return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
				  Rcpp::Named("stderr")       = std_err,
				  Rcpp::Named("df")           = df
				  );

    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}





