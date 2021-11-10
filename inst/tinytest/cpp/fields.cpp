// fields.cpp: RcppArmadillo unit test code for field types
//
// Copyright (C) 2021         Dirk Eddelbuettel
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

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::field<arma::mat> field11m22() {
    arma::mat A = arma::randn(2,2);

    arma::field<arma::mat> F(1,1);
    F(0,0) = A;

    return F;
}

// [[Rcpp::export]]
arma::field<arma::mat> field12m22() {
    arma::mat A = arma::randn(2,2);
    arma::mat B = arma::randn(2,2);

    arma::field<arma::mat> F(1,2);
    F(0,0) = A;
    F(0,1) = B;

    return F;
}
