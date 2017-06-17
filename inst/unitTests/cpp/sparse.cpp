// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// sparse.cpp: RcppArmadillo unit test code for sparse matrices 
//
// Copyright (C) 2014  Dirk Eddelbuettel
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
arma::sp_mat asSpMat(SEXP S) {
    return Rcpp::as<arma::sp_mat>(S);
}

// [[Rcpp::export]]
arma::sp_mat sparseAddition(arma::sp_mat SM) {
    return SM + SM;
}

// [[Rcpp::export]]
arma::sp_mat sparseMultiplication(arma::sp_mat SM, int k) {
    return k * SM;
}

// [[Rcpp::export]]
arma::sp_mat fromTriplet(arma::urowvec ri, arma::urowvec ci, arma::colvec values) {
    arma::umat loc = arma::join_vert(ri, ci);// form 2*N 'locations' matrix
    arma::sp_mat sp(loc, values);            // create sparse from locations and values
    return sp;
}

// [[Rcpp::export]]
arma::sp_mat sparseTranspose(arma::sp_mat SM) {
    return SM.t();
}

// [[Rcpp::export]]
arma::sp_mat sparseSqrt(arma::sp_mat SM) {
    return arma::sqrt(SM);
}

// [[Rcpp::export]]
arma::sp_mat sparseSquare(arma::sp_mat SM) {
    return arma::square(SM);
}

// [[Rcpp::export]]
arma::sp_mat sparseIterators(arma::sp_mat SM, double val) {
    arma::sp_mat::iterator begin = SM.begin();
    arma::sp_mat::iterator end   = SM.end();
    
    for (arma::sp_mat::iterator it = begin; it != end; ++it)
      (*it) += val;
    
    return SM;
}

// [[Rcpp::export]]
Rcpp::List sparseList(Rcpp::List l) {
    arma::sp_mat mat1 = l[0];
    arma::sp_mat mat2 = l[0];
    
    return Rcpp::List::create(mat1, mat2);
}

// [[Rcpp::export]]
arma::sp_mat dtc2dgc(SEXP S){
    S4 mat = S;
    IntegerVector dims = mat.slot("Dim");
    int nrow = dims[0];
    int ncol = dims[1];
    
    // Creating an empty SpMat
    arma::sp_mat res(static_cast<unsigned>(nrow), static_cast<unsigned>(ncol));
    
    // Get the type of sparse matrix
    std::string type = Rcpp::as<std::string>(mat.slot("class"));
    if (type == "dtCMatrix") {
        IntegerVector i = mat.slot("i");
        IntegerVector p = mat.slot("p");
        IntegerVector x = mat.slot("x");
        std::string diag = Rcpp::as<std::string>(mat.slot("diag"));
        
        // Making space for the elements
        res.mem_resize(static_cast<unsigned>(x.size()));
        
        // Copying elements
        std::copy(i.begin(), i.end(), arma::access::rwp(res.row_indices));
        std::copy(p.begin(), p.end(), arma::access::rwp(res.col_ptrs));
        std::copy(x.begin(), x.end(), arma::access::rwp(res.values));
        
        if (diag == "U") {
            res.diag().ones();
        }
    }
    
    // Setting the sentinel
    arma::access::rw(res.col_ptrs[static_cast<unsigned>(ncol + 1)]) =
    std::numeric_limits<arma::uword>::max();
    
    return res;
}
