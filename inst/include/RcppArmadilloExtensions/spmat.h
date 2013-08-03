// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// spmat.h: Conversion between Armadillo sp_mat and the dgCMatrix from Matrix
//
// Copyright (C)  2013  Dirk Eddelbuettel and Romain Francois
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

#ifndef RcppArmadillo__SpMat__h
#define RcppArmadillo__SpMat__h

#include <RcppArmadillo.h>

namespace Rcpp {

    
    // converts an SEXP object from R which was created as a sparse
    // matrix via the Matrix package) into an Armadillo sp_mat matrix
    //
    // TODO: template'ize to allow for types other than double, though
    //       realistically this is all we need
    template <> arma::sp_mat as(SEXP sx) {
        S4 mat(sx);  
        IntegerVector dims = mat.slot("Dim");
        arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
        arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));     
        arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));
        
        int nrow = dims[0], ncol = dims[1];
        arma::sp_mat res(nrow, ncol);

        // create space for values, and copy
        arma::access::rw(res.values) = arma::memory::acquire_chunked<double>(x.size() + 1);
        arma::arrayops::copy(arma::access::rwp(res.values), x.begin(), x.size() + 1);

        // create space for row_indices, and copy 
        arma::access::rw(res.row_indices) = 
            arma::memory::acquire_chunked<arma::uword>(i.size() + 1);
        arma::arrayops::copy(arma::access::rwp(res.row_indices), i.begin(), i.size() + 1);
    
        // create space for col_ptrs, and copy 
        arma::access::rw(res.col_ptrs) = arma::memory::acquire<arma::uword>(p.size() + 2);
        arma::arrayops::copy(arma::access::rwp(res.col_ptrs), p.begin(), p.size() + 1);

        // important: set the sentinel as well
        arma::access::rwp(res.col_ptrs)[p.size()+1] = std::numeric_limits<arma::uword>::max();
    
        // set the number of non-zero elements
        arma::access::rw(res.n_nonzero) = x.size();

        return res;
    }


    // convert an Armadillo sp_mat into a corresponding R sparse matrix
    // we copy to STL vectors as the Matrix package expects vectors whereas the
    // default wrap in Armadillo returns matrix with one row (or col) 
    SEXP wrap(arma::sp_mat sm) {

        IntegerVector dim(2);
        dim[0] = sm.n_rows; 
        dim[1] = sm.n_cols;

        arma::vec  x(sm.n_nonzero);        // create space for values, and copy
        arma::arrayops::copy(x.begin(), sm.values, sm.n_nonzero);
        std::vector<double> vx = arma::conv_to< std::vector< double > >::from(x);

        arma::urowvec i(sm.n_nonzero);	// create space for row_indices, and copy & cast
        arma::arrayops::copy(i.begin(), sm.row_indices, sm.n_nonzero);
        std::vector<int> vi = arma::conv_to< std::vector< int > >::from(i);
 
        arma::urowvec p(sm.n_cols+1);	// create space for col_ptrs, and copy 
        arma::arrayops::copy(p.begin(), sm.col_ptrs, sm.n_cols+1);
        // do not copy sentinel for returning R
        std::vector<int> vp = arma::conv_to< std::vector< int > >::from(p);

        S4 s("dgCMatrix");
        s.slot("i")   = vi;
        s.slot("p")   = vp;
        s.slot("x")   = vx;
        s.slot("Dim") = dim;
        return s;
    }

}

#endif
