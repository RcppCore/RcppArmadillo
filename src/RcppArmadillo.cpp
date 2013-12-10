// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// RcppArmadillo.cpp: Rcpp/Armadillo glue
//
// Copyright (C)  2010 - 2013  Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

using namespace Rcpp;

const unsigned int arma::arma_version::major;
const unsigned int arma::arma_version::minor;
const unsigned int arma::arma_version::patch;

// [[Rcpp::export]]
IntegerVector armadillo_version(bool single) {

    if( single ){
        return wrap( 10000*arma::arma_version::major +
            100*arma::arma_version::minor + 
            arma::arma_version::patch ) ;
    }

    IntegerVector version = IntegerVector::create(_["major"] = arma::arma_version::major,
                                                  _["minor"] = arma::arma_version::minor,
                                                  _["patch"] = arma::arma_version::patch);
    
   return version ;
}

