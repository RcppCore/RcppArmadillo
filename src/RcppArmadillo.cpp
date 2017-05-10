// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// RcppArmadillo.cpp: Rcpp/Armadillo glue
//
// Copyright (C)  2010 - 2014  Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

//' Report the version of Armadillo 
//' 
//' @details The version is defined by Armadillo in the header \code{arma_version.hpp}.
//' @param single A logical vector indicating whether a single return values is requested,
//' or a named vector with three elements \code{major}, \code{minor} and \code{patch}.
//' @return Depending on the value of \code{single}, either a single number describing the Armadillo version
//' or a named vector with three elements \code{major}, \code{minor} and \code{patch}.
//' @seealso Armadillo header file \code{arma_version.hpp}.
// [[Rcpp::export]]
IntegerVector armadillo_version(bool single) {

    if (single) {
        return wrap(10000*arma::arma_version::major +
                    100*arma::arma_version::minor + 
                    arma::arma_version::patch) ;
    }

    IntegerVector version = IntegerVector::create(_["major"] = arma::arma_version::major,
                                                  _["minor"] = arma::arma_version::minor,
                                                  _["patch"] = arma::arma_version::patch);
    
    return version ;
}

// Per request of Gábor Csárdi in https://github.com/RcppCore/RcppArmadillo/issues/11
//
//' Set the Armadillo Random Number Generator to a random value
//' 
//' @details
//' Depending on whether RcppArmadillo was compiled for the C++98 standard
//' (currently the default) or for C++11 (optional), two different RNGs may be used.
//' This function resets either. For C++98, the R programming language's RNG is used.
//' For C++11, the RNG included in the \code{<random>} library is used only when
//'  \code{#define ARMA_USE_CXX11_RNG} is placed before \code{#include <RcppArmadillo.h>}.
//'  Otherwise, the R programming language's RNG will be used.
//' @return The function is invoked for its side effect and has no return value.
//' @note This has been found to not work as espected in \pkg{RStudio}
//' as its code also uses the system RNG library. You may have to either
//' not run within \pkg{RStudio} or change your code to use a different RNG such 
//' as the one from R.
//' @seealso The R documentation on its RNGs all of which are accessible via \pkg{Rcpp}.  
// [[Rcpp::export]]
void armadillo_set_seed_random() {
    arma::arma_rng::set_seed_random();  			// set the seed to a random value
} 

//' Set the Armadillo Random Number Generator to the given value
//' 
//' @param val The seed used to initialize Armadillo's random number generator.
//' @details 
//' Depending on whether RcppArmadillo was compiled for the C++98 standard
//' (currently the default) or for C++11 (optional), two different RNGs may be used.
//' This function resets either. For C++98, the R programming language's RNG is used.
//' For C++11, the RNG included in the \code{<random>} library is used only when
//'  \code{#define ARMA_USE_CXX11_RNG} is placed before \code{#include <RcppArmadillo.h>}.
//'  Otherwise, the R programming language's RNG will be used.
//' @return The function is invoked for its side effect and has no return value. 
//' @note This has been found to not work as espected in \pkg{RStudio}
//' as its code also uses the system RNG library. You may have to either
//' not run within \pkg{RStudio} or change your code to use a different RNG such 
//' as the one from R.
//' @seealso The R documentation on its RNGs all of which are accessible via \pkg{Rcpp}.  
// [[Rcpp::export]]
void armadillo_set_seed(unsigned int val) {
    //Rcpp::Rcout << "Setting value " << val << std::endl;
    arma::arma_rng::set_seed(val);  			    // set the seed to given value
} 
