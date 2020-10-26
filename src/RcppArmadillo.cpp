
// RcppArmadillo.cpp: Rcpp/Armadillo glue
//
// Copyright (C)  2010 - 2020  Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

//' Report the version of Armadillo
//'
//' @details The version is defined by Armadillo in the header \code{arma_version.hpp}.
//' @param single A logical vector indicating whether a single return values is requested,
//' or a named vector with three elements \code{major}, \code{minor} and \code{patch}.
//' @return Depending on the value of \code{single}, either a single number describing
//' the Armadillo version or a named vector with three elements \code{major}, \code{minor}
//' and \code{patch}.
//' @seealso Armadillo header file \code{arma_version.hpp}.
// [[Rcpp::export]]
Rcpp::IntegerVector armadillo_version(bool single) {

    // These are declared as constexpr in Armadillo which actually does not define them
    // They are also defined as macros in arma_version.hpp so we just use that
    const unsigned int major = ARMA_VERSION_MAJOR;
    const unsigned int minor = ARMA_VERSION_MINOR;
    const unsigned int patch = ARMA_VERSION_PATCH;

    if (single) {
        return Rcpp::wrap(10000 * major + 100 * minor + patch) ;
    } else {
        return Rcpp::IntegerVector::create(Rcpp::Named("major") = major,
                                           Rcpp::Named("minor") = minor,
                                           Rcpp::Named("patch") = patch);
    }
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
