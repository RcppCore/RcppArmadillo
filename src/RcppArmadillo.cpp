
// RcppArmadillo.cpp: Rcpp/Armadillo glue
//
// Copyright (C)  2010 - 2025  Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

#include <RcppArmadillo/Lighter>

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
//' @rdname armadillo_set_seed
// [[Rcpp::export]]
void armadillo_set_seed_random() {
    arma::arma_rng::set_seed_random();  			// set the seed to a random value
}

//' @title Set the Armadillo Random Number Generator to given or random value
//' @description Setter functions for the internal Armadillo random number generator
//' @param val The seed used to initialize Armadillo's random number generator.
//' @details
//' Armadillo can switch between two random number generator implementations depending
//' on the compilation standard used. Under normal circumstances RcppArmadillo will connect
//' Armadillo to the R random number generator which also implies that \code{set.seed()}
//' should be used from R. To use this function, one also needs to undefine \code{ARMA_RNG_ALT}
//' so that the Armadillo generators are used.
//' @return The function is invoked for its side effect and has no return value.
//' @note This has been found to not work as expected in \pkg{RStudio}
//' as its code also uses the system RNG library. You may have to either
//' not run within \pkg{RStudio} or change your code to use a different RNG such
//' as the one from R.
//' @seealso The R documentation on its RNGs all of which are accessible via \pkg{Rcpp}.
// [[Rcpp::export]]
void armadillo_set_seed(unsigned int val) {
    arma::arma_rng::set_seed(val);  			    // set the seed to given value
}

//' Report (or Set) Maximum Number of OpenMP Threads
//'
//' @param n Number of threads to be set
//' @return For the getter, and on a system with OpenMP, the maximum
//' number of threads that OpenMP may be using and on systems without it,
//' one. The setter does not return a value.
// [[Rcpp::export]]
int armadillo_get_number_of_omp_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

//' @rdname armadillo_get_number_of_omp_threads
// [[Rcpp::export]]
void armadillo_set_number_of_omp_threads(int n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)(n);                  // prevent unused variable warning
#endif
}
