## init.R: Startup
##
## Copyright (C)  2023  Dirk Eddelbuettel
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

.pkgenv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {
    .pkgenv[["omp_threads"]] <- armadillo_get_number_of_omp_threads()
}

##' Throttle (or Reset) (Rcpp)Armadillo to Two Cores
##'
##' Helper functions to throttle use of cores by RcppArmadillo-internal
##' code on systems with OpenMP. On package load, the initial value is
##' saved and used to reset the value.
##' @param n Integer value of desired cores, default is two
armadillo_throttle_cores <- function(n = 2) {
    armadillo_set_number_of_omp_threads(n)
}

##' @rdname armadillo_throttle_cores
armadillo_reset_cores <- function() {
    n <- .pkgenv[["omp_threads"]]
    armadillo_set_number_of_omp_threads(n)
}
