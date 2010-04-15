
## fastLm.R: Rcpp/Armadillo implementation of lm()
##
## Copyright (C)  2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

fastLm <- function(y, X) {

    stopifnot(is.matrix(X))
    stopifnot(nrow(y)==nrow(X))

    res <- .Cpp("fastLm", y, X, package="RcppArmadillo")
}

## What would be nice here:
##
##   fastLm <- function(x, ...) UseMethod("fastLm")
##
##   fastLm.formula <- ...
##
##   fastLm.default <- ...
##
##   print.fastLm <- function(x, ...)
##
##   summary.fastLm <- function(object, ...
##
##   print.summary.fastLm <- ...
##
## but on the other hand lm.fit() does not do any of this either
