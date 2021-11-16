
# Copyright (C) 2021         Dirk Eddelbuettel
#
# This file is part of RcppArmadillo.
#
# RcppArmadillo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# RcppArmadillo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

library(RcppArmadillo)

Rcpp::sourceCpp("cpp/fields_new.cpp")

.onWindows <- .Platform$OS.type == "windows"

f1m22 <- field1m22n()
expect_true(inherits(f1m22, "array"))
expect_true(inherits(f1m22[[1]], "matrix"))
expect_equal(dim(f1m22), c(1,1,1))
expect_equal(dim(f1m22[[1]]), c(2,2))

f11m22 <- field11m22n()
expect_true(inherits(f11m22, "array"))
expect_true(inherits(f11m22[[1]], "matrix"))
expect_equal(dim(f11m22), c(1,1,1))
expect_equal(dim(f11m22[[1]]), c(2,2))

f12m22 <- field12m22n()
expect_true(inherits(f12m22, "array"))
expect_true(inherits(f12m22[[1]], "matrix"))
expect_equal(dim(f12m22), c(1,2,1))
expect_equal(dim(f12m22[[1]]), c(2,2))
expect_equal(dim(f12m22[[2]]), c(2,2))

f21m22 <- field21m22n()
expect_true(inherits(f21m22, "array"))
expect_true(inherits(f21m22[[1]], "matrix"))
expect_equal(dim(f21m22), c(2,1,1))
expect_equal(dim(f21m22[[1]]), c(2,2))
expect_equal(dim(f21m22[[2]]), c(2,2))

f22m2233 <- field22m2233n()
expect_true(inherits(f22m2233, "array"))
expect_true(inherits(f22m2233[[1]], "matrix"))
expect_equal(dim(f22m2233), c(2,2,1))
expect_equal(dim(f22m2233[[2]]), c(3,3))
expect_equal(dim(f22m2233[[3]]), c(2,2))

f222m223344 <- field222m223344n()
expect_true(inherits(f222m223344, "array"))
expect_true(inherits(f222m223344[[1]], "matrix"))
expect_equal(dim(f222m223344), c(2,2,2))
expect_equal(dim(f222m223344[[2]]), c(3,3))
expect_equal(dim(f222m223344[[3]]), c(2,2))
expect_equal(dim(f222m223344[[5]]), c(4,4))


v <- infield1m22n( field1m22n() )
expect_equal(v, matrix(c(1L, 1L, 1L),3,1))

v <- infield11m22n( field11m22n() )
expect_equal(v, matrix(c(1L, 1L, 1L),3,1))

v <- infield12m22n( field12m22n() )
expect_equal(v, matrix(c(1L, 2L, 1L),3,1))

v <- infield21m22n( field21m22n() )
expect_equal(v, matrix(c(2L, 1L, 1L),3,1))

v <- infield22m2233n( field22m2233n() )
expect_equal(v, matrix(c(2L, 2L, 1L),3,1))

v <- infield222m223344n( field222m223344n() )
expect_equal(v, matrix(c(2L, 2L, 2L),3,1))
