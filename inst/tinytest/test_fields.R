
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

Rcpp::sourceCpp("cpp/fields.cpp")

.onWindows <- .Platform$OS.type == "windows"

f11m22 <- field11m22()
expect_true(inherits(f11m22, "array"))
expect_true(inherits(f11m22[[1]], "matrix"))
expect_equal(dim(f11m22), c(1,1))
expect_equal(dim(f11m22[[1]]), c(2,2))

f12m22 <- field12m22()
expect_true(inherits(f12m22, "array"))
expect_true(inherits(f12m22[[1]], "matrix"))
expect_equal(dim(f12m22), c(1,2))
expect_equal(dim(f12m22[[1]]), c(2,2))
expect_equal(dim(f12m22[[2]]), c(2,2))
