#!/usr/bin/r -t
#
# Copyright (C) 2010	Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

test.wrap.R <- function(){
	res <- .Call( "RcppArmadillo_wrap", PACKAGE = "RcppArmadillo" )

	checkEquals( res[[1]][[1]], matrix(as.integer((diag(3))),nr=3), msg = "eye<imat>(3,3)" )
	checkEquals( res[[1]][[2]], diag(3), msg = "eye<mat>(3,3)" )
	checkEquals( res[[1]][[3]], diag(3), msg = "eye<fmat>(3,3)" )
	checkEquals( res[[1]][[4]], matrix(as.integer((diag(3))),nr=3), msg = "eye<umat>(3,3)" )

	checkEquals( res[[2]][[1]], matrix(0, ncol = 5, nrow=1), msg = "zeros<mat>(5,1)" )
	checkEquals( res[[2]][[2]], matrix(0, ncol = 5, nrow=1), msg = "zeros<fmat>(5,1)" )

	checkEquals( res[[3]][[1]], matrix(0, ncol = 1, nrow=5), msg = "zeros<mat>(1,5)" )
	checkEquals( res[[3]][[2]], matrix(0, ncol = 1, nrow=5), msg = "zeros<mat>(1,5)" )

	checkEquals( res[[4]][[1]], matrix(0:3, ncol = 2, nrow=2), msg = "field<int>" )
	checkEquals( res[[4]][[2]], matrix(letters[1:4], ncol = 2, nrow=2), msg = "field<std::string>" )
}

test.wrap.Glue <- function(){
	res <- .Call( "RcppArmadillo_wrap_Glue", PACKAGE = "RcppArmadillo" )
	checkEquals( res[[1]], 2*diag(3), msg = "wrap(Glue)" )
}

test.wrap.Op <- function(){
	res <- .Call( "RcppArmadillo_wrap_Op", PACKAGE = "RcppArmadillo" )
	checkEquals( res[[1]], -1*diag(3), msg = "wrap(Op)" )
}

test.as.Mat <- function(){

	integer_mat <- matrix( as.integer(diag(4)), ncol = 4, nrow = 4 )
	numeric_mat <- diag(5)
	res <- .Call( "RcppArmadillo_as_Mat",
		list( integer_mat, numeric_mat ),
		PACKAGE = "RcppArmadillo" )
	checkEquals( unlist( res), c(4L, 5L, 4L, 5L ), msg = "as<Mat>" )
}

test.as.Col <- function(){
	res <- .Call( "RcppArmadillo_as_Col",
	      list( 1:10, as.numeric(1:10) ),
	      PACKAGE = "RcppArmadillo" )
	checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Col>" )
}

test.as.Row <- function(){
	res <- .Call( "RcppArmadillo_as_Row",
	      list( 1:10, as.numeric(1:10) ),
	      PACKAGE = "RcppArmadillo" )
	checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Row>" )
}

test.fastLm <- function() {
    data(trees)
    flm <- .Call("fastLm",
                 log(trees$Volume),
                 cbind(rep(1,31), log(trees$Girth)),
                 PACKAGE="RcppArmadillo")
    fit <- lm(log(Volume) ~ log(Girth), data=trees)

    checkEquals( as.numeric(flm$coef), as.numeric(coef(fit)),
                msg="fastLm.coef")
    checkEquals( as.numeric(flm$stderr), as.numeric(coef(summary(fit))[,2]),
                msg="fastLm.stderr")
}

