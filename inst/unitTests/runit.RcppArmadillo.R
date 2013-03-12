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

.setUp <- function(){
    suppressMessages(require(inline))
    sourceCpp(file.path(pathRcppArmadilloTests, "cpp/armadillo.cpp"))
}

test.wrap.R <- function(){
    fx <- wrap_
    res <- fx()

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
    fx <- wrapGlue_
    res <- fx()
    checkEquals( res[[1]], 2*diag(3), msg = "wrap(Glue)" )
}

test.wrap.Op <- function(){
    fx <- wrapOp_
    res <- fx()
    checkEquals( res[[1]], -1*diag(3), msg = "wrap(Op)" )
}

test.as.Mat <- function(){
    fx <- asMat_
    integer_mat <- matrix( as.integer(diag(4)), ncol = 4, nrow = 4 )
    numeric_mat <- diag(5)
    res <- fx( list( integer_mat, numeric_mat ) )
    checkEquals( unlist( res), c(4L, 5L, 4L, 5L ), msg = "as<Mat>" )
}

test.as.Col <- function(){
    fx <- asCol_
    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Col>" )
}

test.as.Row <- function(){
    fx <- asRow_
    res <- fx( list( 1:10, as.numeric(1:10) ) )
    checkEquals( unlist( res ), rep(55.0, 4 ), msg = "as<Row>" )
}

test.cxmat <- function(){
    fx <- cxMat_
    checkEquals(fx(),
		list( double = (1+0i)*diag(3), float = (1+0i)*diag(3) ),
		msg = "support for complex matrices" )

}

test.mtOp <- function(){
    fx <- mtOp_
    checkEquals(fx(),
		(1+2i)*diag(3),
		msg = "support for mtOp" )

}

test.mtGlue <- function(){
    fx <- mtGlue_
    checkEquals(fx(),
		2.0 * diag(3) ,
		msg = "support for mtGlue" )

}


test.sugar <- function(){
    fx <- sugar_
    checkEquals(fx(1:10),
		matrix( 2*(1:10), nrow = 10 ) ,
		msg = "RcppArmadillo and sugar" )

}

test.sugar.cplx <- function(){
    fx <- sugarCplx_
    x <- 1:10*(1+1i)
    checkEquals(fx(x),
		matrix( exp(x), nrow = 10 ) ,
		msg = "RcppArmadillo and sugar (complex)" )

}

test.armadillo.sugar.ctor <- function(){
    fx <- sugarCtor_
    checkEquals(fx(1:10),
		list(mat = matrix( 4*(1:10), nrow = 10 ),
                     rowvec = matrix( 1:10, nrow = 1 ),
                     colvec = matrix( 1:10, ncol = 1 )) ,
		msg = "Mat( sugar expression )" )

}


test.armadillo.sugar.matrix.ctor <- function(){
    fx <- sugarMatrixCtor_
    res <- fx(1:10)
    norm <- function(x, y) sqrt( x*x + y*y )
    checkEquals(res,
		list(mat = diag(2*(1:10)),
                     rowvec = outer( 1, 1:10, norm ),
                     colvec = outer( 1:10, 1, norm )),
		msg = "Mat( sugar expression )" )

}

## test.armadillo.rtti.check <- function() {

##     inc <- '
##     void blah(arma::mat& X) {
##          X.set_size(5,5);
##     }
##     '
##     src <- '
##     arma::vec V;
##     blah(V); // if blah() worked, we have a problem
##     '
##     fun <- cxxfunction(signature(), body=src, inc=inc, plugin = "RcppArmadillo")

##     checkException(fun(), msg="RTTI check on matrix constructor exception")

## }
