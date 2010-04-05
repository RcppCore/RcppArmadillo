## Copyright (C)       2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

RcppArmadilloCxxFlags <- function(Rcpp = TRUE, ...){
    rcpp <- ifelse(Rcpp, Rcpp:::RcppCxxFlags(...), "")

    arma <- system.file("include", package = "RcppArmadillo")
    if (.Platform$OS.type=="windows")
        arma <- paste('"', arma, '"', sep="")
    
    res <- sprintf('%s -I%s -I/usr/include ', rcpp, arma)
}

CxxFlags <- function(Rcpp = TRUE, ...){
    cat( RcppArmadilloCxxFlags(Rcpp = Rcpp, ...), sep = " " )
}


RcppArmadilloLdFlags <- function(Rcpp = TRUE, ...){
    rcpp <- ifelse( Rcpp, Rcpp:::RcppLdFlags(static= !grepl("^linux",R.version$os) ) , "" )
    arma <- c("-L/usr/lib -larmadillo" )
    paste( rcpp, arma, "", sep = " " )
}

LdFlags <- function(Rcpp = TRUE, ...){
    cat( RcppArmadilloLdFlags(Rcpp, ... ), sep = "" )
}
## no ldflags as user packages would not need to link to RcppArmadillo.so
