
## Copyright (C)       2010 Dirk Eddelbuettel and Romain Francois
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

.__CxxFlags <- function(Rcpp = TRUE, ...){
    rcpp <- ifelse(Rcpp, Rcpp:::RcppCxxFlags(...), "")

    arma <- system.file("include", package = "RcppArmadillo")
    if (.Platform$OS.type=="windows")
        arma <- paste('"', arma, '"', sep="")

    res <- sprintf('%s -I%s', rcpp, arma)
}

.__LdFlags <- function(Rcpp = TRUE, ...){
    rcpp <- ifelse( Rcpp, Rcpp:::LdFlags(...) , "" )
    arma <- c("-larmadillo" )
    paste( rcpp, arma, sep = " " )
}

CxxFlags <- function(Rcpp = TRUE, ...){
    cat( .__CxxFlags(Rcpp = Rcpp, ...), sep = " " )
}

LdFlags <- function(Rcpp = TRUE, ...){
    cat( .__LdFlags(Rcpp, ... ), sep = "" )
}
## no ldflags as user packages would not need to link to RcppArmadillo.so
