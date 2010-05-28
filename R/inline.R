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

Makevars.RcppArmadillo <- '
PKG_LIBS = $(shell $(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()" ) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
'
Makevars.win.RcppArmadillo <- '
PKG_LIBS = $(shell $(R_HOME)/bin${R_ARCH_BIN}/Rscript.exe -e "Rcpp:::LdFlags()") $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
'

inline_cxx_plugin <- Rcpp:::Rcpp.plugin.maker(
	include.before = "#include <RcppArmadillo.h>", 
	LinkingTo = c("Rcpp", "RcppArmadillo"), 
	Depends = c("Rcpp", "RcppArmadillo"),
	libs = "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)", 
	Makevars = Makevars.RcppArmadillo, 
	Makevars.win = Makevars.win.RcppArmadillo
)
