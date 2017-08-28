## Copyright (C)       2010 - 2017  Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

inlineCxxPlugin <- function(...) {
    ismacos <- Sys.info()[["sysname"]] == "Darwin"
    openmpflag <- if (ismacos) "" else "$(SHLIB_OPENMP_CFLAGS)"
    plugin <-
        Rcpp::Rcpp.plugin.maker(
                  include.before = "#include <RcppArmadillo.h>",
                  libs           = paste(openmpflag, "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"),
                  package        = "RcppArmadillo"
              )
    settings <- plugin()
    settings$env$PKG_CPPFLAGS <- paste("-I../inst/include", openmpflag)
    if (!ismacos) settings$env$USE_CXX11 <- "yes"
    settings
}

