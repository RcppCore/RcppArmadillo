// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppArmadilloDefines.h: Configure-generated defines for Rcpp/Armadillo glue
//
// Copyright (C)  2010 Dirk Eddelbuettel and Romain Francois
//
// This file is part of RcppArmadillo.
//
// RcppArmadillo is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppArmadillo is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RcppArmadillo__RcppArmadilloDefines__H
#define RcppArmadillo__RcppArmadilloDefines__H

//#define ARMA_VERSION_GE_070 @ARMA_VERSION_GE_070@
#define ARMA_VERSION_GE_070 1
//#define ARMA_VERSION_GE_090 @ARMA_VERSION_GE_090@
#define ARMA_VERSION_GE_090 1

#if ARMA_VERSION_GE_090
#define SCALAR(X) ::arma::as_scalar(X)
#else
#define SCALAR(X) X
#endif

#endif
