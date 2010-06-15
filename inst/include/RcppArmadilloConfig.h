// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppArmadilloConfig.h: Rcpp/Armadillo glue
//
// Copyright (C)  2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

#ifndef RcppArmadillo__RcppArmadilloConfig__h
#define RcppArmadillo__RcppArmadilloConfig__h

/* TODO: we might need to undef this on other platforms as well */
#if defined(__GNUC__) && defined(_WIN64)
#undef ARMA_HAVE_STD_SNPRINTF
#endif

/* 
   suncc does not have std::isfinite (which is not standard)
   so we tell armadillo not to use it, and comment out a few 
   others while we are at it
*/
#if defined(__SUNPRO_CC)
#undef ARMA_HAVE_STD_ISFINITE
#undef ARMA_HAVE_STD_SNPRINTF
#undef ARMA_HAVE_LOG1P
#undef ARMA_HAVE_STD_ISINF
#undef ARMA_HAVE_STD_ISNAN
#endif

#endif

