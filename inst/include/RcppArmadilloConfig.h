// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/* :tabSize=4:indentSize=4:noTabs=false:folding=explicit:collapseFolds=1: */
//
// RcppArmadilloConfig.h: Rcpp/Armadillo glue
//
// Copyright (C)  2010 - 2013  Dirk Eddelbuettel, Romain Francois and Douglas Bates
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

#if !defined(ARMA_USE_LAPACK)
#define ARMA_USE_LAPACK
#endif

#if !defined(ARMA_USE_BLAS)
#define ARMA_USE_BLAS
#endif

#define ARMA_HAVE_STD_ISFINITE
#define ARMA_HAVE_STD_ISINF
#define ARMA_HAVE_STD_ISNAN
#define ARMA_HAVE_STD_SNPRINTF


/* TODO: we might need to undef this on other platforms as well */
#if defined(__GNUC__) && defined(_WIN64) || defined(__FreeBSD__)
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

// Let's be careful for now and undef this as not all compilers support this
//#if defined(ARMA_USE_CXX11)
//#undef ARMA_USE_CXX11
//#endif

// If C++11 has been selected at the R package level, use it for Armadillo too
// This is actually not needed, if the proper switch is set via -std=... then
// Armadillo will know (cf compilation with -DARMA_EXTRA_DEBUG set)
// #if defined(USE_CXX1X)
// #define ARMA_USE_CXX11
// #endif

// Rcpp has its own stream object which cooperates more nicely with R's i/o
// And as of Armadillo 2.4.3, we can use this stream object as well 
#if !defined(ARMA_DEFAULT_OSTREAM)
#define ARMA_DEFAULT_OSTREAM Rcpp::Rcout
#endif

// R now defines NDEBUG which suppresses a number of useful Armadillo tests
// Users can still defined it later, and/or define ARMA_NO_DEBUG
#if defined(NDEBUG)
#undef NDEBUG
#endif

// R can be built with its own Rlapack library, or use an external
// one. Only the latter has zgesdd, a complex-valued SVD using divide-and-conquer 
#if defined(WIN32) || defined(_WIN32)
  // on Windows we do not assume ZGESDD
  #define ARMA_DONT_USE_CX_GESDD 1
#else
  // on the other OSs we test via LAPACK_LIBS (in configure) which
  // updates this include file
  #include <RcppArmadilloLapack.h>
#endif

#endif

