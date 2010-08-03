// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)

// RcppArmadillo: add the R defines
#include <Rconfig.h>

#if !defined(ARMA_USE_LAPACK)
#define ARMA_USE_LAPACK
#endif

#if !defined(ARMA_USE_BLAS)
#define ARMA_USE_BLAS
#endif

// #define ARMA_USE_ATLAS
// #define ARMA_ATLAS_INCLUDE_DIR /usr/include/
//// If you're using ATLAS and the compiler can't find cblas.h and/or clapack.h
//// uncomment the above define and specify the appropriate include directory.
//// Make sure the directory has a trailing /

// #define ARMA_USE_BOOST
// #define ARMA_USE_BOOST_DATE

#define ARMA_HAVE_STD_ISFINITE
#define ARMA_HAVE_STD_ISINF
#define ARMA_HAVE_STD_ISNAN
#define ARMA_HAVE_STD_SNPRINTF

// #define ARMA_HAVE_LOG1P
// #define ARMA_HAVE_GETTIMEOFDAY

// #define ARMA_EXTRA_DEBUG
// #define ARMA_NO_DEBUG

// RcppArmadillo: add the RcppArmadillo configuration
#include <RcppArmadilloConfig.h>
