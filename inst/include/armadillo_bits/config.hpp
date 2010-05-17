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

#include <Rconfig.h>

#define ARMA_HAVE_STD_ISFINITE
#define ARMA_HAVE_STD_ISINF
#define ARMA_HAVE_STD_ISNAN
#define ARMA_HAVE_STD_SNPRINTF

#define ARMA_HAVE_LOG1P
#define ARMA_HAVE_GETTIMEOFDAY

/* #undef ARMA_USE_ATLAS */
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef HAVE_F77_UNDERSCORE
#pragma message (** armadillo calls to BLAS and Lapack disabled **)
#undef ARMA_USE_LAPACK
#undef ARMA_USE_BLAS
#endif
/* #undef ARMA_USE_BOOST */
/* #undef ARMA_USE_BOOST_DATE */

/* #undef ARMA_EXTRA_DEBUG */
/* #undef ARMA_NO_DEBUG */

#if defined(ARMA_USE_ATLAS)
  #if !defined(ARMA_ATLAS_INCLUDE_DIR)
    #define ARMA_ATLAS_INCLUDE_DIR  /usr/include
  #endif
#endif

#include <RcppArmadilloConfig.h>

