// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup typedef
//! @{


#if   defined(ARMA_BLAS_LONG_LONG)
  typedef long long blas_int;
#elif defined(ARMA_BLAS_LONG)
  typedef long      blas_int;
#else
  typedef int       blas_int;
#endif


//! @}
