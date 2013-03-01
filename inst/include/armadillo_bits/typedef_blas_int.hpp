// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


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
