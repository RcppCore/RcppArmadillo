// Copyright (C) 2010-2015 Conrad Sanderson
// Copyright (C) 2010-2015 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_conv
//! @{



//! Convolution, which is also equivalent to polynomial multiplication and FIR digital filtering.

template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const Glue<T1, T2, glue_conv>
  >::result
conv(const T1& A, const T2& B)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_conv>(A, B);
  }



//! @}
