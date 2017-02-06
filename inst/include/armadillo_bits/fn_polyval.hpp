// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_polyval
//! @{



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  (is_supported_blas_type<typename T1::elem_type>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type,typename T2::elem_type>::value),
  const Glue<T1, T2, glue_polyval>
  >::result
polyval(const Base<typename T1::elem_type, T1>& P, const T2& X)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_polyval>(P.get_ref(), X);
  }



//! @}
