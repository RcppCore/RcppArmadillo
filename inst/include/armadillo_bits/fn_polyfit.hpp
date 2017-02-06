// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_polyfit
//! @{



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  is_supported_blas_type<typename T1::elem_type>::value,
  bool
  >::result
polyfit(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& X, const Base<typename T1::elem_type, T2>& Y, const uword N)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_polyfit::apply_direct(out, X.get_ref(), Y.get_ref(), N);
  
  if(status == false)
    {
    out.reset();
    arma_debug_warn("polyfit(): failed");
    }
  
  return status;
  }



template<typename T1, typename T2>
arma_warn_unused
inline
typename
enable_if2
  <
  is_supported_blas_type<typename T1::elem_type>::value,
  const Glue<T1, T2, glue_polyfit>
  >::result
polyfit(const Base<typename T1::elem_type, T1>& X, const Base<typename T1::elem_type, T2>& Y, const uword N)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_polyfit>(X.get_ref(), Y.get_ref(), N);
  }



//! @}
