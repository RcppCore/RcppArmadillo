// Copyright (C) 2009-2011 Conrad Sanderson
// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_pinv
//! @{



template<typename T1>
inline
const Op<T1, op_pinv>
pinv
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename T1::elem_type tol = 0.0,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Op<T1, op_pinv>(X.get_ref(), tol);
  }



template<typename T1>
inline
bool
pinv
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& X,
  const typename T1::elem_type tol = 0.0,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  try
    {
    out = pinv(X,tol);
    }
  catch(std::runtime_error&)
    {
    return false;
    }
  
  return true;
  }



//! @}
