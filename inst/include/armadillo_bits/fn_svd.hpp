// Copyright (C) 2009-2013 Conrad Sanderson
// Copyright (C) 2009-2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_svd
//! @{



template<typename T1>
inline
bool
svd
  (
         Col<typename T1::pod_type>&     S,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  // it doesn't matter if X is an alias of S, as auxlib::svd() makes a copy of X
  
  const bool status = auxlib::svd(S, X);
  
  if(status == false)
    {
    S.reset();
    arma_bad("svd(): failed to converge", false);
    }
  
  return status;
  }



template<typename T1>
inline
Col<typename T1::pod_type>
svd
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  Col<typename T1::pod_type> out;
  
  const bool status = auxlib::svd(out, X);
  
  if(status == false)
    {
    out.reset();
    arma_bad("svd(): failed to converge");
    }
  
  return out;
  }



template<typename T1>
inline
bool
svd
  (
         Mat<typename T1::elem_type>&    U,
         Col<typename T1::pod_type >&    S,
         Mat<typename T1::elem_type>&    V,
  const Base<typename T1::elem_type,T1>& X,
  const char*                            method = "standard",
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check
    (
    ( ((void*)(&U) == (void*)(&S)) || (&U == &V) || ((void*)(&S) == (void*)(&V)) ),
    "svd(): two or more output objects are the same object"
    );
  
  const char sig = method[0];
  
  arma_debug_check( ((sig != 's') && (sig != 'd')), "svd(): unknown method specified" );
  
  // auxlib::svd() makes an internal copy of X
#if ARMA_CAN_USE_ZGESDD  
  const bool status = (sig == 'd') ? auxlib::svd_dc(U, S, V, X) : auxlib::svd(U, S, V, X);
#else
  const bool status = auxlib::svd(U, S, V, X);
#endif

  if(status == false)
    {
    U.reset();
    S.reset();
    V.reset();
    arma_bad("svd(): failed to converge", false);
    }
  
  return status;
  }



template<typename T1>
inline
bool
svd_econ
  (
         Mat<typename T1::elem_type>&    U,
         Col<typename T1::pod_type >&    S,
         Mat<typename T1::elem_type>&    V,
  const Base<typename T1::elem_type,T1>& X,
  const char                             mode,
  const char*                            method = "standard",
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check
    (
    ( ((void*)(&U) == (void*)(&S)) || (&U == &V) || ((void*)(&S) == (void*)(&V)) ),
    "svd_econ(): two or more output objects are the same object"
    );
  
  arma_debug_check
    (
    ( (mode != 'l') && (mode != 'r') && (mode != 'b') ),
    "svd_econ(): parameter 'mode' or 'side' is incorrect"
    );
  
  const char sig = method[0];
  
  arma_debug_check( ((sig != 's') && (sig != 'd')), "svd_econ(): unknown method specified" );

#if ARMA_CAN_USE_ZGESDD  
  const bool status = ((mode == 'b') && (sig == 'd')) ? auxlib::svd_dc_econ(U, S, V, X) : auxlib::svd_econ(U, S, V, X, mode);
#else
  const bool status = auxlib::svd_econ(U, S, V, X, mode);
#endif

  if(status == false)
    {
    U.reset();
    S.reset();
    V.reset();
    arma_bad("svd_econ(): failed to converge", false);
    }
  
  return status;
  }



template<typename T1>
inline
bool
svd_econ
  (
         Mat<typename T1::elem_type>&    U,
         Col<typename T1::pod_type >&    S,
         Mat<typename T1::elem_type>&    V,
  const Base<typename T1::elem_type,T1>& X,
  const char*                            side   = "both",
  const char*                            method = "standard",
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return svd_econ(U,S,V,X,side[0],method);
  }



// TODO: remove this function in version 4.0
template<typename T1>
arma_deprecated
inline
bool
svd_thin
  (
         Mat<typename T1::elem_type>&    U,
         Col<typename T1::pod_type >&    S,
         Mat<typename T1::elem_type>&    V,
  const Base<typename T1::elem_type,T1>& X,
  const char                             mode = 'b',
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_ignore(junk);
  
  return svd_econ(U,S,V,X,mode);
  }



//! @}
