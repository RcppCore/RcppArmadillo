// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


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
  const char* method =                   "",
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
  
  bool use_divide_and_conquer = false;
  
  const char sig = method[0];
  
  switch(sig)
    {
    case '\0':
    case 's':
      break;
      
    case 'd':
      use_divide_and_conquer = true;
      break;
    
    default:
      {
      arma_stop("svd(): unknown method specified");
      return false;
      }
    }
  
  // auxlib::svd() makes an internal copy of X
  const bool status = (use_divide_and_conquer == false) ? auxlib::svd(U, S, V, X) : auxlib::svd_dc(U, S, V, X);
  
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
  const char                             mode = 'b',
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
    "svd_econ(): parameter 'mode' is incorrect"
    );
  
  
  // auxlib::svd_econ() makes an internal copy of X
  const bool status = auxlib::svd_econ(U, S, V, X, mode);
  
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
