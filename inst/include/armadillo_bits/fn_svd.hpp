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


//! \addtogroup fn_svd
//! @{



template<typename T1>
inline
bool
svd
  (
  Col<typename T1::pod_type>& S,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  // unwrap_check not used as T1::elem_type and T1::pod_type may not be the same.
  // furthermore, it doesn't matter if A is an alias of S, as auxlib::svd() makes a copy of A
  
  const unwrap<T1> tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  const bool status = auxlib::svd(S, A);
    
  if(status == false)
    {
    arma_print("svd(): singular value decomposition failed");
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
  
  Col<typename T1::pod_type> out;
  
  const bool status = svd(out, X);
  
  if(status == false)
    {
    out.set_size(0);
    }
  
  return out;
  }



template<typename T1>
inline
bool
svd
  (
         Mat<typename T1::elem_type>&    U,
         Col<typename T1::pod_type>&     S,
         Mat<typename T1::elem_type>&    V,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( ( ((void*)(&U) == (void*)(&S)) || (&U == &V) || ((void*)(&S) == (void*)(&V)) ), "svd(): two or more output objects are the same object" );
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  // auxlib::svd() makes an internal copy of A
  const bool status = auxlib::svd(U, S, V, A);
  
  if(status == false)
    {
    arma_print("svd(): singular value decomposition failed");
    }
  
  return status;
  }



//! @}
