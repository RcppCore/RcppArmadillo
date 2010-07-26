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


//! \addtogroup fn_log_det
//! @{



//! log determinant of mat
template<typename T1>
inline
void
log_det
  (
  typename T1::elem_type& out_val,
  typename T1::pod_type& out_sign,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( !A.is_square(), "log_det(): matrix must be square" );
  
  auxlib::log_det(out_val, out_sign, A);
  }



template<typename T1>
inline
void
log_det
  (
  typename T1::elem_type& out_val,
  typename T1::pod_type& out_sign,
  const Op<T1,op_diagmat>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const diagmat_proxy<T1> A(X.m);
  
  const u32 N = A.n_elem;
  
  arma_debug_check( (N == 0), "log_det(): given matrix has no elements" );
  
  const eT x = A[0];
  
  T  sign = (is_complex<eT>::value == false) ? ( (access::tmp_real(x) < T(0)) ? -1 : +1 ) : +1;
  eT val  = (is_complex<eT>::value == false) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x ) : std::log(x);

  for(u32 i=1; i<N; ++i)
    {
    const eT x = A[i];
    
    sign *= (is_complex<eT>::value == false) ? ( (access::tmp_real(x) < T(0)) ? -1 : +1 ) : +1;
    val  += (is_complex<eT>::value == false) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x ) : std::log(x);
    }
  
  out_val  = val;
  out_sign = sign;
  }



//! @}
