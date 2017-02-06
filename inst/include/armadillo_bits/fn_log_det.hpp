// Copyright (C) 2010-2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_log_det
//! @{



//! log determinant of mat
template<typename T1>
inline
void
log_det
  (
        typename T1::elem_type&          out_val,
        typename T1::pod_type&           out_sign,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const bool status = auxlib::log_det(out_val, out_sign, X);
  
  if(status == false)
    {
    out_val  = eT(Datum<T>::nan);
    out_sign = T(0);
    
    arma_warn("log_det(): failed to find determinant");
    }
  }



template<typename T1>
inline
void
log_det
  (
        typename T1::elem_type& out_val,
        typename T1::pod_type&  out_sign,
  const Op<T1,op_diagmat>&      X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const diagmat_proxy<T1> A(X.m);
  
  arma_debug_check( (A.n_rows != A.n_cols), "log_det(): given matrix must be square sized" );
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  if(N == 0)
    {
    out_val  = eT(0);
    out_sign =  T(1);
    
    return;
    }
  
  eT x = A[0];
  
  T  sign = (is_cx<eT>::no) ? ( (access::tmp_real(x) < T(0)) ? -1 : +1 ) : +1;
  eT val  = (is_cx<eT>::no) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x ) : std::log(x);
  
  for(uword i=1; i<N; ++i)
    {
    x = A[i];
    
    sign *= (is_cx<eT>::no) ? ( (access::tmp_real(x) < T(0)) ? -1 : +1 ) : +1;
    val  += (is_cx<eT>::no) ? std::log( (access::tmp_real(x) < T(0)) ? x*T(-1) : x ) : std::log(x);
    }
  
  out_val  = val;
  out_sign = sign;
  }



template<typename T1>
inline
arma_warn_unused
std::complex<typename T1::pod_type>
log_det
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  eT out_val  = eT(0);
   T out_sign =  T(0);
  
  log_det(out_val, out_sign, X.get_ref());
  
  return (is_cx<eT>::yes) ? std::complex<T>(out_val) : ( (out_sign >= T(1)) ? std::complex<T>(access::tmp_real(out_val),T(0)) : std::complex<T>(access::tmp_real(out_val),Datum<T>::pi) );
  }



//! @}
