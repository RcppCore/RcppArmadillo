// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup fn_eps
//! @{



//! \brief
//! eps version for non-complex matrices and vectors
template<typename T1>
inline
const eOp<T1, eop_eps>
eps(const Base<typename T1::elem_type, T1>& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  return eOp<T1, eop_eps>(X.get_ref());
  }



//! \brief
//! eps version for complex matrices and vectors
template<typename T1>
inline
Mat< typename T1::pod_type >
eps(const Base< std::complex<typename T1::pod_type>, T1>& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  typedef typename T1::pod_type   T;
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  Mat<T> out(A.n_rows, A.n_cols);
  
         T* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const uword n_elem  = A.n_elem;
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = eop_aux::direct_eps( A_mem[i] );
    }
  
  
  return out;
  }



template<typename eT>
arma_inline
arma_warn_unused
typename arma_integral_only<eT>::result
eps(const eT& x)
  {
  arma_ignore(x);
  
  return eT(0);
  }



template<typename eT>
arma_inline
arma_warn_unused
typename arma_float_only<eT>::result
eps(const eT& x)
  {
  return eop_aux::direct_eps(x);
  }



template<typename T>
arma_inline
arma_warn_unused
typename arma_float_only<T>::result
eps(const std::complex<T>& x)
  {
  return eop_aux::direct_eps(x);
  }



//! @}
