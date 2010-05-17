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


//! \addtogroup operator_div
//! @{



//! Base / scalar
template<typename T1>
arma_inline
const eOp<T1, eop_scalar_div_post>
operator/
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename T1::elem_type           k
  )
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_div_post>(X.get_ref(), k);
  }



//! scalar / Base
template<typename T1>
arma_inline
const eOp<T1, eop_scalar_div_pre>
operator/
  (
  const typename T1::elem_type           k,
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_div_pre>(X.get_ref(), k);
  }



//! non-complex Base / complex scalar (experimental)
template<typename T1>
arma_inline
Mat<typename std::complex<typename T1::pod_type> >
operator/
  (
  const Base<typename T1::pod_type, T1>&     X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.get_ref());
  
  Mat<eT> out(A.n_rows, A.n_cols);
  
  const u32 n_elem  = A.n_elem;
        eT* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = A[i] / k;
    }
  
  return out;
  }



//! complex scalar / non-complex Base (experimental)
template<typename T1>
arma_inline
Mat<typename std::complex<typename T1::pod_type> >
operator/
  (
  const std::complex<typename T1::pod_type>& k,
  const Base<typename T1::pod_type, T1>&     X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<typename T1::pod_type> eT;
  typedef typename T1::pod_type                         T;
  
  const Proxy<T1> A(X.get_ref());
  
  Mat<eT> out(A.n_rows, A.n_cols);
  
  const u32 n_elem  = A.n_elem;
        eT* out_mem = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = k / A[i];
    }
  
  return out;
  }



//! element-wise division of Base objects with same element type
template<typename T1, typename T2>
arma_inline
const eGlue<T1, T2, eglue_div>
operator/
  (
  const Base<typename T1::elem_type,T1>& X,
  const Base<typename T1::elem_type,T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return eGlue<T1, T2, eglue_div>(X.get_ref(), Y.get_ref());
  }



//! element-wise division of Base objects with different element types
template<typename T1, typename T2>
arma_inline
Mat<typename promote_type<typename T1::elem_type, typename T2::elem_type>::result>
operator/
  (
  const Base< typename force_different_type<typename T1::elem_type, typename T2::elem_type>::T1_result, T1>& X,
  const Base< typename force_different_type<typename T1::elem_type, typename T2::elem_type>::T2_result, T2>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
  
  arma_debug_assert_same_size(A, B, "element-wise matrix division");
  
  Mat<out_eT> out(A.n_rows, A.n_cols);

        out_eT* out_mem = out.memptr();
  const u32     n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(A[i]) / upgrade_val<eT1,eT2>::apply(B[i]);
    }
  
  return out;
  }



//! @}
