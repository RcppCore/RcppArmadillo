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


//! \addtogroup eglue_core
//! @{



template<typename eglue_type>
template<typename T1, typename T2>
arma_inline
typename T1::elem_type
eglue_core<eglue_type>::get_elem(const eGlue<T1, T2, eglue_type>& x, const u32 i)
  {
  // arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return x.P1[i] + x.P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return x.P1[i] - x.P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return x.P1[i] / x.P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return x.P1[i] * x.P2[i]; }
  else
    {
    arma_stop("eglue_core::get_elem(): unhandled eglue_type");
    return eT(0);
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_inline
typename T1::elem_type
eglue_core<eglue_type>::get_elem(const eGlue<T1, T2, eglue_type>& x, const u32 row, const u32 col)
  {
  // arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { return x.P1.at(row,col) + x.P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { return x.P1.at(row,col) - x.P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { return x.P1.at(row,col) / x.P2.at(row,col); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { return x.P1.at(row,col) * x.P2.at(row,col); }
  else
    {
    arma_stop("eglue_core::get_elem(): unhandled eglue_type");
    return eT(0);
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_inline
void
eglue_core<eglue_type>::apply(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  if( (is_Mat<T1>::value == true) && (is_Mat<T2>::value == true) )
    {
    eglue_core<eglue_type>::apply_unwrap(out, x);
    }
  else
    {
    eglue_core<eglue_type>::apply_proxy(out, x);
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_proxy(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  // eglue_type::apply_proxy() function is not allowed to unwrap things
  // (in order to get the input into a common format).
  // the proxy class is already providing objects with element access
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P1 = x.P1;
  const Proxy<T2>& P2 = x.P2;
  
  out.set_size(P1.n_rows, P1.n_cols);
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = P1[i] * P2[i]; }
  else
    {
    arma_stop("eglue_core::apply_proxy(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_unwrap(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<typename Proxy<T1>::stored_type> tmp1(x.P1.Q);
  const unwrap<typename Proxy<T2>::stored_type> tmp2(x.P2.Q);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const eT* B_mem   = B.memptr();
  const u32 n_elem  = A.n_elem;
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = A_mem[i] + B_mem[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = A_mem[i] - B_mem[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = A_mem[i] / B_mem[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = A_mem[i] * B_mem[i]; }
  else
    {
    arma_stop("eglue_core::apply_unwrap(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_plus(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P1 = x.P1;
  const Proxy<T2>& P2 = x.P2;
  
  arma_assert_same_size(out.n_rows, out.n_cols, P1.n_rows, P1.n_cols, "matrix addition");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] += P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] += P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] += P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] += P1[i] * P2[i]; }
  else
    {
    arma_stop("eglue_core::apply_inplace_plus(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_minus(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P1 = x.P1;
  const Proxy<T2>& P2 = x.P2;
  
  arma_assert_same_size(out.n_rows, out.n_cols, P1.n_rows, P1.n_cols, "matrix subtraction");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] -= P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] -= P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] -= P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] -= P1[i] * P2[i]; }
  else 
    {
    arma_stop("eglue_core::apply_inplace_minus(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_schur(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P1 = x.P1;
  const Proxy<T2>& P2 = x.P2;
  
  arma_assert_same_size(out.n_rows, out.n_cols, P1.n_rows, P1.n_cols, "element-wise matrix multiplication");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] *= P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] *= P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] *= P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] *= P1[i] * P2[i]; }
  else
    {
    arma_stop("eglue_core::apply_inplace_schur(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_div(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1>& P1 = x.P1;
  const Proxy<T2>& P2 = x.P2;
  
  arma_assert_same_size(out.n_rows, out.n_cols, P1.n_rows, P1.n_cols, "element-wise matrix division");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] /= P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] /= P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] /= P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] /= P1[i] * P2[i]; }
  else
    {
    arma_stop("eglue_core::apply_inplace_div(): unhandled eglue_type");
    }
  }



//! @}
