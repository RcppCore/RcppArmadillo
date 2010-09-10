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


//! \addtogroup eglue_cube_core
//! @{



template<typename eglue_type>
template<typename T1, typename T2>
arma_inline
typename T1::elem_type
eglue_cube_core<eglue_type>::get_elem(const eGlueCube<T1, T2, eglue_type>& x, const u32 i)
  {
  typedef typename T1::elem_type eT;
  
  // the optimiser will keep only one return statement
  
       if(is_same_type<eglue_type, eglue_cube_plus >::value == true) { return x.P1[i] + x.P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_minus>::value == true) { return x.P1[i] - x.P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_div  >::value == true) { return x.P1[i] / x.P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_schur>::value == true) { return x.P1[i] * x.P2[i]; }
  else
    {
    arma_stop("eglue_cube_core::get_elem(): unhandled eglue_type");
    return eT(0);
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_inline
typename T1::elem_type
eglue_cube_core<eglue_type>::get_elem(const eGlueCube<T1, T2, eglue_type>& x, const u32 row, const u32 col, const u32 slice)
  {
  typedef typename T1::elem_type eT;
  
       if(is_same_type<eglue_type, eglue_cube_plus >::value == true) { return x.P1.at(row,col,slice) + x.P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_cube_minus>::value == true) { return x.P1.at(row,col,slice) - x.P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_cube_div  >::value == true) { return x.P1.at(row,col,slice) / x.P2.at(row,col,slice); }
  else if(is_same_type<eglue_type, eglue_cube_schur>::value == true) { return x.P1.at(row,col,slice) * x.P2.at(row,col,slice); }
  else
    {
    arma_stop("eglue_cube_core::get_elem(): unhandled eglue_type");
    return eT(0);
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_inline
void
eglue_cube_core<eglue_type>::apply(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  if( (is_Cube<T1>::value == true) && (is_Cube<T2>::value == true) )
    {
    eglue_cube_core<eglue_type>::apply_unwrap(out, x);
    }
  else
    {
    eglue_cube_core<eglue_type>::apply_proxy(out, x);
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_cube_core<eglue_type>::apply_proxy(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  // eglue_type::apply_proxy() function is not allowed to unwrap things
  // (in order to get the input into a common format).
  // the proxy class is already providing objects with element access
  
  typedef typename T1::elem_type eT;
  
  const ProxyCube<T1>& P1 = x.P1;
  const ProxyCube<T2>& P2 = x.P2;
  
  out.set_size(P1.n_rows, P1.n_cols, P1.n_slices);
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_cube_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = P1[i] * P2[i]; }
  else
    {
    arma_stop("eglue_cube_core::apply_proxy(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_cube_core<eglue_type>::apply_unwrap(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<typename ProxyCube<T1>::stored_type> tmp1(x.P1.Q);
  const unwrap_cube<typename ProxyCube<T2>::stored_type> tmp2(x.P2.Q);
  
  const Cube<eT>& A = tmp1.M;
  const Cube<eT>& B = tmp2.M;
  
  out.set_size(A.n_rows, A.n_cols, A.n_slices);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const eT* B_mem   = B.memptr();
  const u32 n_elem  = A.n_elem;
  
       if(is_same_type<eglue_type, eglue_cube_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = A_mem[i] + B_mem[i]; }
  else if(is_same_type<eglue_type, eglue_cube_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = A_mem[i] - B_mem[i]; }
  else if(is_same_type<eglue_type, eglue_cube_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = A_mem[i] / B_mem[i]; }
  else if(is_same_type<eglue_type, eglue_cube_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] = A_mem[i] * B_mem[i]; }
  else
    {
    arma_stop("eglue_cube_core::apply_unwrap(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_cube_core<eglue_type>::apply_inplace_plus(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const ProxyCube<T1>& P1 = x.P1;
  const ProxyCube<T2>& P2 = x.P2;
  
  arma_assert_same_size(out.n_rows, out.n_cols, out.n_slices, P1.n_rows, P1.n_cols, P1.n_slices, "cube addition");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_cube_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] += P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] += P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] += P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] += P1[i] * P2[i]; }
  else
    {
    arma_stop("eglue_cube_core::apply_inplace_plus(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_cube_core<eglue_type>::apply_inplace_minus(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const ProxyCube<T1>& P1 = x.P1;
  const ProxyCube<T2>& P2 = x.P2;
  
  arma_assert_same_size(out.n_rows, out.n_cols, out.n_slices, P1.n_rows, P1.n_cols, P1.n_slices, "cube subtraction");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_cube_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] -= P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] -= P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] -= P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] -= P1[i] * P2[i]; }
  else 
    {
    arma_stop("eglue_cube_core::apply_inplace_minus(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_cube_core<eglue_type>::apply_inplace_schur(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const ProxyCube<T1>& P1 = x.P1;
  const ProxyCube<T2>& P2 = x.P2;
  
  arma_assert_same_size(out.n_rows, out.n_cols, out.n_slices, P1.n_rows, P1.n_cols, P1.n_slices, "element-wise cube multiplication");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_cube_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] *= P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] *= P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] *= P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] *= P1[i] * P2[i]; }
  else
    {
    arma_stop("eglue_cube_core::apply_inplace_schur(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_cube_core<eglue_type>::apply_inplace_div(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const ProxyCube<T1>& P1 = x.P1;
  const ProxyCube<T2>& P2 = x.P2;
  
  arma_assert_same_size(out.n_rows, out.n_cols, out.n_slices, P1.n_rows, P1.n_cols, P1.n_slices, "element-wise cube division");
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = P1.n_elem;
  
       if(is_same_type<eglue_type, eglue_cube_plus >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] /= P1[i] + P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_minus>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] /= P1[i] - P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_div  >::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] /= P1[i] / P2[i]; }
  else if(is_same_type<eglue_type, eglue_cube_schur>::value == true) for(u32 i=0; i<n_elem; ++i) { out_mem[i] /= P1[i] * P2[i]; }
  else
    {
    arma_stop("eglue_cube_core::apply_inplace_div(): unhandled eglue_type");
    }
  }



//! @}
