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


//! \addtogroup fn_as_scalar
//! @{



template<u32 N>
struct as_scalar_redirect
  {
  template<typename T1>
  inline static typename T1::elem_type apply(const T1& X);
  };



template<>
struct as_scalar_redirect<2>
  {
  template<typename T1, typename T2>
  inline static typename T1::elem_type apply(const Glue<T1,T2,glue_times>& X);
  };


template<>
struct as_scalar_redirect<3>
  {
  template<typename T1, typename T2, typename T3>
  inline static typename T1::elem_type apply(const Glue< Glue<T1, T2, glue_times>, T3, glue_times>& X);
  };



template<u32 N>
template<typename T1>
inline
typename T1::elem_type
as_scalar_redirect<N>::apply(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.n_elem != 1), "as_scalar(): expression doesn't evaluate to exactly one element" );
  
  return A.mem[0];
  }



template<typename T1, typename T2>
inline
typename T1::elem_type
as_scalar_redirect<2>::apply(const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // T1 must result in a matrix with one row
  // T2 must result in a matrix with one column
  
  const partial_unwrap<T1> tmp1(X.A);
  const partial_unwrap<T2> tmp2(X.B);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  const u32 A_n_rows = (tmp1.do_trans == false) ? A.n_rows : A.n_cols;
  const u32 A_n_cols = (tmp1.do_trans == false) ? A.n_cols : A.n_rows;
  
  const u32 B_n_rows = (tmp2.do_trans == false) ? B.n_rows : B.n_cols;
  const u32 B_n_cols = (tmp2.do_trans == false) ? B.n_cols : B.n_rows;
  
  const eT val = tmp1.val * tmp2.val;
  
  arma_debug_check( (A_n_rows != 1) || (B_n_cols != 1) || (A_n_cols != B_n_rows), "as_scalar(): incompatible dimensions" );
  
  return val * op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  }



template<typename T1, typename T2, typename T3>
inline
typename T1::elem_type
as_scalar_redirect<3>::apply(const Glue< Glue<T1, T2, glue_times>, T3, glue_times >& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // T1 * T2 must result in a matrix with one row
  // T3 must result in a matrix with one column
  
  typedef typename strip_inv    <T2           >::stored_type T2_stripped_1;
  typedef typename strip_diagmat<T2_stripped_1>::stored_type T2_stripped_2;
  
  const strip_inv    <T2>            strip1(X.A.B);
  const strip_diagmat<T2_stripped_1> strip2(strip1.M);
  
  const bool tmp2_do_inv     = strip1.do_inv;
  const bool tmp2_do_diagmat = strip2.do_diagmat;
  
  const partial_unwrap<T1>            tmp1(X.A.A);
  const partial_unwrap<T2_stripped_2> tmp2(strip2.M);
  const partial_unwrap<T3>            tmp3(X.B);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  const Mat<eT>& C = tmp3.M;
  
  
  if(tmp2_do_diagmat == false)
    {
    const u32 A_n_rows = (tmp1.do_trans == false) ? A.n_rows : A.n_cols;
    const u32 A_n_cols = (tmp1.do_trans == false) ? A.n_cols : A.n_rows;
    
    const u32 B_n_rows = (tmp2.do_trans == false) ? B.n_rows : B.n_cols;
    const u32 B_n_cols = (tmp2.do_trans == false) ? B.n_cols : B.n_rows;
    
    const u32 C_n_rows = (tmp3.do_trans == false) ? C.n_rows : C.n_cols;
    const u32 C_n_cols = (tmp3.do_trans == false) ? C.n_cols : C.n_rows;
    
    const eT val = tmp1.val * tmp2.val * tmp3.val;
    
    arma_debug_check
      (
      (A_n_rows != 1)        ||
      (C_n_cols != 1)        ||
      (A_n_cols != B_n_rows) ||
      (B_n_cols != C_n_rows)
      ,
      "as_scalar(): incompatible dimensions"
      );
    
    
    if(tmp2_do_inv == true)
      {
      arma_debug_check( (B.is_square() == false), "as_scalar(): incompatible dimensions" );
      
      Mat<eT> B_inv;
      
      if(tmp2.do_trans == false)
        {
        op_inv::apply(B_inv, B);
        }
      else
        {
        const Mat<eT> B_trans = trans(B);
        op_inv::apply(B_inv, B_trans);
        }
      
      return val * op_dotext::direct_rowvec_mat_colvec(A.mem, B_inv, C.mem);
      }
    else
      {
      if(tmp2.do_trans == false)
        {
        return val * op_dotext::direct_rowvec_mat_colvec(A.mem, B, C.mem);
        }
      else
        {
        return val * op_dotext::direct_rowvec_transmat_colvec(A.mem, B, C.mem);
        }
      }
    }
  else
    {
    const u32 A_n_rows = (tmp1.do_trans == false) ? A.n_rows : A.n_cols;
    const u32 A_n_cols = (tmp1.do_trans == false) ? A.n_cols : A.n_rows;
    
    const bool B_is_vec = B.is_vec();
    
    const u32 B_n_rows = (B_is_vec == true) ? B.n_elem : ( (tmp2.do_trans == false) ? B.n_rows : B.n_cols );
    const u32 B_n_cols = (B_is_vec == true) ? B.n_elem : ( (tmp2.do_trans == false) ? B.n_cols : B.n_rows );
    
    const u32 C_n_rows = (tmp3.do_trans == false) ? C.n_rows : C.n_cols;
    const u32 C_n_cols = (tmp3.do_trans == false) ? C.n_cols : C.n_rows;
    
    const eT val = tmp1.val * tmp2.val * tmp3.val;
    
    arma_debug_check
      (
      (A_n_rows != 1)        ||
      (C_n_cols != 1)        ||
      (A_n_cols != B_n_rows) ||
      (B_n_cols != C_n_rows)
      ,
      "as_scalar(): incompatible dimensions"
      );
    
    
    if(B_is_vec == true)
      {
      if(tmp2_do_inv == true)
        {
        return val * op_dotext::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
        }
      else
        {
        return val * op_dot::direct_dot(A.n_elem, A.mem, B.mem, C.mem);
        }
      }
    else
      {
      if(tmp2_do_inv == true)
        {
        return val * op_dotext::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
        }
      else
        {
        return val * op_dotext::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
        }
      }
    }
  }



template<typename T1>
inline
typename T1::elem_type
as_scalar_diag(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.n_elem != 1), "as_scalar(): expression doesn't evaluate to exactly one element" );
  
  return A.mem[0];
  }



template<typename T1, typename T2, typename T3>
inline
typename T1::elem_type
as_scalar_diag(const Glue< Glue<T1, T2, glue_times_diag>, T3, glue_times >& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // T1 * T2 must result in a matrix with one row
  // T3 must result in a matrix with one column
  
  typedef typename strip_diagmat<T2>::stored_type T2_stripped;
  
  const strip_diagmat<T2> strip(X.A.B);
  
  const partial_unwrap<T1>          tmp1(X.A.A);
  const partial_unwrap<T2_stripped> tmp2(strip.M);
  const partial_unwrap<T3>          tmp3(X.B);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  const Mat<eT>& C = tmp3.M;
  
  
  const u32 A_n_rows = (tmp1.do_trans == false) ? A.n_rows : A.n_cols;
  const u32 A_n_cols = (tmp1.do_trans == false) ? A.n_cols : A.n_rows;
  
  const bool B_is_vec = B.is_vec();
  
  const u32 B_n_rows = (B_is_vec == true) ? B.n_elem : ( (tmp2.do_trans == false) ? B.n_rows : B.n_cols );
  const u32 B_n_cols = (B_is_vec == true) ? B.n_elem : ( (tmp2.do_trans == false) ? B.n_cols : B.n_rows );
  
  const u32 C_n_rows = (tmp3.do_trans == false) ? C.n_rows : C.n_cols;
  const u32 C_n_cols = (tmp3.do_trans == false) ? C.n_cols : C.n_rows;
  
  const eT val = tmp1.val * tmp2.val * tmp3.val;
  
  arma_debug_check
    (
    (A_n_rows != 1)        ||
    (C_n_cols != 1)        ||
    (A_n_cols != B_n_rows) ||
    (B_n_cols != C_n_rows)
    ,
    "as_scalar(): incompatible dimensions"
    );
  
  
  if(B_is_vec == true)
    {
    return val * op_dot::direct_dot(A.n_elem, A.mem, B.mem, C.mem);
    }
  else
    {
    return val * op_dotext::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
    }
  }



template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename T1::elem_type
as_scalar(const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  if(is_glue_times_diag<T1>::value == false)
    {
    const s32 N_mat = 1 + depth_lhs< glue_times, Glue<T1,T2,glue_times> >::num;
    
    arma_extra_debug_print(arma_boost::format("N_mat = %d") % N_mat);
    
    return as_scalar_redirect<N_mat>::apply(X);
    }
  else
    {
    return as_scalar_diag(X);
    }
  }



template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
as_scalar(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.n_elem != 1), "as_scalar(): expression doesn't evaluate to exactly one element" );
  
  return A.mem[0];
  }



template<typename T1>
arma_inline
arma_warn_unused
typename T1::elem_type
as_scalar(const eOp<T1, eop_neg>& X)
  {
  arma_extra_debug_sigprint();
  
  return -(as_scalar(X.P.Q));
  }



template<typename T1>
inline
arma_warn_unused
typename T1::elem_type
as_scalar(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> tmp(X.get_ref());
  const Cube<eT>& A   = tmp.M;
  
  arma_debug_check( (A.n_elem != 1), "as_scalar(): expression doesn't evaluate to exactly one element" );
  
  return A.mem[0];
  }



template<typename T>
arma_inline
arma_warn_unused
const typename arma_scalar_only<T>::result &
as_scalar(const T& x)
  {
  return x;
  }



//! @}
