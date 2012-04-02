// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup glue_times
//! @{



template<bool is_eT_blas_type>
template<typename T1, typename T2>
arma_hot
inline
void
glue_times_redirect2_helper<is_eT_blas_type>::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap_check<T1> tmp1(X.A, out);
  const partial_unwrap_check<T2> tmp2(X.B, out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  const bool use_alpha = partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times;
  const eT       alpha = use_alpha ? (tmp1.get_val() * tmp2.get_val()) : eT(0);
  
  glue_times::apply
    <
    eT,
    partial_unwrap_check<T1>::do_trans,
    partial_unwrap_check<T2>::do_trans,
    (partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times)
    >
    (out, A, B, alpha);
  }



template<typename T1, typename T2>
arma_hot
inline
void
glue_times_redirect2_helper<true>::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(strip_inv<T1>::do_inv == false)
    {
    const partial_unwrap_check<T1> tmp1(X.A, out);
    const partial_unwrap_check<T2> tmp2(X.B, out);
    
    const Mat<eT>& A = tmp1.M;
    const Mat<eT>& B = tmp2.M;
    
    const bool use_alpha = partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times;
    const eT       alpha = use_alpha ? (tmp1.get_val() * tmp2.get_val()) : eT(0);
    
    glue_times::apply
      <
      eT,
      partial_unwrap_check<T1>::do_trans,
      partial_unwrap_check<T2>::do_trans,
      (partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times)
      >
      (out, A, B, alpha);
    }
  else
    {
    arma_extra_debug_print("glue_times_redirect<2>::apply(): detected inv(A)*B");
    
    typedef typename strip_inv<T1>::stored_type T1_stripped;
    
    const strip_inv<T1> A_strip(X.A);
    
    Mat<eT> A = A_strip.M;
    
    arma_debug_check( (A.is_square() == false), "inv(): given matrix is not square" );
    
    const unwrap_check<T2> B_tmp(X.B, out);
    const Mat<eT>& B = B_tmp.M;
    
    glue_solve::solve_direct( out, A, B, A_strip.slow );
    }
  }



template<uword N>
template<typename T1, typename T2>
arma_hot
inline
void
glue_times_redirect<N>::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap_check<T1> tmp1(X.A, out);
  const partial_unwrap_check<T2> tmp2(X.B, out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  const bool use_alpha = partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times;
  const eT       alpha = use_alpha ? (tmp1.get_val() * tmp2.get_val()) : eT(0);
  
  glue_times::apply
    <
    eT,
    partial_unwrap_check<T1>::do_trans,
    partial_unwrap_check<T2>::do_trans,
    (partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times)
    >
    (out, A, B, alpha);
  }



template<typename T1, typename T2>
arma_hot
inline
void
glue_times_redirect<2>::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  glue_times_redirect2_helper< is_supported_blas_type<eT>::value >::apply(out, X);
  }



template<typename T1, typename T2, typename T3>
arma_hot
inline
void
glue_times_redirect<3>::apply(Mat<typename T1::elem_type>& out, const Glue< Glue<T1,T2,glue_times>, T3, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // there is exactly 3 objects
  // hence we can safely expand X as X.A.A, X.A.B and X.B
  
  const partial_unwrap_check<T1> tmp1(X.A.A, out);
  const partial_unwrap_check<T2> tmp2(X.A.B, out);
  const partial_unwrap_check<T3> tmp3(X.B,   out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  const Mat<eT>& C = tmp3.M;
  
  const bool use_alpha = partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times || partial_unwrap_check<T3>::do_times;
  const eT       alpha = use_alpha ? (tmp1.get_val() * tmp2.get_val() * tmp3.get_val()) : eT(0);
  
  glue_times::apply
    <
    eT,
    partial_unwrap_check<T1>::do_trans,
    partial_unwrap_check<T2>::do_trans,
    partial_unwrap_check<T3>::do_trans,
    (partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times || partial_unwrap_check<T3>::do_times)
    >
    (out, A, B, C, alpha);
  }



template<typename T1, typename T2, typename T3, typename T4>
arma_hot
inline
void
glue_times_redirect<4>::apply(Mat<typename T1::elem_type>& out, const Glue< Glue< Glue<T1,T2,glue_times>, T3, glue_times>, T4, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // there is exactly 4 objects
  // hence we can safely expand X as X.A.A.A, X.A.A.B, X.A.B and X.B
  
  const partial_unwrap_check<T1> tmp1(X.A.A.A, out);
  const partial_unwrap_check<T2> tmp2(X.A.A.B, out);
  const partial_unwrap_check<T3> tmp3(X.A.B,   out);
  const partial_unwrap_check<T4> tmp4(X.B,     out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  const Mat<eT>& C = tmp3.M;
  const Mat<eT>& D = tmp4.M;
  
  const bool use_alpha = partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times || partial_unwrap_check<T3>::do_times || partial_unwrap_check<T4>::do_times;
  const eT       alpha = use_alpha ? (tmp1.get_val() * tmp2.get_val() * tmp3.get_val() * tmp4.get_val()) : eT(0);
  
  glue_times::apply
    <
    eT,
    partial_unwrap_check<T1>::do_trans,
    partial_unwrap_check<T2>::do_trans,
    partial_unwrap_check<T3>::do_trans,
    partial_unwrap_check<T4>::do_trans,
    (partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times || partial_unwrap_check<T3>::do_times || partial_unwrap_check<T4>::do_times)
    >
    (out, A, B, C, D, alpha);
  }



template<typename T1, typename T2>
arma_hot
inline
void
glue_times::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const sword N_mat = 1 + depth_lhs< glue_times, Glue<T1,T2,glue_times> >::num;

  arma_extra_debug_print(arma_boost::format("N_mat = %d") % N_mat);

  glue_times_redirect<N_mat>::apply(out, X);
  }



template<typename T1>
arma_hot
inline
void
glue_times::apply_inplace(Mat<typename T1::elem_type>& out, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(X, out);
  const Mat<eT>& B     = tmp.M;
  
  arma_debug_assert_mul_size(out, B, "matrix multiplication");
  
  const uword out_n_rows = out.n_rows;
  const uword out_n_cols = out.n_cols;
  
  if(out_n_cols == B.n_cols)
    {
    // size of resulting matrix is the same as 'out'
    
    podarray<eT> tmp(out_n_cols);
    
    eT* tmp_rowdata = tmp.memptr();
    
    for(uword row=0; row < out_n_rows; ++row)
      {
      tmp.copy_row(out, row);
      
      for(uword col=0; col < out_n_cols; ++col)
        {
        out.at(row,col) = op_dot::direct_dot( out_n_cols, tmp_rowdata, B.colptr(col) );
        }
      }
    
    }
  else
    {
    const Mat<eT> tmp(out);
    
    glue_times::apply<eT, false, false, false>(out, tmp, B, eT(1));
    }
  
  }



template<typename T1, typename T2>
arma_hot
inline
void
glue_times::apply_inplace_plus(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_times>& X, const sword sign)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap_check<T1> tmp1(X.A, out);
  const partial_unwrap_check<T2> tmp2(X.B, out);
  
  const Mat<eT>& A     = tmp1.M;
  const Mat<eT>& B     = tmp2.M;
  
  const bool do_trans_A = partial_unwrap_check<T1>::do_trans;
  const bool do_trans_B = partial_unwrap_check<T2>::do_trans;
  
  const bool use_alpha = partial_unwrap_check<T1>::do_times || partial_unwrap_check<T2>::do_times || (sign < sword(0));
  const eT       alpha = use_alpha ? ( tmp1.get_val() * tmp2.get_val() * ( (sign > sword(0)) ? eT(1) : eT(-1) ) ) : eT(0);
  
  arma_debug_assert_mul_size(A, B, do_trans_A, do_trans_B, "matrix multiplication");
  
  const uword result_n_rows = (do_trans_A == false) ? A.n_rows : A.n_cols;
  const uword result_n_cols = (do_trans_B == false) ? B.n_cols : B.n_rows;
  
  arma_assert_same_size(out.n_rows, out.n_cols, result_n_rows, result_n_cols, "addition");
  
  if(out.n_elem > 0)
    {
    if( (do_trans_A == false) && (do_trans_B == false) && (use_alpha == false) )
      {
      if( (A.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,         false, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1));
        }
      else
      if(B.n_cols == 1)
        {
        gemv<false,        false, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1));
        }
      else
        {
        gemm<false, false, false, true>::apply(out, A, B, alpha, eT(1));
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == false) && (use_alpha == true) )
      {
      if( (A.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,         true, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1));
        }
      else
      if(B.n_cols == 1)
        {
        gemv<false,        true, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1));
        }
      else
        {
        gemm<false, false, true, true>::apply(out, A, B, alpha, eT(1));
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == false) && (use_alpha == false) )
      {
      if( (A.n_cols == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,        false, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1));
        }
      else
      if(B.n_cols == 1)
        {
        gemv<true,        false, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1));
        }
      else
        {
        gemm<true, false, false, true>::apply(out, A, B, alpha, eT(1));
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == false) && (use_alpha == true) )
      {
      if( (A.n_cols == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,        true, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1));
        }
      else
      if(B.n_cols == 1)
        {
        gemv<true,        true, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1));
        }
      else
        {
        gemm<true, false, true, true>::apply(out, A, B, alpha, eT(1));
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == true) && (use_alpha == false) )
      {
      if( (A.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,       false, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1));
        }
      else
      if( (B.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,       false, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1));
        }
      else
        {
        gemm<false, true, false, true>::apply(out, A, B, alpha, eT(1));
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == true) && (use_alpha == true) )
      {
      if( (A.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,       true, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1));
        }
      else
      if( (B.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,       true, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1));
        }
      else
        {
        gemm<false, true, true, true>::apply(out, A, B, alpha, eT(1));
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == true) && (use_alpha == false) )
      {
      if( (A.n_cols == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,      false, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1));
        }
      else
      if( (B.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,       false, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1));
        }
      else
        {
        gemm<true, true, false, true>::apply(out, A, B, alpha, eT(1));
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == true) && (use_alpha == true) )
      {
      if( (A.n_cols == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,      true, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1));
        }
      else
      if( (B.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,       true, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1));
        }
      else
        {
        gemm<true, true, true, true>::apply(out, A, B, alpha, eT(1));
        }
      }
    }
  
  
  }



template<typename eT, const bool do_trans_A, const bool do_trans_B>
arma_inline
uword
glue_times::mul_storage_cost(const Mat<eT>& A, const Mat<eT>& B)
  {
  const uword final_A_n_rows = (do_trans_A == false) ? A.n_rows : A.n_cols;
  const uword final_B_n_cols = (do_trans_B == false) ? B.n_cols : B.n_rows;
  
  return final_A_n_rows * final_B_n_cols;
  }



template
  <
  typename   eT,
  const bool do_trans_A,
  const bool do_trans_B,
  const bool use_alpha
  >
arma_hot
inline
void
glue_times::apply
  (
        Mat<eT>& out,
  const Mat<eT>& A,
  const Mat<eT>& B,
  const eT       alpha
  )
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_mul_size(A, B, do_trans_A, do_trans_B, "matrix multiplication");
  
  const uword final_n_rows = (do_trans_A == false) ? A.n_rows : A.n_cols;
  const uword final_n_cols = (do_trans_B == false) ? B.n_cols : B.n_rows;
  
  out.set_size(final_n_rows, final_n_cols);
  
  if( (A.n_elem > 0) && (B.n_elem > 0) )
    {
    if( (do_trans_A == false) && (do_trans_B == false) && (use_alpha == false) )
      {
      if( (A.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,         false, false>::apply(out.memptr(), B, A.memptr());
        }
      else
      if(B.n_cols == 1)
        {
        gemv<false,        false, false>::apply(out.memptr(), A, B.memptr());
        }
      else
        {
        gemm<false, false, false, false>::apply(out, A, B);
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == false) && (use_alpha == true) )
      {
      if( (A.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,         true, false>::apply(out.memptr(), B, A.memptr(), alpha);
        }
      else
      if(B.n_cols == 1)
        {
        gemv<false,        true, false>::apply(out.memptr(), A, B.memptr(), alpha);
        }
      else
        {
        gemm<false, false, true, false>::apply(out, A, B, alpha);
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == false) && (use_alpha == false) )
      {
      if( (A.n_cols == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,        false, false>::apply(out.memptr(), B, A.memptr());
        }
      else
      if(B.n_cols == 1)
        {
        gemv<true,        false, false>::apply(out.memptr(), A, B.memptr());
        }
      else
        {
        gemm<true, false, false, false>::apply(out, A, B);
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == false) && (use_alpha == true) )
      {
      if( (A.n_cols == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,        true, false>::apply(out.memptr(), B, A.memptr(), alpha);
        }
      else
      if(B.n_cols == 1)
        {
        gemv<true,        true, false>::apply(out.memptr(), A, B.memptr(), alpha);
        }
      else
        {
        gemm<true, false, true, false>::apply(out, A, B, alpha);
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == true) && (use_alpha == false) )
      {
      if( (A.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,       false, false>::apply(out.memptr(), B, A.memptr());
        }
      else
      if( (B.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,       false, false>::apply(out.memptr(), A, B.memptr());
        }
      else
        {
        gemm<false, true, false, false>::apply(out, A, B);
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == true) && (use_alpha == true) )
      {
      if( (A.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,       true, false>::apply(out.memptr(), B, A.memptr(), alpha);
        }
      else
      if( (B.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,       true, false>::apply(out.memptr(), A, B.memptr(), alpha);
        }
      else
        {
        gemm<false, true, true, false>::apply(out, A, B, alpha);
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == true) && (use_alpha == false) )
      {
      if( (A.n_cols == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,      false, false>::apply(out.memptr(), B, A.memptr());
        }
      else
      if( (B.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,       false, false>::apply(out.memptr(), A, B.memptr());
        }
      else
        {
        gemm<true, true, false, false>::apply(out, A, B);
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == true) && (use_alpha == true) )
      {
      if( (A.n_cols == 1) && (is_complex<eT>::value == false) )
        {
        gemv<false,      true, false>::apply(out.memptr(), B, A.memptr(), alpha);
        }
      else
      if( (B.n_rows == 1) && (is_complex<eT>::value == false) )
        {
        gemv<true,       true, false>::apply(out.memptr(), A, B.memptr(), alpha);
        }
      else
        {
        gemm<true, true, true, false>::apply(out, A, B, alpha);
        }
      }
    }
  else
    {
    out.zeros();
    }
  }



template
  <
  typename   eT,
  const bool do_trans_A,
  const bool do_trans_B,
  const bool do_trans_C,
  const bool use_alpha
  >
arma_hot
inline
void
glue_times::apply
  (
        Mat<eT>& out,
  const Mat<eT>& A,
  const Mat<eT>& B,
  const Mat<eT>& C,
  const eT       alpha
  )
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  
  if( glue_times::mul_storage_cost<eT, do_trans_A, do_trans_B>(A, B) <= glue_times::mul_storage_cost<eT, do_trans_B, do_trans_C>(B, C) )
    {
    // out = (A*B)*C
    
    glue_times::apply<eT, do_trans_A, do_trans_B, use_alpha>(tmp, A,   B, alpha);
    glue_times::apply<eT, false,      do_trans_C, false    >(out, tmp, C, eT(0));
    }
  else
    {
    // out = A*(B*C)
    
    glue_times::apply<eT, do_trans_B, do_trans_C, use_alpha>(tmp, B, C,   alpha);
    glue_times::apply<eT, do_trans_A, false,      false    >(out, A, tmp, eT(0));
    }
  }



template
  <
  typename   eT,
  const bool do_trans_A,
  const bool do_trans_B,
  const bool do_trans_C,
  const bool do_trans_D,
  const bool use_alpha
  >
arma_hot
inline
void
glue_times::apply
  (
        Mat<eT>& out,
  const Mat<eT>& A,
  const Mat<eT>& B,
  const Mat<eT>& C,
  const Mat<eT>& D,
  const eT       alpha
  )
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp;
  
  if( glue_times::mul_storage_cost<eT, do_trans_A, do_trans_C>(A, C) <= glue_times::mul_storage_cost<eT, do_trans_B, do_trans_D>(B, D) )
    {
    // out = (A*B*C)*D
    
    glue_times::apply<eT, do_trans_A, do_trans_B, do_trans_C, use_alpha>(tmp, A, B, C, alpha);
    
    glue_times::apply<eT, false, do_trans_D, false>(out, tmp, D, eT(0));
    }
  else
    {
    // out = A*(B*C*D)
    
    glue_times::apply<eT, do_trans_B, do_trans_C, do_trans_D, use_alpha>(tmp, B, C, D, alpha);
    
    glue_times::apply<eT, do_trans_A, false, false>(out, A, tmp, eT(0));
    }
  }



//
// glue_times_diag


template<typename T1, typename T2>
arma_hot
inline
void
glue_times_diag::apply(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_times_diag>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const strip_diagmat<T1> S1(X.A);
  const strip_diagmat<T2> S2(X.B);
  
  typedef typename strip_diagmat<T1>::stored_type T1_stripped;
  typedef typename strip_diagmat<T2>::stored_type T2_stripped;
  
  if( (strip_diagmat<T1>::do_diagmat == true) && (strip_diagmat<T2>::do_diagmat == false) )
    {
    const diagmat_proxy_check<T1_stripped> A(S1.M, out);
    
    const unwrap_check<T2> tmp(X.B, out);
    const Mat<eT>& B     = tmp.M;
    
    const uword A_n_elem = A.n_elem;
    const uword B_n_rows = B.n_rows;
    const uword B_n_cols = B.n_cols;
    
    arma_debug_assert_mul_size(A_n_elem, A_n_elem, B_n_rows, B_n_cols, "matrix multiplication");
    
    out.set_size(A_n_elem, B_n_cols);
    
    for(uword col=0; col < B_n_cols; ++col)
      {
            eT* out_coldata = out.colptr(col);
      const eT* B_coldata   = B.colptr(col);
      
      uword i,j;
      for(i=0, j=1; j < B_n_rows; i+=2, j+=2)
        {
        eT tmp_i = A[i];
        eT tmp_j = A[j];
        
        tmp_i *= B_coldata[i];
        tmp_j *= B_coldata[j];
        
        out_coldata[i] = tmp_i;
        out_coldata[j] = tmp_j;
        }
      
      if(i < B_n_rows)
        {
        out_coldata[i] = A[i] * B_coldata[i];
        }
      }
    }
  else
  if( (strip_diagmat<T1>::do_diagmat == false) && (strip_diagmat<T2>::do_diagmat == true) )
    {
    const unwrap_check<T1> tmp(X.A, out);
    const Mat<eT>& A     = tmp.M;
    
    const diagmat_proxy_check<T2_stripped> B(S2.M, out);
    
    const uword A_n_rows = A.n_rows;
    const uword A_n_cols = A.n_cols;
    const uword B_n_elem = B.n_elem;
    
    arma_debug_assert_mul_size(A_n_rows, A_n_cols, B_n_elem, B_n_elem, "matrix multiplication");
    
    out.set_size(A_n_rows, B_n_elem);
    
    for(uword col=0; col < A_n_cols; ++col)
      {
      const eT  val = B[col];
      
            eT* out_coldata = out.colptr(col);
      const eT*   A_coldata =   A.colptr(col);
      
      uword i,j;
      for(i=0, j=1; j < A_n_rows; i+=2, j+=2)
        {
        const eT tmp_i = A_coldata[i] * val;
        const eT tmp_j = A_coldata[j] * val;
        
        out_coldata[i] = tmp_i;
        out_coldata[j] = tmp_j;
        }
      
      if(i < A_n_rows)
        {
        out_coldata[i] = A_coldata[i] * val;
        }
      }
    }
  else
  if( (strip_diagmat<T1>::do_diagmat == true) && (strip_diagmat<T2>::do_diagmat == true) )
    {
    const diagmat_proxy_check<T1_stripped> A(S1.M, out);
    const diagmat_proxy_check<T2_stripped> B(S2.M, out);
    
    const uword A_n_elem = A.n_elem;
    const uword B_n_elem = B.n_elem;
    
    arma_debug_assert_mul_size(A_n_elem, A_n_elem, B_n_elem, B_n_elem, "matrix multiplication");
    
    out.zeros(A_n_elem, A_n_elem);
    
    for(uword i=0; i < A_n_elem; ++i)
      {
      out.at(i,i) = A[i] * B[i];
      }
    }
  }



//! @}
