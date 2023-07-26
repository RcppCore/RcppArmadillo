// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup glue_times_misc
//! @{



template<typename eT>
arma_inline
typename arma_not_cx<eT>::result
dense_sparse_helper::dot(const eT* A_mem, const SpMat<eT>& B, const uword col)
  {
  arma_extra_debug_sigprint();
  
        uword      col_offset = B.col_ptrs[col    ];
  const uword next_col_offset = B.col_ptrs[col + 1];
  
  const uword* start_ptr = &(B.row_indices[     col_offset]);
  const uword*   end_ptr = &(B.row_indices[next_col_offset]);
  
  const eT* B_values = B.values;
  
  eT acc = eT(0);
  
  for(const uword* ptr = start_ptr; ptr != end_ptr; ++ptr)
    {
    const uword index = (*ptr);
    
    acc += A_mem[index] * B_values[col_offset];
    
    ++col_offset;
    }
  
  return acc;
  }



template<typename eT>
arma_inline
typename arma_cx_only<eT>::result
dense_sparse_helper::dot(const eT* A_mem, const SpMat<eT>& B, const uword col)
  {
  arma_extra_debug_sigprint();
  
  typedef typename get_pod_type<eT>::result T;
  
        uword      col_offset = B.col_ptrs[col    ];
  const uword next_col_offset = B.col_ptrs[col + 1];
  
  const uword* start_ptr = &(B.row_indices[     col_offset]);
  const uword*   end_ptr = &(B.row_indices[next_col_offset]);
  
  const eT* B_values = B.values;
  
  T acc_real = T(0);
  T acc_imag = T(0);
  
  for(const uword* ptr = start_ptr; ptr != end_ptr; ++ptr)
    {
    const uword index = (*ptr);
    
    const std::complex<T>& X = A_mem[index];
    const std::complex<T>& Y = B_values[col_offset];
    
    const T a = X.real();
    const T b = X.imag();
    
    const T c = Y.real();
    const T d = Y.imag();
    
    acc_real += (a*c) - (b*d);
    acc_imag += (a*d) + (b*c);
    
    ++col_offset;
    }
  
  return std::complex<T>(acc_real, acc_imag);
  }



template<typename T1, typename T2>
inline
void
glue_times_dense_sparse::apply(Mat<typename T1::elem_type>& out, const SpToDGlue<T1,T2,glue_times_dense_sparse>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_op_diagmat<T1>::value)  { out = SpMat<eT>(expr.A) * expr.B; return; }  // SpMat has specialised handling for op_diagmat
  
  const quasi_unwrap<T1> UA(expr.A);
  
  if(UA.is_alias(out))
    {
    Mat<eT> tmp;
    
    glue_times_dense_sparse::apply_noalias(tmp, UA.M, expr.B);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_times_dense_sparse::apply_noalias(out, UA.M, expr.B);
    }
  }



template<typename T1, typename T2>
inline
void
glue_times_dense_sparse::apply_noalias(Mat<typename T1::elem_type>& out, const T1& x, const T2& y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UA(x);
  const Mat<eT>&     A = UA.M;
  
  const unwrap_spmat<T2> UB(y);
  const SpMat<eT>&   B = UB.M;
  
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  out.set_size(A.n_rows, B.n_cols);
  
  if((A.n_elem == 0) || (B.n_nonzero == 0))  { out.zeros(); return; }
  
  if((resolves_to_rowvector<T1>::value) || (A.n_rows == 1))
    {
    arma_extra_debug_print("using row vector specialisation");
    
    if( (arma_config::openmp) && (mp_thread_limit::in_parallel() == false) && (B.n_cols >= 2) && mp_gate<eT>::eval(B.n_nonzero) )
      {
      #if defined(ARMA_USE_OPENMP)
        {
        arma_extra_debug_print("openmp implementation");
        
              eT* out_mem = out.memptr();
        const eT*   A_mem =   A.memptr();
        
        const uword B_n_cols  = B.n_cols;
        const int   n_threads = mp_thread_limit::get();
        
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for(uword col=0; col < B_n_cols; ++col)
          {
          out_mem[col] = dense_sparse_helper::dot(A_mem, B, col);
          }
        }
      #endif
      }
    else
      {
      arma_extra_debug_print("serial implementation");
      
            eT* out_mem = out.memptr();
      const eT*   A_mem =   A.memptr();
      
      const uword B_n_cols = B.n_cols;
      
      for(uword col=0; col < B_n_cols; ++col)
        {
        out_mem[col] = dense_sparse_helper::dot(A_mem, B, col);
        }
      }
    }
  else
  if( (arma_config::openmp) && (mp_thread_limit::in_parallel() == false) && (A.n_rows <= (A.n_cols / uword(100))) )
    {
    #if defined(ARMA_USE_OPENMP)
      {
      arma_extra_debug_print("using parallelised multiplication");
      
      const uword B_n_cols  = B.n_cols;
      const int   n_threads = mp_thread_limit::get();
      
      #pragma omp parallel for schedule(static) num_threads(n_threads)
      for(uword i=0; i < B_n_cols; ++i)
        {
        const uword col_offset_1 = B.col_ptrs[i  ];
        const uword col_offset_2 = B.col_ptrs[i+1];
        
        const uword col_offset_delta = col_offset_2 - col_offset_1;
        
        const uvec    indices(const_cast<uword*>(&(B.row_indices[col_offset_1])), col_offset_delta, false, false);
        const Col<eT>   B_col(const_cast<   eT*>(&(     B.values[col_offset_1])), col_offset_delta, false, false);
        
        out.col(i) = A.cols(indices) * B_col;
        }
      }
    #endif
    }
  else
    {
    arma_extra_debug_print("using standard multiplication");
    
    out.zeros();
    
    typename SpMat<eT>::const_iterator B_it = B.begin();
    
    const uword nnz        = B.n_nonzero;
    const uword out_n_rows = out.n_rows;
    
    for(uword count = 0; count < nnz; ++count, ++B_it)
      {
      const eT    B_it_val = (*B_it);
      const uword B_it_col = B_it.col();
      const uword B_it_row = B_it.row();
      
      const eT*   A_col =   A.colptr(B_it_row);
            eT* out_col = out.colptr(B_it_col);
      
      for(uword row = 0; row < out_n_rows; ++row)
        {
        out_col[row] += A_col[row] * B_it_val;
        }
      }
    }
  }



template<typename T1, typename T2>
inline
void
glue_times_dense_sparse::apply_mixed(Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >& out, const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  if( (is_same_type<eT1,out_eT>::no) && (is_same_type<eT2,out_eT>::yes) )
    {
    // upgrade T1
    
    const quasi_unwrap<T1> UA(X);
    const unwrap_spmat<T2> UB(Y);
    
    const   Mat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    const Mat<out_eT> AA = conv_to< Mat<out_eT> >::from(A);
    
    const SpMat<out_eT>& BB = reinterpret_cast< const SpMat<out_eT>& >(B);
    
    glue_times_dense_sparse::apply_noalias(out, AA, BB);
    }
  else
  if( (is_same_type<eT1,out_eT>::yes) && (is_same_type<eT2,out_eT>::no) )
    {
    // upgrade T2 
    
    const quasi_unwrap<T1> UA(X);
    const unwrap_spmat<T2> UB(Y);
    
    const   Mat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    const Mat<out_eT>& AA = reinterpret_cast< const Mat<out_eT>& >(A);
    
    SpMat<out_eT> BB(arma_layout_indicator(), B);
    
    for(uword i=0; i < B.n_nonzero; ++i)  { access::rw(BB.values[i]) = out_eT(B.values[i]); }
    
    glue_times_dense_sparse::apply_noalias(out, AA, BB);
    }
  else
    {
    // upgrade T1 and T2
    
    const quasi_unwrap<T1> UA(X);
    const unwrap_spmat<T2> UB(Y);
    
    const   Mat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    const Mat<out_eT> AA = conv_to< Mat<out_eT> >::from(A);
    
    SpMat<out_eT> BB(arma_layout_indicator(), B);
    
    for(uword i=0; i < B.n_nonzero; ++i)  { access::rw(BB.values[i]) = out_eT(B.values[i]); }
    
    glue_times_dense_sparse::apply_noalias(out, AA, BB);
    }
  }



//



template<typename T1, typename T2>
inline
void
glue_times_sparse_dense::apply(Mat<typename T1::elem_type>& out, const SpToDGlue<T1,T2,glue_times_sparse_dense>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_op_diagmat<T2>::value)  { out = expr.A * SpMat<eT>(expr.B); return; }  // SpMat has specialised handling for op_diagmat
  
  const quasi_unwrap<T2> UB(expr.B);
  
  if((sp_strip_trans<T1>::do_htrans && is_cx<eT>::no) || (sp_strip_trans<T1>::do_strans))
    {
    arma_extra_debug_print("detected non-conjugate transpose of A");
    
    const sp_strip_trans<T1> x_strip(expr.A);
    
    if(UB.is_alias(out))
      {
      Mat<eT> tmp;
      
      glue_times_sparse_dense::apply_noalias_trans(tmp, x_strip.M, UB.M);
      
      out.steal_mem(tmp);
      }
    else
      {
      glue_times_sparse_dense::apply_noalias_trans(out, x_strip.M, UB.M);
      }
    }
  else
    {
    if(UB.is_alias(out))
      {
      Mat<eT> tmp;
      
      glue_times_sparse_dense::apply_noalias(tmp, expr.A, UB.M);
      
      out.steal_mem(tmp);
      }
    else
      {
      glue_times_sparse_dense::apply_noalias(out, expr.A, UB.M);
      }
    }
  }



template<typename T1, typename T2>
inline
void
glue_times_sparse_dense::apply_noalias(Mat<typename T1::elem_type>& out, const T1& x, const T2& y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(x);
  const SpMat<eT>&   A = UA.M;
  
  const quasi_unwrap<T2> UB(y);
  const Mat<eT>&     B = UB.M;
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  arma_debug_assert_mul_size(A_n_rows, A_n_cols, B_n_rows, B_n_cols, "matrix multiplication");
  
  if((resolves_to_colvector<T2>::value) || (B_n_cols == 1))
    {
    arma_extra_debug_print("using column vector specialisation");
    
    out.zeros(A_n_rows, 1);
    
          eT* out_mem = out.memptr();
    const eT*   B_mem =   B.memptr();
    
    typename SpMat<eT>::const_iterator A_it = A.begin();
    
    const uword nnz = A.n_nonzero;
    
    for(uword count = 0; count < nnz; ++count, ++A_it)
      {
      const eT    A_it_val = (*A_it);
      const uword A_it_row = A_it.row();
      const uword A_it_col = A_it.col();
      
      out_mem[A_it_row] += A_it_val * B_mem[A_it_col];
      }
    }
  else
  if(B_n_cols >= (B_n_rows / uword(100)))
    {
    arma_extra_debug_print("using transpose-based multiplication");
    
    const SpMat<eT> At = A.st();
    const   Mat<eT> Bt = B.st();
    
    if(A_n_rows == B_n_cols)
      {
      glue_times_dense_sparse::apply_noalias(out, Bt, At);
      
      op_strans::apply_mat(out, out);  // since 'out' is square-sized, this will do an inplace transpose
      }
    else
      {
      Mat<eT> tmp;
      
      glue_times_dense_sparse::apply_noalias(tmp, Bt, At);
      
      op_strans::apply_mat(out, tmp);
      }
    }
  else
    {
    arma_extra_debug_print("using standard multiplication");
    
    out.zeros(A_n_rows, B_n_cols);
    
    typename SpMat<eT>::const_iterator A_it = A.begin();
    
    const uword nnz = A.n_nonzero;
    
    for(uword count = 0; count < nnz; ++count, ++A_it)
      {
      const eT    A_it_val = (*A_it);
      const uword A_it_row = A_it.row();
      const uword A_it_col = A_it.col();
      
      for(uword col = 0; col < B_n_cols; ++col)
        {
        out.at(A_it_row, col) += A_it_val * B.at(A_it_col, col);
        }
      }
    }
  }



template<typename T1, typename T2>
inline
void
glue_times_sparse_dense::apply_noalias_trans(Mat<typename T1::elem_type>& out, const T1& x, const T2& y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(x);
  const SpMat<eT>&   A = UA.M;  // NOTE: this is the given matrix without the transpose operation applied
  
  const quasi_unwrap<T2> UB(y);
  const Mat<eT>&     B = UB.M;
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  arma_debug_assert_mul_size(A_n_cols, A_n_rows, B_n_rows, B_n_cols, "matrix multiplication");
  
  if((resolves_to_colvector<T2>::value) || (B_n_cols == 1))
    {
    arma_extra_debug_print("using column vector specialisation (avoiding transpose of A)");
    
    if( (arma_config::openmp) && (mp_thread_limit::in_parallel() == false) && (A_n_cols >= 2) && mp_gate<eT>::eval(A.n_nonzero) )
      {
      arma_extra_debug_print("opemp implementation");
      
      #if defined(ARMA_USE_OPENMP)
        {
        out.zeros(A_n_cols, 1);
        
              eT* out_mem = out.memptr();
        const eT*   B_mem =   B.memptr();
        
        const int n_threads = mp_thread_limit::get();
        
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for(uword col=0; col < A_n_cols; ++col)
          {
          out_mem[col] = dense_sparse_helper::dot(B_mem, A, col);
          }
        }
      #endif
      }
    else
      {
      arma_extra_debug_print("serial implementation");
      
      out.zeros(A_n_cols, 1);
      
            eT* out_mem = out.memptr();
      const eT*   B_mem =   B.memptr();
      
      for(uword col=0; col < A_n_cols; ++col)
        {
        out_mem[col] = dense_sparse_helper::dot(B_mem, A, col);
        }
      }
    }
  else
  if(B_n_cols >= (B_n_rows / uword(100)))
    {
    arma_extra_debug_print("using transpose-based multiplication (avoiding transpose of A)");
    
    const Mat<eT> Bt = B.st();
    
    if(A_n_cols == B_n_cols)
      {
      glue_times_dense_sparse::apply_noalias(out, Bt, A);
      
      op_strans::apply_mat(out, out);  // since 'out' is square-sized, this will do an inplace transpose
      }
    else
      {
      Mat<eT> tmp;
      
      glue_times_dense_sparse::apply_noalias(tmp, Bt, A);
      
      op_strans::apply_mat(out, tmp);
      }
    }
  else
    {
    arma_extra_debug_print("using standard multiplication (avoiding transpose of A)");
    
    out.zeros(A_n_cols, B_n_cols);
    
    typename SpMat<eT>::const_iterator A_it = A.begin();
    
    const uword nnz = A.n_nonzero;
    
    for(uword count = 0; count < nnz; ++count, ++A_it)
      {
      const eT    A_it_val = (*A_it);
      const uword A_it_row = A_it.row();
      const uword A_it_col = A_it.col();
      
      for(uword col = 0; col < B_n_cols; ++col)
        {
        out.at(A_it_col, col) += A_it_val * B.at(A_it_row, col);
        }
      }
    }
  }



template<typename T1, typename T2>
inline
void
glue_times_sparse_dense::apply_mixed(Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >& out, const T1& X, const T2& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  if( (is_same_type<eT1,out_eT>::no) && (is_same_type<eT2,out_eT>::yes) )
    {
    // upgrade T1
    
    const unwrap_spmat<T1> UA(X);
    const quasi_unwrap<T2> UB(Y);
    
    const SpMat<eT1>& A = UA.M;
    const   Mat<eT2>& B = UB.M;
    
    SpMat<out_eT> AA(arma_layout_indicator(), A);
    
    for(uword i=0; i < A.n_nonzero; ++i)  { access::rw(AA.values[i]) = out_eT(A.values[i]); }
    
    const Mat<out_eT>& BB = reinterpret_cast< const Mat<out_eT>& >(B);
    
    glue_times_sparse_dense::apply_noalias(out, AA, BB);
    }
  else
  if( (is_same_type<eT1,out_eT>::yes) && (is_same_type<eT2,out_eT>::no) )
    {
    // upgrade T2 
    
    const unwrap_spmat<T1> UA(X);
    const quasi_unwrap<T2> UB(Y);
    
    const SpMat<eT1>& A = UA.M;
    const   Mat<eT2>& B = UB.M;
    
    const SpMat<out_eT>& AA = reinterpret_cast< const SpMat<out_eT>& >(A);
    
    const Mat<out_eT> BB = conv_to< Mat<out_eT> >::from(B);
    
    glue_times_sparse_dense::apply_noalias(out, AA, BB);
    }
  else
    {
    // upgrade T1 and T2
    
    const unwrap_spmat<T1> UA(X);
    const quasi_unwrap<T2> UB(Y);
    
    const SpMat<eT1>& A = UA.M;
    const   Mat<eT2>& B = UB.M;
    
    SpMat<out_eT> AA(arma_layout_indicator(), A);
    
    for(uword i=0; i < A.n_nonzero; ++i)  { access::rw(AA.values[i]) = out_eT(A.values[i]); }
    
    const Mat<out_eT> BB = conv_to< Mat<out_eT> >::from(B);
    
    glue_times_sparse_dense::apply_noalias(out, AA, BB);
    }
  }



//! @}
