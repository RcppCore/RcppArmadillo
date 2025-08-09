// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------



//! \addtogroup op_orth_null
//! @{



template<typename T1>
inline
void
op_orth::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_orth>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const T tol = access::tmp_real(expr.aux);
  
  const bool status = op_orth::apply_direct(out, expr.m, tol);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("orth(): svd failed");
    }
  }



template<typename T1>
inline
bool
op_orth::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, typename T1::pod_type tol)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  arma_conform_check((tol < T(0)), "orth(): tolerance must be >= 0");
  
  Mat<eT> U;
  Col< T> s;
  Mat<eT> V;
  
  Mat<eT> A(expr.get_ref());
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  const uword N_limit = (is_cx<eT>::yes) ? uword(20000) : uword(23000);
  
  const bool allow_dc = (sizeof(blas_int) >= std::size_t(8)) ? true : (N <= N_limit);
  
  if(allow_dc == false)
    {
    arma_warn(3, "orth(): matrix size too large for divide-and-conquer algorithm; using standard algorithm instead");
    }
  
  const bool status = (allow_dc) ? auxlib::svd_dc(U, s, V, A) : auxlib::svd(U, s, V, A);
  
  V.reset();
  
  if(status == false)  { return false; }
  
  if(s.is_empty())  { out.reset(); return true; }
  
  const uword s_n_elem = s.n_elem;
  const T*    s_mem    = s.memptr();
  
  // set tolerance to default if it hasn't been specified
  if(tol == T(0))  { tol = (std::max)(A.n_rows, A.n_cols) * s_mem[0] * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i < s_n_elem; ++i)  { count += (s_mem[i] > tol) ? uword(1) : uword(0); }
  
  if(count > 0)
    {
    out = U.head_cols(count);  // out *= eT(-1);
    }
  else
    {
    out.set_size(A.n_rows, 0);
    }
  
  return true;
  }



//



template<typename T1>
inline
void
op_null::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_null>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const T tol = access::tmp_real(expr.aux);
  
  const bool status = op_null::apply_direct(out, expr.m, tol);
  
  if(status == false)
    {
    out.soft_reset();
    arma_stop_runtime_error("null(): svd failed");
    }
  }



template<typename T1>
inline
bool
op_null::apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& expr, typename T1::pod_type tol)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  arma_conform_check((tol < T(0)), "null(): tolerance must be >= 0");
  
  Mat<eT> U;
  Col< T> s;
  Mat<eT> V;
  
  Mat<eT> A(expr.get_ref());
  
  const uword N = (std::min)(A.n_rows, A.n_cols);
  
  const uword N_limit = (is_cx<eT>::yes) ? uword(20000) : uword(23000);
  
  const bool allow_dc = (sizeof(blas_int) >= std::size_t(8)) ? true : (N <= N_limit);
  
  if(allow_dc == false)
    {
    arma_warn(3, "null(): matrix size too large for divide-and-conquer algorithm; using standard algorithm instead");
    }
  
  const bool status = (allow_dc) ? auxlib::svd_dc(U, s, V, A) : auxlib::svd(U, s, V, A);
  
  U.reset();
  
  if(status == false)  { return false; }
  
  if(s.is_empty())  { out.reset(); return true; }
  
  const uword s_n_elem = s.n_elem;
  const T*    s_mem    = s.memptr();
  
  // set tolerance to default if it hasn't been specified
  if(tol == T(0))  { tol = (std::max)(A.n_rows, A.n_cols) * s_mem[0] * std::numeric_limits<T>::epsilon(); }
  
  uword count = 0;
  
  for(uword i=0; i < s_n_elem; ++i)  { count += (s_mem[i] > tol) ? uword(1) : uword(0); }
  
  if(count < A.n_cols)
    {
    out = V.tail_cols(A.n_cols - count);
    
    const uword out_n_elem = out.n_elem;
          eT*   out_mem    = out.memptr();
    
    for(uword i=0; i<out_n_elem; ++i)
      {
      if(std::abs(out_mem[i]) < std::numeric_limits<T>::epsilon())  { out_mem[i] = eT(0); }
      }
    }
  else
    {
    out.set_size(A.n_cols, 0);
    }
  
  return true;
  }



//! @}
