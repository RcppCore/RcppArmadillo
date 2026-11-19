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


//! \addtogroup fn_solve_lse
//! @{



template<typename T1, typename T2, typename T3, typename T4>
inline
bool
solve_lse_worker
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A_expr,
  const Base<typename T1::elem_type,T2>& B_expr,
  const Base<typename T1::elem_type,T3>& C_expr,
  const Base<typename T1::elem_type,T4>& D_expr
  )
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // A, B, C, D are overwritten by auxlib::solve_lse()
  
  Mat<eT> A = A_expr.get_ref();
  Mat<eT> B = B_expr.get_ref();
  Mat<eT> C = C_expr.get_ref();
  Mat<eT> D = D_expr.get_ref();
  
  arma_conform_check( (B.n_cols != uword(1)), "solve_lse(): B must be a column vector" );
  arma_conform_check( (D.n_cols != uword(1)), "solve_lse(): D must be a column vector" );
  
  arma_conform_check( (A.n_rows != B.n_rows), "solve_lse(): number of rows in A and B must be the same" );
  arma_conform_check( (C.n_rows != D.n_rows), "solve_lse(): number of rows in C and D must be the same" );
  
  arma_conform_check( (A.n_cols != C.n_cols), "solve_lse(): number of columns in A and C must be the same" );
  
  arma_conform_check( (C.n_rows > C.n_cols), "solve_lse(): C.n_rows must be <= C.n_cols" );
  
  arma_conform_check( (A.n_cols > (A.n_rows + C.n_cols)), "solve_lse(): A.n_cols must be <= (A.n_rows + C.n_cols)" );
  
  out.set_size(A.n_cols, 1);
  
  if(out.is_empty() || A.is_empty() || B.is_empty() || C.is_empty() || D.is_empty())  { out.zeros(); return true; }
    
  const bool status = auxlib::solve_lse(out, A, C, B, D);  // B and C are deliberately swapped
  
  if(status == false)  { out.soft_reset(); }
  
  return status;
  }



template<typename T1, typename T2, typename T3, typename T4>
inline
typename enable_if2< is_blas_type<typename T1::elem_type>::value, bool >::result
solve_lse
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A_expr,
  const Base<typename T1::elem_type,T2>& B_expr,
  const Base<typename T1::elem_type,T3>& C_expr,
  const Base<typename T1::elem_type,T4>& D_expr
  )
  {
  arma_debug_sigprint();
  
  const bool status = solve_lse_worker(out, A_expr, B_expr, C_expr, D_expr);
  
  if(status == false)  { arma_warn(3, "solve_lse(): solution not found"); }
  
  return status;
  }



template<typename T1, typename T2, typename T3, typename T4>
arma_warn_unused
inline
Col<typename T1::elem_type>
solve_lse
  (
  const Base<typename T1::elem_type,T1>& A_expr,
  const Base<typename T1::elem_type,T2>& B_expr,
  const Base<typename T1::elem_type,T3>& C_expr,
  const Base<typename T1::elem_type,T4>& D_expr
  )
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Col<eT> out;
  
  const bool status = solve_lse_worker(out, A_expr, B_expr, C_expr, D_expr);
  
  if(status == false)  { arma_stop_runtime_error("solve_lse(): solution not found"); }
  
  return out;
  }



//! @}
