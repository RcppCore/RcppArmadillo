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


//! \addtogroup glue_polyval
//! @{



template<typename eT>
inline
void
glue_polyval::apply_noalias(Mat<eT>& out, const Mat<eT>& P, const Mat<eT>& X)
  {
  arma_debug_sigprint();
  
  out.set_size(X.n_rows, X.n_cols);
  
  const eT*   P_mem    = P.memptr();
  const uword P_n_elem = P.n_elem;
  
  out.fill(P_mem[0]);

  for(uword i=1; i < P_n_elem; ++i)
    {
    out = out % X + P_mem[i];
    }
  }



template<typename T1, typename T2>
inline
void
glue_polyval::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_polyval>& expr)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UP(expr.A);
  const quasi_unwrap<T2> UX(expr.B);
  
  arma_conform_check( ((UP.M.is_vec() == false) && (UP.M.is_empty() == false)), "polyval(): argument P must be a vector" );
  
  if(UP.M.is_empty() || UX.M.is_empty())  { out.zeros(UX.M.n_rows, UX.M.n_cols); return; }
  
  if(UP.is_alias(out) || UX.is_alias(out))
    {
    Mat<eT> tmp;
    
    glue_polyval::apply_noalias(tmp, UP.M, UX.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_polyval::apply_noalias(out, UP.M, UX.M);
    }
  }



template<typename T1, typename T2>
inline
void
glue_polyval::apply(Mat_noalias<typename T1::elem_type>& out, const Glue<T1,T2,glue_polyval>& expr)
  {
  arma_debug_sigprint();
  
  const quasi_unwrap<T1> UP(expr.A);
  const quasi_unwrap<T2> UX(expr.B);
  
  arma_conform_check( ((UP.M.is_vec() == false) && (UP.M.is_empty() == false)), "polyval(): argument P must be a vector" );
  
  if(UP.M.is_empty() || UX.M.is_empty())  { out.zeros(UX.M.n_rows, UX.M.n_cols); return; }
  
  glue_polyval::apply_noalias(out, UP.M, UX.M);
  }



//! @}
