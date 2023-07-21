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


//! \addtogroup fn_diags_spdiags
//! @{



template<typename T1, typename T2>
inline
Mat<typename T1::elem_type>
diags(const Base<typename T1::elem_type, T1>& V_expr, const Base<sword,T2>& D_expr, const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UV(V_expr.get_ref());
  const Mat<eT>& V     = UV.M;
  
  const quasi_unwrap<T2> UD(D_expr.get_ref());
  const Mat<sword>& D  = UD.M;
  
  arma_debug_check( ((D.is_vec() == false) && (D.is_empty() == false)), "D must be a vector" );
  
  arma_debug_check( (V.n_cols != D.n_elem), "number of colums in matrix V must match the length of vector D" );
  
  Mat<eT> out(n_rows, n_cols, fill::zeros);
  
  for(uword i=0; i < D.n_elem; ++i)
    {
    const sword diag_id = D[i];
    
    const uword row_offset = (diag_id < 0) ? uword(-diag_id) : 0;
    const uword col_offset = (diag_id > 0) ? uword( diag_id) : 0;
    
    arma_debug_check_bounds
      (
      ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
      "diags(): requested diagonal out of bounds"
      );
    
    const uword diag_len = (std::min)(n_rows - row_offset, n_cols - col_offset);
    
    const uword V_start = (diag_id < 0) ? uword(0) : uword(diag_id);
    
    const eT* V_colmem = V.colptr(i);
    
    for(uword j=0; j < diag_len; ++j)
      {
      const uword V_index = V_start + j;
      
      if(V_index >= V.n_rows)  { break; }
      
      out.at(j + row_offset, j + col_offset) = V_colmem[V_index];
      }
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
SpMat<typename T1::elem_type>
spdiags(const Base<typename T1::elem_type, T1>& V_expr, const Base<sword,T2>& D_expr, const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UV(V_expr.get_ref());
  const Mat<eT>& V     = UV.M;
  
  const quasi_unwrap<T2> UD(D_expr.get_ref());
  const Mat<sword>& D  = UD.M;
  
  arma_debug_check( ((D.is_vec() == false) && (D.is_empty() == false)), "D must be a vector" );
  
  arma_debug_check( (V.n_cols != D.n_elem), "number of colums in matrix V must match the length of vector D" );
  
  MapMat<eT> tmp(n_rows, n_cols);
  
  for(uword i=0; i < D.n_elem; ++i)
    {
    const sword diag_id = D[i];
    
    const uword row_offset = (diag_id < 0) ? uword(-diag_id) : 0;
    const uword col_offset = (diag_id > 0) ? uword( diag_id) : 0;
    
    arma_debug_check_bounds
      (
      ((row_offset > 0) && (row_offset >= n_rows)) || ((col_offset > 0) && (col_offset >= n_cols)),
      "diags(): requested diagonal out of bounds"
      );
    
    const uword diag_len = (std::min)(n_rows - row_offset, n_cols - col_offset);
    
    const uword V_start = (diag_id < 0) ? uword(0) : uword(diag_id);
    
    const eT* V_colmem = V.colptr(i);
    
    for(uword j=0; j < diag_len; ++j)
      {
      const uword V_index = V_start + j;
      
      if(V_index >= V.n_rows)  { break; }
      
      tmp.at(j + row_offset, j + col_offset) = V_colmem[V_index];
      }
    }
  
  return SpMat<eT>(tmp);
  }



//! @}
