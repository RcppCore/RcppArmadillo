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



//! \addtogroup spop_shift
//! @{



template<typename eT>
inline
void
spop_shift::apply_noalias(SpMat<eT>& out, const SpMat<eT>& X, const uword len, const uword neg, const uword dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check_bounds( ((dim == 0) && (len >= X.n_rows)), "shift(): shift amount out of bounds" );
  arma_debug_check_bounds( ((dim == 1) && (len >= X.n_cols)), "shift(): shift amount out of bounds" );
  
  if(X.n_nonzero == 0)  { out.zeros(X.n_rows, X.n_cols); return; }
  
  if(len == 0)  { out = X; return; }
  
  typename SpMat<eT>::const_iterator it     = X.begin();
  typename SpMat<eT>::const_iterator it_end = X.end();
  
  Mat<uword> locs(2, X.n_nonzero, arma_nozeros_indicator());
  
  uword* locs_mem = locs.memptr();
  
  if(dim == 0)
    {
    const uword X_row_threshold = X.n_rows - len;
    
    for(; it != it_end; ++it)
      {
      const uword X_row = it.row();
            uword Y_row = 0;
      
           if(neg == 0)  { Y_row = (X_row <  X_row_threshold) ? (X_row + len) : (X_row - X_row_threshold); }
      else if(neg == 1)  { Y_row = (X_row >= len            ) ? (X_row - len) : (X_row + X_row_threshold); }
      
      locs_mem[0] = Y_row;
      locs_mem[1] = it.col();
      
      locs_mem += 2;
      }
    }
  else
  if(dim == 1)
    {
    const uword X_col_threshold = X.n_cols - len;
    
    for(; it != it_end; ++it)
      {
      const uword X_col = it.col();
            uword Y_col = 0;
      
           if(neg == 0)  { Y_col = (X_col <  X_col_threshold) ? (X_col + len) : (X_col - X_col_threshold); }
      else if(neg == 1)  { Y_col = (X_col >= len            ) ? (X_col - len) : (X_col + X_col_threshold); }
      
      locs_mem[0] = it.row();
      locs_mem[1] = Y_col;
      
      locs_mem += 2;
      }
    }
  
  const Col<eT> vals(const_cast<eT*>(X.values), X.n_nonzero, false);
  
  SpMat<eT> Y(locs, vals, X.n_rows, X.n_cols, true, false);
  
  out.steal_mem(Y);
  }



//! @}
