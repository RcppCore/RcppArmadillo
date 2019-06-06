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


//! \addtogroup spglue_merge
//! @{



template<typename eT>
arma_hot
inline
void
spglue_merge::apply(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  const bool is_alias = A.is_alias(out) || B.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_merge::apply_noalias(out, A, B);
    }
  else
    {
    SpMat<eT> tmp;
    
    spglue_merge::apply_noalias(tmp, A, B);
    
    out.steal_mem(tmp);
    }
  }



//! merge under the assumption that there is no overlap of non-zero elements between A and B
//! this is a simplified form of spglue_plus::apply_noalias()
template<typename eT>
arma_hot
inline
void
spglue_merge::apply_noalias(SpMat<eT>& out, const SpMat<eT>& A, const SpMat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "merge");
  
  const uword merge_n_nonzero = A.n_nonzero + B.n_nonzero;
  
  if(merge_n_nonzero == 0)  { out.zeros(A.n_rows, A.n_cols); return; }
  
  if(A.n_nonzero == 0)  { out = B; return; }
  if(B.n_nonzero == 0)  { out = A; return; }
  
  out.reserve(A.n_rows, A.n_cols, merge_n_nonzero);
  
  typename SpMat<eT>::const_iterator x_it  = A.begin();
  typename SpMat<eT>::const_iterator x_end = A.end();
  
  typename SpMat<eT>::const_iterator y_it  = B.begin();
  typename SpMat<eT>::const_iterator y_end = B.end();
  
  uword count = 0;
  
  while( (x_it != x_end) || (y_it != y_end) )
    {
    eT out_val;
    
    const uword x_it_col = x_it.col();
    const uword x_it_row = x_it.row();
    
    const uword y_it_col = y_it.col();
    const uword y_it_row = y_it.row();
    
    bool use_y_loc = false;
    
    if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
      {
      out_val = (*x_it);
      
      ++x_it;
      }
    else
      {
      out_val = (*y_it);
      
      ++y_it;
      
      use_y_loc = true;
      }
    
    access::rw(out.values[count]) = out_val;
    
    const uword out_row = (use_y_loc == false) ? x_it_row : y_it_row;
    const uword out_col = (use_y_loc == false) ? x_it_col : y_it_col;
    
    access::rw(out.row_indices[count]) = out_row;
    access::rw(out.col_ptrs[out_col + 1])++;
    ++count;
    }
  
  arma_check( (count != merge_n_nonzero), "spglue_merge::apply_noalias(): internal error: count != merge_n_nonzero" );
  
  const uword out_n_cols = out.n_cols;
  
  uword* col_ptrs = access::rwp(out.col_ptrs);
  
  for(uword c = 1; c <= out_n_cols; ++c)
    {
    col_ptrs[c] += col_ptrs[c - 1];
    }
  }



//! @}
