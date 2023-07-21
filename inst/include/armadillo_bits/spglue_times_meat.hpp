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


//! \addtogroup spglue_times
//! @{



template<typename T1, typename T2>
inline
void
spglue_times::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(X.A);
  const unwrap_spmat<T2> UB(X.B);
  
  const bool is_alias = (UA.is_alias(out) || UB.is_alias(out));
  
  if(is_alias == false)
    {
    spglue_times::apply_noalias(out, UA.M, UB.M);
    }
  else
    {
    SpMat<eT> tmp;
    
    spglue_times::apply_noalias(tmp, UA.M, UB.M);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
spglue_times::apply(SpMat<typename T1::elem_type>& out, const SpGlue<SpOp<T1,spop_scalar_times>,T2,spglue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_spmat<T1> UA(X.A.m);
  const unwrap_spmat<T2> UB(X.B);
  
  const bool is_alias = (UA.is_alias(out) || UB.is_alias(out));
  
  if(is_alias == false)
    {
    spglue_times::apply_noalias(out, UA.M, UB.M);
    }
  else
    {
    SpMat<eT> tmp;
    
    spglue_times::apply_noalias(tmp, UA.M, UB.M);
    
    out.steal_mem(tmp);
    }
  
  out *= X.A.aux;
  }



template<typename eT>
inline
void
spglue_times::apply_noalias(SpMat<eT>& c, const SpMat<eT>& x, const SpMat<eT>& y)
  {
  arma_extra_debug_sigprint();
  
  const uword x_n_rows = x.n_rows;
  const uword x_n_cols = x.n_cols;
  const uword y_n_rows = y.n_rows;
  const uword y_n_cols = y.n_cols;
  
  arma_debug_assert_mul_size(x_n_rows, x_n_cols, y_n_rows, y_n_cols, "matrix multiplication");
  
  // First we must determine the structure of the new matrix (column pointers).
  // This follows the algorithm described in 'Sparse Matrix Multiplication
  // Package (SMMP)' (R.E. Bank and C.C. Douglas, 2001).  Their description of
  // "SYMBMM" does not include anything about memory allocation.  In addition it
  // does not consider that there may be elements which space may be allocated
  // for but which evaluate to zero anyway.  So we have to modify the algorithm
  // to work that way.  For the "SYMBMM" implementation we will not determine
  // the row indices but instead just the column pointers.
  
  //SpMat<typename T1::elem_type> c(x_n_rows, y_n_cols); // Initializes col_ptrs to 0.
  c.zeros(x_n_rows, y_n_cols);
  
  //if( (x.n_elem == 0) || (y.n_elem == 0) )  { return; }
  if( (x.n_nonzero == 0) || (y.n_nonzero == 0) )  { return; }
  
  // Auxiliary storage which denotes when items have been found.
  podarray<uword> index(x_n_rows);
  index.fill(x_n_rows); // Fill with invalid links.
  
  typename SpMat<eT>::const_iterator y_it  = y.begin();
  typename SpMat<eT>::const_iterator y_end = y.end();
  
  // SYMBMM: calculate column pointers for resultant matrix to obtain a good
  // upper bound on the number of nonzero elements.
  uword cur_col_length = 0;
  uword last_ind = x_n_rows + 1;
  do
    {
    const uword y_it_row = y_it.row();
    
    // Look through the column that this point (*y_it) could affect.
    typename SpMat<eT>::const_iterator x_it = x.begin_col_no_sync(y_it_row);
    
    while(x_it.col() == y_it_row)
      {
      const uword x_it_row = x_it.row();
      
      // A point at x(i, j) and y(j, k) implies a point at c(i, k).
      if(index[x_it_row] == x_n_rows)
        {
        index[x_it_row] = last_ind;
        last_ind = x_it_row;
        ++cur_col_length;
        }
      
      ++x_it;
      }
    
    const uword old_col = y_it.col();
    ++y_it;
    
    // See if column incremented.
    if(old_col != y_it.col())
      {
      // Set column pointer (this is not a cumulative count; that is done later).
      access::rw(c.col_ptrs[old_col + 1]) = cur_col_length;
      cur_col_length = 0;
      
      // Return index markers to zero.  Use last_ind for traversal.
      while(last_ind != x_n_rows + 1)
        {
        const uword tmp = index[last_ind];
        index[last_ind] = x_n_rows;
        last_ind = tmp;
        }
      }
    }
  while(y_it != y_end);
  
  // Accumulate column pointers.
  for(uword i = 0; i < c.n_cols; ++i)
    {
    access::rw(c.col_ptrs[i + 1]) += c.col_ptrs[i];
    }
  
  // Now that we know a decent bound on the number of nonzero elements,
  // allocate the memory and fill it.
  
  const uword max_n_nonzero = c.col_ptrs[c.n_cols];
  
  c.mem_resize(max_n_nonzero);
  
  // Now the implementation of the NUMBMM algorithm.
  uword cur_pos = 0; // Current position in c matrix.
  podarray<eT> sums(x_n_rows); // Partial sums.
  sums.zeros();
  
  podarray<uword> sorted_indices(x_n_rows);  // upper bound
  
  // last_ind is already set to x_n_rows, and cur_col_length is already set to 0.
  // We will loop through all columns as necessary.
  uword cur_col = 0;
  while(cur_col < c.n_cols)
    {
    // Skip to next column with elements in it.
    while((cur_col < c.n_cols) && (c.col_ptrs[cur_col] == c.col_ptrs[cur_col + 1]))
      {
      // Update current column pointer to actual number of nonzero elements up
      // to this point.
      access::rw(c.col_ptrs[cur_col]) = cur_pos;
      ++cur_col;
      }
    
    if(cur_col == c.n_cols)  { break; }
    
    // Update current column pointer.
    access::rw(c.col_ptrs[cur_col]) = cur_pos;
    
    // Check all elements in this column.
    typename SpMat<eT>::const_iterator y_col_it = y.begin_col_no_sync(cur_col);
    
    while(y_col_it.col() == cur_col)
      {
      const uword y_col_it_row = y_col_it.row();
      
      // Check all elements in the column of the other matrix corresponding to
      // the row of this column.
      typename SpMat<eT>::const_iterator x_col_it = x.begin_col_no_sync(y_col_it_row);
      
      const eT y_value = (*y_col_it);
      
      while(x_col_it.col() == y_col_it_row)
        {
        const uword x_col_it_row = x_col_it.row();
        
        // A point at x(i, j) and y(j, k) implies a point at c(i, k).
        // Add to partial sum.
        const eT x_value = (*x_col_it);
        sums[x_col_it_row] += (x_value * y_value);
        
        // Add point if it hasn't already been marked.
        if(index[x_col_it_row] == x_n_rows)
          {
          index[x_col_it_row] = last_ind;
          last_ind = x_col_it_row;
          }
        
        ++x_col_it;
        }
      
      ++y_col_it;
      }
    
    // Now sort the indices that were used in this column.
    uword cur_index = 0;
    while(last_ind != x_n_rows + 1)
      {
      const uword tmp = last_ind;
      
      // Check that it wasn't a "fake" nonzero element.
      if(sums[tmp] != eT(0))
        {
        // Assign to next open position.
        sorted_indices[cur_index] = tmp;
        ++cur_index;
        }

      last_ind = index[tmp];
      index[tmp] = x_n_rows;
      }
    
    // Now sort the indices.
    if(cur_index != 0)
      {
      op_sort::direct_sort_ascending(sorted_indices.memptr(), cur_index);

      for(uword k = 0; k < cur_index; ++k)
        {
        const uword row = sorted_indices[k];
        access::rw(c.row_indices[cur_pos]) = row;
        access::rw(c.values[cur_pos]) = sums[row];
        sums[row] = eT(0);
        ++cur_pos;
        }
      }

    // Move to next column.
    ++cur_col;
    }
  
  // Update last column pointer and resize to actual memory size.
  
  // access::rw(c.col_ptrs[c.n_cols]) = cur_pos;
  // c.mem_resize(cur_pos);
  
  access::rw(c.col_ptrs[c.n_cols]) = cur_pos;
  
  if(cur_pos < max_n_nonzero)  { c.mem_resize(cur_pos); }
  }



//
//
//



template<typename T1, typename T2>
inline
void
spglue_times_mixed::apply(SpMat<typename eT_promoter<T1,T2>::eT>& out, const mtSpGlue<typename eT_promoter<T1,T2>::eT, T1, T2, spglue_times_mixed>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename eT_promoter<T1,T2>::eT out_eT;
  
  if( (is_same_type<eT1,out_eT>::no) && (is_same_type<eT2,out_eT>::yes) )
    {
    // upgrade T1
    
    const unwrap_spmat<T1> UA(expr.A);
    const unwrap_spmat<T2> UB(expr.B);
    
    const SpMat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    SpMat<out_eT> AA(arma_layout_indicator(), A);
    
    for(uword i=0; i < A.n_nonzero; ++i)  { access::rw(AA.values[i]) = out_eT(A.values[i]); }
    
    const SpMat<out_eT>& BB = reinterpret_cast< const SpMat<out_eT>& >(B);
    
    out = AA * BB;
    }
  else
  if( (is_same_type<eT1,out_eT>::yes) && (is_same_type<eT2,out_eT>::no) )
    {
    // upgrade T2 
    
    const unwrap_spmat<T1> UA(expr.A);
    const unwrap_spmat<T2> UB(expr.B);
    
    const SpMat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    const SpMat<out_eT>& AA = reinterpret_cast< const SpMat<out_eT>& >(A);
    
    SpMat<out_eT> BB(arma_layout_indicator(), B);
    
    for(uword i=0; i < B.n_nonzero; ++i)  { access::rw(BB.values[i]) = out_eT(B.values[i]); }
    
    out = AA * BB;
    }
  else
    {
    // upgrade T1 and T2
    
    const unwrap_spmat<T1> UA(expr.A);
    const unwrap_spmat<T2> UB(expr.B);
    
    const SpMat<eT1>& A = UA.M;
    const SpMat<eT2>& B = UB.M;
    
    SpMat<out_eT> AA(arma_layout_indicator(), A);
    SpMat<out_eT> BB(arma_layout_indicator(), B);
    
    for(uword i=0; i < A.n_nonzero; ++i)  { access::rw(AA.values[i]) = out_eT(A.values[i]); }
    for(uword i=0; i < B.n_nonzero; ++i)  { access::rw(BB.values[i]) = out_eT(B.values[i]); }
    
    out = AA * BB;
    }
  }



//! @}
