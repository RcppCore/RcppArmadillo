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


//! \addtogroup op_sort
//! @{



template<typename eT>
inline 
void
op_sort::direct_sort(eT* X, const uword n_elem, const uword sort_mode)
  {
  arma_debug_sigprint();
  
  if(sort_mode == 0)
    {
    arma_lt_comparator<eT> comparator;
    
    std::sort(&X[0], &X[n_elem], comparator);
    }
  else
    {
    arma_gt_comparator<eT> comparator;
    
    std::sort(&X[0], &X[n_elem], comparator);
    }
  }



template<typename eT>
inline 
void
op_sort::direct_sort_ascending(eT* X, const uword n_elem)
  {
  arma_debug_sigprint();
  
  arma_lt_comparator<eT> comparator;
  
  std::sort(&X[0], &X[n_elem], comparator);
  }



template<typename eT>
inline
void
op_sort::apply_noalias(Mat<eT>& out, const Mat<eT>& X, const uword sort_mode, const uword dim)
  {
  arma_debug_sigprint();
  
  if(X.n_elem <= 1)  { out = X; return; }
  
  if(is_cx<eT>::yes)
    {
    arma_debug_print("op_sort::apply(): complex version");
    
    if(dim == 0)  // sort the contents of each column
      {
      arma_debug_print("op_sort::apply(): dim = 0");
      
      const uword n_rows = X.n_rows;
      const uword n_cols = X.n_cols;
      
      out.set_size(n_rows, n_cols);
      
      uvec indices;
      
      for(uword col=0; col < n_cols; ++col)
        {
        const Col<eT> X_col( const_cast<eT*>(X.colptr(col)), n_rows, false );
        
        op_sort_index::apply_noalias_mat(indices, X_col, sort_mode);
        
        const uword* indices_mem = indices.memptr();
        const    eT*   X_col_mem =   X_col.memptr();
                 eT* out_col_mem =  out.colptr(col);
        
        for(uword i=0; i < n_rows; ++i)  { out_col_mem[i] = X_col_mem[ indices_mem[i] ]; }
        }
      }
    else
    if(dim == 1)  // sort the contents of each row
      {
      arma_debug_print("op_sort::apply(): dim = 1");
      
      Mat<eT> Y;
      
      op_strans::apply_mat_noalias(Y, X);
      
      const uword n_rows = Y.n_rows;
      const uword n_cols = Y.n_cols;
      
      Mat<eT> tmp(n_rows, n_cols);
      
      uvec indices;
      
      for(uword col=0; col < n_cols; ++col)
        {
        const Col<eT> Y_col( const_cast<eT*>(Y.colptr(col)), n_rows, false );
        
        op_sort_index::apply_noalias_mat(indices, Y_col, sort_mode);
        
        const uword* indices_mem = indices.memptr();
        const    eT*   Y_col_mem =   Y_col.memptr();
                 eT* tmp_col_mem =  tmp.colptr(col);
        
        for(uword i=0; i < n_rows; ++i)  { tmp_col_mem[i] = Y_col_mem[ indices_mem[i] ]; }
        }
      
      op_strans::apply_mat_noalias(out, tmp);
      }
    }
  else
    {
    arma_debug_print("op_sort::apply(): plain version");
    
    if(dim == 0)  // sort the contents of each column
      {
      arma_debug_print("op_sort::apply(): dim = 0");
      
      out = X;
      
      const uword n_rows = out.n_rows;
      const uword n_cols = out.n_cols;
      
      for(uword col=0; col < n_cols; ++col)
        {
        op_sort::direct_sort( out.colptr(col), n_rows, sort_mode );
        }
      }
    else
    if(dim == 1)  // sort the contents of each row
      {
      if(X.n_rows == 1)  // a row vector
        {
        arma_debug_print("op_sort::apply(): dim = 1, vector specific");
        
        out = X;
        
        op_sort::direct_sort(out.memptr(), out.n_elem, sort_mode);
        }
      else  // not a row vector
        {
        arma_debug_print("op_sort::apply(): dim = 1, generic");
        
        Mat<eT> Y;
        
        op_strans::apply_mat_noalias(Y, X);
        
        const uword n_rows = Y.n_rows;
        const uword n_cols = Y.n_cols;
        
        for(uword col=0; col < n_cols; ++col)
          {
          op_sort::direct_sort( Y.colptr(col), n_rows, sort_mode );
          }
        
        op_strans::apply_mat_noalias(out, Y);
        }
      }
    }
  }



template<typename T1>
inline
void
op_sort::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sort>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> U(in.m);
  const Mat<eT>&     X = U.M;
  
  const uword sort_mode = in.aux_uword_a;
  const uword dim       = in.aux_uword_b;
  
  arma_conform_check( (sort_mode > 1),        "sort(): parameter 'sort_mode' must be 0 or 1" );
  arma_conform_check( (dim > 1),              "sort(): parameter 'dim' must be 0 or 1"       );
  arma_conform_check( (X.internal_has_nan()), "sort(): detected NaN"                         );
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_sort::apply_noalias(tmp, X, sort_mode, dim);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_sort::apply_noalias(out, X, sort_mode, dim);
    }
  }



template<typename T1>
inline
void
op_sort_vec::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sort_vec>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   U(in.m);  // not using quasi_unwrap, to ensure there is no aliasing with subviews
  const Mat<eT>& X = U.M;
  
  const uword sort_mode = in.aux_uword_a;
  
  arma_conform_check( (sort_mode > 1),        "sort(): parameter 'sort_mode' must be 0 or 1" );
  arma_conform_check( (X.internal_has_nan()), "sort(): detected NaN"                         );
  
  if(X.n_elem <= 1)  { out = X; return; }
  
  if(is_cx<eT>::yes)
    {
    uvec indices;
    
    op_sort_index::apply_noalias_mat(indices, X, sort_mode);
    
    const uword N = indices.n_elem;
    
    arma_check( (N != X.n_elem), "internal error: op_sort_vec::apply(): N != X.n_elem" );
    
    Mat<eT> tmp(X.n_rows, X.n_cols, arma_nozeros_indicator());  // in case there is aliasing
    
    const uword* indices_mem = indices.memptr();
    const    eT*       X_mem =       X.memptr();
             eT*     tmp_mem =     tmp.memptr();
    
    for(uword i=0; i < N; ++i)  { tmp_mem[i] = X_mem[ indices_mem[i] ]; }
    
    out.steal_mem(tmp);
    }
  else
    {
    out = X;  // not checking for aliasing, to allow inplace sorting of vectors
    
    eT* out_mem = out.memptr();
    
    eT* start_ptr =  out_mem;
    eT* endp1_ptr = &out_mem[out.n_elem];
    
    if(sort_mode == 0)
      {
      arma_lt_comparator<eT> comparator;
      
      std::sort(start_ptr, endp1_ptr, comparator);
      }
    else
      {
      arma_gt_comparator<eT> comparator;
      
      std::sort(start_ptr, endp1_ptr, comparator);
      }
    }
  }



//! @}
