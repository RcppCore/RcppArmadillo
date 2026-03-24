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


//! \addtogroup op_flip
//! @{



template<typename T1>
inline
void
op_flipud::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_flipud>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Mat<T1>::value)
    {
    const unwrap<T1> U(in.m);
    
    if(&out == &(U.M))  { op_flipud::apply_mat_inplace(out); return; }
    
    // fallthrough if operation is not inplace
    }
  
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_flipud::apply_mat_noalias(tmp, U.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_flipud::apply_mat_noalias(out, U.M);
    }
  }



template<typename T1>
inline
void
op_flipud::apply(Mat_noalias<typename T1::elem_type>& out, const Op<T1,op_flipud>& in)
  {
  arma_debug_sigprint();
  
  const quasi_unwrap<T1> U(in.m);
  
  op_flipud::apply_mat_noalias(out, U.M);
  }



template<typename eT>
inline
void
op_flipud::apply_mat_inplace(Mat<eT>& X)
  {
  arma_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword X_n_rows_m1 = X_n_rows - 1;
  
  const uword N = X_n_rows / 2;
  
  if(X_n_cols == 1)
    {
    eT* X_mem = X.memptr();
    
    for(uword row=0; row < N; ++row)
      {
      std::swap(X_mem[X_n_rows_m1 - row], X_mem[row]);
      }
    }
  else
    {
    for(uword col=0; col < X_n_cols; ++col)
      {
      eT* X_colmem = X.colptr(col);
      
      for(uword row=0; row < N; ++row)
        {
        std::swap(X_colmem[X_n_rows_m1 - row], X_colmem[row]);
        }
      }
    }
  }



template<typename eT>
inline
void
op_flipud::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword X_n_rows_m1 = X_n_rows - 1;
  
  out.set_size(X_n_rows, X_n_cols);
  
  if(X_n_cols == 1)
    {
    const eT*   X_mem =   X.memptr();
          eT* out_mem = out.memptr();
    
    for(uword row=0; row < X_n_rows; ++row)
      {
      out_mem[X_n_rows_m1 - row] = X_mem[row];
      }
    }
  else
    {
    for(uword col=0; col < X_n_cols; ++col)
      {
      const eT*   X_colmem =   X.colptr(col);
            eT* out_colmem = out.colptr(col);
      
      for(uword row=0; row < X_n_rows; ++row)
        {
        out_colmem[X_n_rows_m1 - row] = X_colmem[row];
        }
      }
    }
  }



//



template<typename T1>
inline
void
op_fliplr::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_fliplr>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(is_Mat<T1>::value)
    {
    const unwrap<T1> U(in.m);
    
    if(&out == &(U.M))  { op_fliplr::apply_mat_inplace(out); return; }
    
    // fallthrough if operation is not inplace
    }
    
  const quasi_unwrap<T1> U(in.m);
  
  if(U.is_alias(out))
    {
    Mat<eT> tmp;
    
    op_fliplr::apply_mat_noalias(tmp, U.M);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_fliplr::apply_mat_noalias(out, U.M);
    }
  }



template<typename T1>
inline
void
op_fliplr::apply(Mat_noalias<typename T1::elem_type>& out, const Op<T1,op_fliplr>& in)
  {
  arma_debug_sigprint();
  
  const quasi_unwrap<T1> U(in.m);
  
  op_fliplr::apply_mat_noalias(out, U.M);
  }



template<typename eT>
inline
void
op_fliplr::apply_mat_inplace(Mat<eT>& X)
  {
  arma_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword X_n_cols_m1 = X_n_cols - 1;
  
  const uword N = X_n_cols / 2;
  
  if(X_n_rows == 1)
    {
    eT* X_mem = X.memptr();
    
    for(uword col=0; col < N; ++col)
      {
      std::swap(X_mem[X_n_cols_m1 - col], X_mem[col]);
      }
    }
  else
    {
    for(uword col=0; col < N; ++col)
      {
      X.swap_cols(X_n_cols_m1 - col, col);
      }
    }
  }



template<typename eT>
inline
void
op_fliplr::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword X_n_cols_m1 = X_n_cols - 1;
  
  out.set_size(X_n_rows, X_n_cols);
  
  if(X_n_rows == 1)
    {
    const eT*   X_mem =   X.memptr();
          eT* out_mem = out.memptr();
    
    for(uword col=0; col < X_n_cols; ++col)
      {
      out_mem[X_n_cols_m1 - col] = X_mem[col];
      }
    }
  else
    {
    for(uword col=0; col < X_n_cols; ++col)
      {
      out.col(X_n_cols_m1 - col) = X.col(col);
      }
    }
  }



//! @}
