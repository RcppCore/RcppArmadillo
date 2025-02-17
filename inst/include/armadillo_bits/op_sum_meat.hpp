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


//! \addtogroup op_sum
//! @{



template<typename T1>
inline
void
op_sum::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sum>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  arma_conform_check( (dim > 1), "sum(): parameter 'dim' must be 0 or 1" );
  
  if((is_Mat<T1>::value) || (is_Mat<typename Proxy<T1>::stored_type>::value) || (arma_config::openmp && Proxy<T1>::use_mp))
    {
    const quasi_unwrap<T1> U(in.m);
    
    if(U.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_sum::apply_mat_noalias(tmp, U.M, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_mat_noalias(out, U.M, dim);
      }
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_sum::apply_proxy_noalias(tmp, P, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_proxy_noalias(out, P, dim);
      }
    }
  }



template<typename T1>
inline
void
op_sum::apply(Mat<typename T1::elem_type>& out, const Op< eOp<T1,eop_square>, op_sum >& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  typedef eOp<T1,eop_square> inner_expr_type;
  
  typedef typename inner_expr_type::proxy_type::stored_type inner_expr_P_stored_type;
  
  const uword dim = in.aux_uword_a;
  
  arma_conform_check( (dim > 1), "sum(): parameter 'dim' must be 0 or 1" );
  
  if(is_Mat<inner_expr_P_stored_type>::value)
    {
    const quasi_unwrap<inner_expr_P_stored_type> U(in.m.P.Q);
    
    if(U.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_sum::apply_mat_square_noalias(tmp, U.M, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_mat_square_noalias(out, U.M, dim);
      }
    }
  else
  if(arma_config::openmp && Proxy<inner_expr_type>::use_mp)
    {
    const quasi_unwrap<inner_expr_type> U(in.m);  // force evaluation of compound inner expression
    
    if(U.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_sum::apply_mat_noalias(tmp, U.M, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_mat_noalias(out, U.M, dim);
      }
    }
  else
    {
    const Proxy<inner_expr_type> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_sum::apply_proxy_noalias(tmp, P, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_proxy_noalias(out, P, dim);
      }
    }
  }



template<typename T1>
inline
void
op_sum::apply(Mat<typename T1::elem_type>& out, const Op< eOp<T1,eop_pow>, op_sum >& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(in.m.aux == eT(2))
    {
    typedef Op< eOp<T1,eop_square>, op_sum > modified_whole_expr_type;
    
    op_sum::apply(out, reinterpret_cast<const modified_whole_expr_type& >(in) );
    
    return;
    }
  
  if((in.m.aux == eT(0.5)) && is_non_integral<eT>::value)
    {
    typedef Op< eOp<T1,eop_sqrt>, op_sum > modified_whole_expr_type;
    
    op_sum::apply(out, reinterpret_cast<const modified_whole_expr_type& >(in) );
    
    return;
    }
  
  typedef eOp<T1,eop_pow> inner_expr_type;
  
  typedef typename inner_expr_type::proxy_type::stored_type inner_expr_P_stored_type;
  
  const uword dim = in.aux_uword_a;
  
  arma_conform_check( (dim > 1), "sum(): parameter 'dim' must be 0 or 1" );
  
  if( (is_Mat<inner_expr_P_stored_type>::value) || (arma_config::openmp && Proxy<inner_expr_type>::use_mp) )
    {
    const quasi_unwrap<inner_expr_type> U(in.m);  // force evaluation of eop_pow
    
    if(U.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_sum::apply_mat_noalias(tmp, U.M, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_mat_noalias(out, U.M, dim);
      }
    }
  else
    {
    const Proxy<inner_expr_type> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_sum::apply_proxy_noalias(tmp, P, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_proxy_noalias(out, P, dim);
      }
    }
  }



template<typename eT>
inline
void
op_sum::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim)
  {
  arma_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword out_n_rows = (dim == 0) ? uword(1) : X_n_rows;
  const uword out_n_cols = (dim == 0) ? X_n_cols : uword(1);
  
  out.set_size(out_n_rows, out_n_cols);
  
  if(X.n_elem == 0)  { out.zeros(); return; }
  
  const eT* X_colptr =   X.memptr();
        eT* out_mem  = out.memptr();
  
  if(dim == 0)
    {
    for(uword col=0; col < X_n_cols; ++col)
      {
      out_mem[col] = arrayops::accumulate( X_colptr, X_n_rows );
      
      X_colptr += X_n_rows;
      }
    }
  else
    {
    arrayops::copy(out_mem, X_colptr, X_n_rows);
    
    X_colptr += X_n_rows;
    
    for(uword col=1; col < X_n_cols; ++col)
      {
      arrayops::inplace_plus( out_mem, X_colptr, X_n_rows );
      
      X_colptr += X_n_rows;
      }
    }
  }



template<typename eT>
inline
void
op_sum::apply_mat_square_noalias(Mat<eT>& out, const Mat<eT>& X, const uword dim)
  {
  arma_debug_sigprint();
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword out_n_rows = (dim == 0) ? uword(1) : X_n_rows;
  const uword out_n_cols = (dim == 0) ? X_n_cols : uword(1);
  
  out.set_size(out_n_rows, out_n_cols);
  
  if(X.n_elem == 0)  { out.zeros(); return; }
  
  const eT* X_colptr =   X.memptr();
        eT* out_mem  = out.memptr();
  
  if(dim == 0)
    {
    for(uword col=0; col < X_n_cols; ++col)
      {
      out_mem[col] = op_dot::direct_dot(X_n_rows, X_colptr, X_colptr);
      
      X_colptr += X_n_rows;
      }
    }
  else
    {
    for(uword row=0; row < X_n_rows; ++row)  { const eT tmp = X_colptr[row]; out_mem[row] = tmp*tmp; }
    
    X_colptr += X_n_rows;
    
    for(uword col=1; col < X_n_cols; ++col)
      {
      for(uword row=0; row < X_n_rows; ++row)  { const eT tmp = X_colptr[row]; out_mem[row] += tmp*tmp; }
      
      X_colptr += X_n_rows;
      }
    }
  }



template<typename T1>
inline
void
op_sum::apply_proxy_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword dim)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword P_n_rows = P.get_n_rows();
  const uword P_n_cols = P.get_n_cols();
  
  const uword out_n_rows = (dim == 0) ? uword(1) : P_n_rows;
  const uword out_n_cols = (dim == 0) ? P_n_cols : uword(1);
  
  out.set_size(out_n_rows, out_n_cols);
  
  if(P.get_n_elem() == 0)  { out.zeros(); return; }
  
  eT* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    if(dim == 0)
      {
      uword count = 0;
      
      for(uword col=0; col < P_n_cols; ++col)
        {
        eT val1 = eT(0);
        eT val2 = eT(0);
        
        uword j;
        for(j=1; j < P_n_rows; j+=2)
          {
          val1 += P[count]; ++count;
          val2 += P[count]; ++count;
          }
        
        if((j-1) < P_n_rows)
          {
          val1 += P[count]; ++count;
          }
        
        out_mem[col] = (val1 + val2);
        }
      }
    else
      {
      uword count = 0;
      
      for(uword row=0; row < P_n_rows; ++row)
        {
        out_mem[row] = P[count]; ++count;
        }
      
      for(uword col=1; col < P_n_cols; ++col)
      for(uword row=0; row < P_n_rows; ++row)
        {
        out_mem[row] += P[count]; ++count;
        }
      }
    }
  else
    {
    if(dim == 0)
      {
      for(uword col=0; col < P_n_cols; ++col)
        {
        eT val1 = eT(0);
        eT val2 = eT(0);
        
        uword i,j;
        for(i=0, j=1; j < P_n_rows; i+=2, j+=2)
          {
          val1 += P.at(i,col);
          val2 += P.at(j,col);
          }
        
        if(i < P_n_rows)
          {
          val1 += P.at(i,col);
          }
        
        out_mem[col] = (val1 + val2);
        }
      }
    else
      {
      for(uword row=0; row < P_n_rows; ++row)
        {
        out_mem[row] = P.at(row,0);
        }
      
      for(uword col=1; col < P_n_cols; ++col)
      for(uword row=0; row < P_n_rows; ++row)
        {
        out_mem[row] += P.at(row,col);
        }
      }
    }
  }



//
// cubes



template<typename T1>
inline
void
op_sum::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_sum>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  arma_conform_check( (dim > 2), "sum(): parameter 'dim' must be 0 or 1 or 2" );
  
  if((is_Cube<T1>::value) || (is_Cube<typename ProxyCube<T1>::stored_type>::value) || (arma_config::openmp && ProxyCube<T1>::use_mp))
    {
    const unwrap_cube<T1> U(in.m);
    
    if(U.is_alias(out))
      {
      Cube<eT> tmp;
      
      op_sum::apply_cube_noalias(tmp, U.M, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_cube_noalias(out, U.M, dim);
      }
    }
  else
    {
    const ProxyCube<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Cube<eT> tmp;
      
      op_sum::apply_proxy_noalias(tmp, P, dim);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_sum::apply_proxy_noalias(out, P, dim);
      }
    }
  }



template<typename eT>
inline
void
op_sum::apply_cube_noalias(Cube<eT>& out, const Cube<eT>& X, const uword dim)
  {
  arma_debug_sigprint();
  
  const uword X_n_rows   = X.n_rows;
  const uword X_n_cols   = X.n_cols;
  const uword X_n_slices = X.n_slices;
  
  if(dim == 0)
    {
    out.set_size(1, X_n_cols, X_n_slices);
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      eT* out_mem = out.slice_memptr(slice);
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        out_mem[col] = arrayops::accumulate( X.slice_colptr(slice,col), X_n_rows );
        }
      }
    }
  else
  if(dim == 1)
    {
    out.zeros(X_n_rows, 1, X_n_slices);
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      eT* out_mem = out.slice_memptr(slice);
      
      for(uword col=0; col < X_n_cols; ++col)
        {
        arrayops::inplace_plus( out_mem, X.slice_colptr(slice,col), X_n_rows );
        }
      }
    }
  else
  if(dim == 2)
    {
    out.zeros(X_n_rows, X_n_cols, 1);
    
    eT* out_mem = out.memptr();
    
    for(uword slice=0; slice < X_n_slices; ++slice)
      {
      arrayops::inplace_plus(out_mem, X.slice_memptr(slice), X.n_elem_slice );
      }
    }
  }



template<typename T1>
inline
void
op_sum::apply_proxy_noalias(Cube<typename T1::elem_type>& out, const ProxyCube<T1>& P, const uword dim)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword P_n_rows   = P.get_n_rows();
  const uword P_n_cols   = P.get_n_cols();
  const uword P_n_slices = P.get_n_slices();
  
  if(dim == 0)
    {
    out.set_size(1, P_n_cols, P_n_slices);
    
    for(uword slice=0; slice < P_n_slices; ++slice)
      {
      eT* out_mem = out.slice_memptr(slice);
      
      for(uword col=0; col < P_n_cols; ++col)
        {
        eT val1 = eT(0);
        eT val2 = eT(0);
        
        uword i,j;
        for(i=0, j=1; j < P_n_rows; i+=2, j+=2)
          {
          val1 += P.at(i,col,slice);
          val2 += P.at(j,col,slice);
          }
        
        if(i < P_n_rows)
          {
          val1 += P.at(i,col,slice);
          }
        
        out_mem[col] = (val1 + val2);
        }
      }
    }
  else
  if(dim == 1)
    {
    out.zeros(P_n_rows, 1, P_n_slices);
    
    for(uword slice=0; slice < P_n_slices; ++slice)
      {
      eT* out_mem = out.slice_memptr(slice);
      
      for(uword col=0; col < P_n_cols; ++col)
      for(uword row=0; row < P_n_rows; ++row)
        {
        out_mem[row] += P.at(row,col,slice);
        }
      }
    }
  else
  if(dim == 2)
    {
    out.zeros(P_n_rows, P_n_cols, 1);
    
    for(uword slice=0; slice < P_n_slices; ++slice)
      {
      for(uword col=0; col < P_n_cols; ++col)
        {
        eT* out_mem = out.slice_colptr(0,col);
        
        for(uword row=0; row < P_n_rows; ++row)
          {
          out_mem[row] += P.at(row,col,slice);
          }
        }
      }
    }
  }



//! @}
