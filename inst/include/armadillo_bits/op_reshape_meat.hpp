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



//! \addtogroup op_reshape
//! @{



template<typename T1>
inline
void
op_reshape::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_reshape>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword new_n_rows = in.aux_uword_a;
  const uword new_n_cols = in.aux_uword_b;
  
  if(is_Mat<T1>::value)
    {
    const unwrap<T1>   U(in.m);
    const Mat<eT>& A = U.M;
    
    if(&out == &A)
      {
      op_reshape::apply_mat_inplace(out, new_n_rows, new_n_cols);
      }
    else
      {
      op_reshape::apply_mat_noalias(out, A, new_n_rows, new_n_cols);
      }
    }
  else
  if((is_Mat<typename Proxy<T1>::stored_type>::value) || (arma_config::openmp && Proxy<T1>::use_mp))
    {
    const quasi_unwrap<T1> U(in.m);
    
    if(U.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_reshape::apply_mat_noalias(tmp, U.M, new_n_rows, new_n_cols);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_reshape::apply_mat_noalias(out, U.M, new_n_rows, new_n_cols);
      }
    }
  else
    {
    const Proxy<T1> P(in.m);
    
    if(P.is_alias(out))
      {
      Mat<eT> tmp;
      
      op_reshape::apply_proxy_noalias(tmp, P, new_n_rows, new_n_cols);
      
      out.steal_mem(tmp);
      }
    else
      {
      op_reshape::apply_proxy_noalias(out, P, new_n_rows, new_n_cols);
      }
    }
  }



template<typename eT>
inline
void
op_reshape::apply_mat_inplace(Mat<eT>& A, const uword new_n_rows, const uword new_n_cols)
  {
  arma_debug_sigprint();
  
  if( (A.n_rows == new_n_rows) && (A.n_cols == new_n_cols) )  { return; }
  
  arma_conform_check( (A.vec_state == 1) && (new_n_cols != 1), "reshape(): requested size is not compatible with column vector layout" );
  arma_conform_check( (A.vec_state == 2) && (new_n_rows != 1), "reshape(): requested size is not compatible with row vector layout"    );
  
  if(A.is_empty())  { A.zeros(new_n_rows, new_n_cols); return; }
  
  const bool is_into_empty  = ( (new_n_cols == uword(0)) || (new_n_rows == uword(0)) );
  const bool is_into_colvec = ( (new_n_cols == uword(1)) && (new_n_rows == A.n_elem) );
  const bool is_into_rowvec = ( (new_n_rows == uword(1)) && (new_n_cols == A.n_elem) );
  const bool is_rowcol_swap = ( (new_n_cols == A.n_rows) && (new_n_rows == A.n_cols) );
  
  if(is_into_empty || is_into_colvec || is_into_rowvec || is_rowcol_swap)  { A.set_size(new_n_rows, new_n_cols); return; }
  
  Mat<eT> B(new_n_rows, new_n_cols, arma_nozeros_indicator());
  
  op_reshape::apply_mat_noalias(B, A, new_n_rows, new_n_cols);
  
  A.steal_mem(B);
  }



template<typename eT>
inline
void
op_reshape::apply_mat_noalias(Mat<eT>& out, const Mat<eT>& A, const uword new_n_rows, const uword new_n_cols)
  {
  arma_debug_sigprint();
  
  out.set_size(new_n_rows, new_n_cols);
  
  const uword n_elem_to_copy = (std::min)(A.n_elem, out.n_elem);
  
  eT* out_mem = out.memptr();
  
  arrayops::copy( out_mem, A.memptr(), n_elem_to_copy );
  
  if(n_elem_to_copy < out.n_elem)
    {
    const uword n_elem_leftover = out.n_elem - n_elem_to_copy;
    
    arrayops::fill_zeros(&(out_mem[n_elem_to_copy]), n_elem_leftover);
    }
  }



template<typename T1>
inline
void
op_reshape::apply_proxy_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword new_n_rows, const uword new_n_cols)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out.set_size(new_n_rows, new_n_cols);
  
  const uword n_elem_to_copy = (std::min)(P.get_n_elem(), out.n_elem);
  
  eT* out_mem = out.memptr();
  
  if(Proxy<T1>::use_at == false)
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
    
    for(uword i=0; i < n_elem_to_copy; ++i)  { out_mem[i] = Pea[i]; }
    }
  else
    {
    uword i = 0;
    
    const uword P_n_rows = P.get_n_rows();
    const uword P_n_cols = P.get_n_cols();
    
    for(uword col=0; col < P_n_cols; ++col)
    for(uword row=0; row < P_n_rows; ++row)
      {
      if(i >= n_elem_to_copy)  { goto nested_loop_end; }
      
      out_mem[i] = P.at(row,col);
      
      ++i;
      }
    
    nested_loop_end: ;
    }
  
  if(n_elem_to_copy < out.n_elem)
    {
    const uword n_elem_leftover = out.n_elem - n_elem_to_copy;
    
    arrayops::fill_zeros(&(out_mem[n_elem_to_copy]), n_elem_leftover);
    }
  }



template<typename T1>
inline
void
op_reshape::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_reshape>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> U(in.m);
  const Cube<eT>& A   = U.M;
  
  const uword new_n_rows   = in.aux_uword_a;
  const uword new_n_cols   = in.aux_uword_b;
  const uword new_n_slices = in.aux_uword_c;
  
  if(&out == &A)
    {
    op_reshape::apply_cube_inplace(out, new_n_rows, new_n_cols, new_n_slices);
    }
  else
    {
    op_reshape::apply_cube_noalias(out, A, new_n_rows, new_n_cols, new_n_slices);
    }
  }
  


template<typename eT>
inline
void
op_reshape::apply_cube_inplace(Cube<eT>& A, const uword new_n_rows, const uword new_n_cols, const uword new_n_slices)
  {
  arma_debug_sigprint();
  
  if( (A.n_rows == new_n_rows) && (A.n_cols == new_n_cols) && (A.n_slices == new_n_slices) )  { return; }
  
  if(A.is_empty())  { A.zeros(new_n_rows, new_n_cols, new_n_slices); return; }
  
  const bool is_into_empty  = ( (new_n_cols == uword(0)) || (new_n_rows == uword(0)) || (new_n_slices == uword(0)  ) );
  const bool is_into_colvec = ( (new_n_cols == uword(1)) && (new_n_rows == A.n_elem) && (new_n_slices == uword(1)  ) );
  const bool is_into_rowvec = ( (new_n_rows == uword(1)) && (new_n_cols == A.n_elem) && (new_n_slices == uword(1)  ) );
  const bool is_rowcol_swap = ( (new_n_cols == A.n_rows) && (new_n_rows == A.n_cols) && (new_n_slices == A.n_slices) );
  
  if(is_into_empty || is_into_colvec || is_into_rowvec || is_rowcol_swap)  { A.set_size(new_n_rows, new_n_cols, new_n_slices); return; }
  
  Cube<eT> B(new_n_rows, new_n_cols, new_n_slices, arma_nozeros_indicator());
  
  op_reshape::apply_cube_noalias(B, A, new_n_rows, new_n_cols, new_n_slices);
  
  A.steal_mem(B);
  }



template<typename eT>
inline
void
op_reshape::apply_cube_noalias(Cube<eT>& out, const Cube<eT>& A, const uword new_n_rows, const uword new_n_cols, const uword new_n_slices)
  {
  arma_debug_sigprint();
  
  out.set_size(new_n_rows, new_n_cols, new_n_slices);
  
  const uword n_elem_to_copy = (std::min)(A.n_elem, out.n_elem);
  
  eT* out_mem = out.memptr();
  
  arrayops::copy( out_mem, A.memptr(), n_elem_to_copy );
  
  if(n_elem_to_copy < out.n_elem)
    {
    const uword n_elem_leftover = out.n_elem - n_elem_to_copy;
    
    arrayops::fill_zeros(&(out_mem[n_elem_to_copy]), n_elem_leftover);
    }
  }



//! @}
