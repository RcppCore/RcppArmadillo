// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_reshape
//! @{



template<typename T1>
inline
void
op_reshape::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_reshape>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   A_tmp(in.m);
  const Mat<eT>& A = A_tmp.M;
  
  const bool is_alias = (&out == &A);
  
  const uword in_n_rows = in.aux_uword_a;
  const uword in_n_cols = in.aux_uword_b;
  const uword in_dim    = in.aux_uword_c;
  
  const uword in_n_elem = in_n_rows * in_n_cols;
  
  if(A.n_elem == in_n_elem)
    {
    if(in_dim == 0)
      {
      if(is_alias == false)
        {
        out.set_size(in_n_rows, in_n_cols);
        arrayops::copy( out.memptr(), A.memptr(), out.n_elem );
        }
      else  // &out == &A, i.e. inplace resize
        {
        const bool same_size = ( (out.n_rows == in_n_rows) && (out.n_cols == in_n_cols) );
        
        if(same_size == false)
          {
          arma_debug_check
            (
            (out.mem_state == 3),
            "reshape(): size can't be changed as template based size specification is in use"
            );
          
          access::rw(out.n_rows) = in_n_rows;
          access::rw(out.n_cols) = in_n_cols;
          }
        }
      }
    else
      {
      unwrap_check< Mat<eT> > B_tmp(A, is_alias);
      const Mat<eT>& B      = B_tmp.M;
      
      out.set_size(in_n_rows, in_n_cols);
      
      eT* out_mem = out.memptr();
      uword i = 0;
      
      const uword B_n_rows = B.n_rows;
      const uword B_n_cols = B.n_cols;
      
      for(uword row=0; row<B_n_rows; ++row)
      for(uword col=0; col<B_n_cols; ++col)
        {
        out_mem[i] = B.at(row,col);
        ++i;
        }
      }
    }
  else
    {
    const unwrap_check< Mat<eT> > B_tmp(A, is_alias);
    const Mat<eT>& B            = B_tmp.M;
    
    const uword n_elem_to_copy = (std::min)(B.n_elem, in_n_elem);
    
    out.set_size(in_n_rows, in_n_cols);
    
    eT* out_mem = out.memptr();
    
    if(in_dim == 0)
      {
      arrayops::copy( out_mem, B.memptr(), n_elem_to_copy );
      }
    else
      {
      uword row = 0;
      uword col = 0;
      
      const uword B_n_cols = B.n_cols;
      
      for(uword i=0; i<n_elem_to_copy; ++i)
        {
        out_mem[i] = B.at(row,col);
        
        ++col;
        
        if(col >= B_n_cols)
          {
          col = 0;
          ++row;
          }
        }
      }
    
    for(uword i=n_elem_to_copy; i<in_n_elem; ++i)
      {
      out_mem[i] = eT(0);
      }
    
    }
  }



template<typename T1>
inline
void
op_reshape::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_reshape>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> A_tmp(in.m);
  const Cube<eT>& A   = A_tmp.M;
  
  const uword in_n_rows   = in.aux_uword_a;
  const uword in_n_cols   = in.aux_uword_b;
  const uword in_n_slices = in.aux_uword_c;
  const uword in_dim      = in.aux_uword_d;
  
  const uword in_n_elem = in_n_rows * in_n_cols * in_n_slices;
  
  if(A.n_elem == in_n_elem)
    {
    if(in_dim == 0)
      {
      if(&out != &A)
        {
        out.set_size(in_n_rows, in_n_cols, in_n_slices);
        arrayops::copy( out.memptr(), A.memptr(), out.n_elem );
        }
      else  // &out == &A, i.e. inplace resize
        {
        const bool same_size = ( (out.n_rows == in_n_rows) && (out.n_cols == in_n_cols) && (out.n_slices == in_n_slices) );
        
        if(same_size == false)
          {
          arma_debug_check
            (
            (out.mem_state == 3),
            "reshape(): size can't be changed as template based size specification is in use"
            );
          
          out.delete_mat();
          
          access::rw(out.n_rows)       = in_n_rows;
          access::rw(out.n_cols)       = in_n_cols;
          access::rw(out.n_elem_slice) = in_n_rows * in_n_cols;
          access::rw(out.n_slices)     = in_n_slices;
          
          out.create_mat();
          }
        }
      }
    else
      {
      unwrap_cube_check< Cube<eT> > B_tmp(A, out);
      const Cube<eT>& B           = B_tmp.M;
      
      out.set_size(in_n_rows, in_n_cols, in_n_slices);
      
      eT* out_mem = out.memptr();
      uword i = 0;
      
      const uword B_n_rows   = B.n_rows;
      const uword B_n_cols   = B.n_cols;
      const uword B_n_slices = B.n_slices;
      
      for(uword slice=0; slice<B_n_slices; ++slice)
      for(uword row=0; row<B_n_rows; ++row)
      for(uword col=0; col<B_n_cols; ++col)
        {
        out_mem[i] = B.at(row,col,slice);
        ++i;
        }
      }
    }
  else
    {
    const unwrap_cube_check< Cube<eT> > B_tmp(A, out);
    const Cube<eT>& B                 = B_tmp.M;
    
    const uword n_elem_to_copy = (std::min)(B.n_elem, in_n_elem);
    
    out.set_size(in_n_rows, in_n_cols, in_n_slices);
    
    eT* out_mem = out.memptr();
    
    if(in_dim == 0)
      {
      arrayops::copy( out_mem, B.memptr(), n_elem_to_copy );
      }
    else
      {
      uword row   = 0;
      uword col   = 0;
      uword slice = 0;
      
      const uword B_n_rows = B.n_rows;
      const uword B_n_cols = B.n_cols;
      
      for(uword i=0; i<n_elem_to_copy; ++i)
        {
        out_mem[i] = B.at(row,col,slice);
        
        ++col;
        
        if(col >= B_n_cols)
          {
          col = 0;
          ++row;
          
          if(row >= B_n_rows)
            {
            row = 0;
            ++slice;
            }
          }
        }
      }
    
    for(uword i=n_elem_to_copy; i<in_n_elem; ++i)
      {
      out_mem[i] = eT(0);
      }
    
    }
  }



//! @}
