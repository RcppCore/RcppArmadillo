// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_reshape
//! @{



template<typename T1>
inline
void
op_reshape::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_reshape>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  const u32 in_n_rows = in.aux_u32_a;
  const u32 in_n_cols = in.aux_u32_b;
  
  const u32 in_n_elem = in_n_rows * in_n_cols;
  
  if(A.n_elem == in_n_elem)
    {
    if(in.aux == eT(0))
      {
      if(&out != &A)
        {
        out.set_size(in_n_rows, in_n_cols);
        syslib::copy_elem( out.memptr(), A.memptr(), out.n_elem );
        }
      else
        {
        access::rw(out.n_rows) = in_n_rows;
        access::rw(out.n_cols) = in_n_cols;
        }
      }
    else
      {
      unwrap_check< Mat<eT> > tmp(A, out);
      const Mat<eT>& B      = tmp.M;

      out.set_size(in_n_rows, in_n_cols);
      
      eT* out_mem = out.memptr();
      u32 i = 0;
      
      for(u32 row=0; row<B.n_rows; ++row)
        {
        for(u32 col=0; col<B.n_cols; ++col)
          {
          out_mem[i] = B.at(row,col);
          ++i;
          }
        }
        
      }
    }
  else
    {
    const unwrap_check< Mat<eT> > tmp(A, out);
    const Mat<eT>& B            = tmp.M;
    
    const u32 n_elem_to_copy = (std::min)(B.n_elem, in_n_elem);
    
    out.set_size(in_n_rows, in_n_cols);
    
    eT* out_mem = out.memptr();
    
    if(in.aux == eT(0))
      {
      syslib::copy_elem( out_mem, B.memptr(), n_elem_to_copy );
      }
    else
      {
      u32 row = 0;
      u32 col = 0;
      
      for(u32 i=0; i<n_elem_to_copy; ++i)
        {
        out_mem[i] = B.at(row,col);
        
        ++col;
        
        if(col >= B.n_cols)
          {
          col = 0;
          ++row;
          }
        }
      }
    
    for(u32 i=n_elem_to_copy; i<in_n_elem; ++i)
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
  
  const unwrap_cube<T1> tmp(in.m);
  const Cube<eT>& A   = tmp.M;
  
  const u32 in_n_rows   = in.aux_u32_a;
  const u32 in_n_cols   = in.aux_u32_b;
  const u32 in_n_slices = in.aux_u32_c;
  
  const u32 in_n_elem = in_n_rows * in_n_cols * in_n_slices;
  
  if(A.n_elem == in_n_elem)
    {
    if(in.aux == eT(0))
      {
      if(&out != &A)
        {
        out.set_size(in_n_rows, in_n_cols, in_n_slices);
        syslib::copy_elem( out.memptr(), A.memptr(), out.n_elem );
        }
      else
        {
        access::rw(out.n_rows)   = in_n_rows;
        access::rw(out.n_cols)   = in_n_cols;
        access::rw(out.n_slices) = in_n_slices;
        }
      }
    else
      {
      unwrap_cube_check< Cube<eT> > tmp(A, out);
      const Cube<eT>& B           = tmp.M;
      
      out.set_size(in_n_rows, in_n_cols, in_n_slices);
      
      eT* out_mem = out.memptr();
      u32 i = 0;
      
      for(u32 slice=0; slice<B.n_slices; ++slice)
        {
        for(u32 row=0; row<B.n_rows; ++row)
          {
          for(u32 col=0; col<B.n_cols; ++col)
            {
            out_mem[i] = B.at(row,col,slice);
            ++i;
            }
          }
        }
        
      }
    }
  else
    {
    const unwrap_cube_check< Cube<eT> > tmp(A, out);
    const Cube<eT>& B                 = tmp.M;
    
    const u32 n_elem_to_copy = (std::min)(B.n_elem, in_n_elem);
    
    out.set_size(in_n_rows, in_n_cols, in_n_slices);
    
    eT* out_mem = out.memptr();
    
    if(in.aux == eT(0))
      {
      syslib::copy_elem( out_mem, B.memptr(), n_elem_to_copy );
      }
    else
      {
      u32 row   = 0;
      u32 col   = 0;
      u32 slice = 0;
      
      for(u32 i=0; i<n_elem_to_copy; ++i)
        {
        out_mem[i] = B.at(row,col,slice);
        
        ++col;
        
        if(col >= B.n_cols)
          {
          col = 0;
          ++row;
          
          if(row >= B.n_rows)
            {
            row = 0;
            ++slice;
            }
          }
        }
      }
    
    for(u32 i=n_elem_to_copy; i<in_n_elem; ++i)
      {
      out_mem[i] = eT(0);
      }
    
    }
  }



//! @}
