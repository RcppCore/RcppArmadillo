// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_flip
//! @{



template<typename T1>
inline
void
op_flipud::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_flipud>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>  tmp(in.m);
  const Mat<eT> X = tmp.M;
  
  if(&out != &X)
    {
    out.copy_size(X);
    
    for(u32 i=0; i<X.n_rows; ++i)
      {
      out.row(i) = X.row(X.n_rows-1 - i);
      }
    }
  else
    {
    const u32 N = X.n_rows / 2;
    
    for(u32 i=0; i<N; ++i)
      {
      out.swap_rows(i, X.n_rows-1 - i);
      }
    }
  }



template<typename T1>
inline
void
op_fliplr::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_fliplr>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>  tmp(in.m);
  const Mat<eT> X = tmp.M;
  
  if(&out != &X)
    {
    out.copy_size(X);
    
    for(u32 i=0; i<X.n_cols; ++i)
      {
      out.col(i) = X.col(X.n_cols-1 - i);
      }
    }
  else
    {
    const u32 N = X.n_cols / 2;
    
    for(u32 i=0; i<N; ++i)
      {
      out.swap_cols(i, X.n_cols-1 - i);
      }
    }
  }



//! @}
