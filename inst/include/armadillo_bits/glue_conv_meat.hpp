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


//! \addtogroup glue_conv
//! @{


//! rudimentary implementation of the convolution operation

template<typename T1, typename T2>
inline
void
glue_conv::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_conv>& X)
  {
  arma_extra_debug_sigprint();
  
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> A_tmp(X.A, out);
  const unwrap_check<T2> B_tmp(X.B, out);
  
  const Mat<eT>& A = A_tmp.M;
  const Mat<eT>& B = B_tmp.M;
  
  arma_debug_check( ( (A.is_vec() == false) || (B.is_vec() == false) ), "conv(): inputs must be vectors"  );
  arma_debug_check( ( (A.n_elem   == 0    ) || (B.n_elem   == 0    ) ), "conv(): zero-length input given" );
  
  const Mat<eT>& h = (A.n_elem <= B.n_elem) ? A : B;
  const Mat<eT>& x = (A.n_elem <= B.n_elem) ? B : A;
  
  
  const u32   h_n_elem = h.n_elem;
  const u32   x_n_elem = x.n_elem;
  const u32 out_n_elem = h_n_elem + x_n_elem - 1;
  
  
  (A.n_cols == 1) ? out.set_size(out_n_elem, 1) : out.set_size(1, out_n_elem);
  
  
  const eT*   h_mem = h.memptr();
  const eT*   x_mem = x.memptr();
        eT* out_mem = out.memptr();
  
  
  for(u32 out_i = 0; out_i < (h_n_elem-1); ++out_i)
    {
    eT acc = eT(0);
    
    u32 h_i = out_i;
    
    for(u32 x_i = 0; x_i <= out_i; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
    
    out_mem[out_i] = acc;
    }
  
  
  for(u32 out_i = h_n_elem-1; out_i < out_n_elem - (h_n_elem-1); ++out_i)
    {
    eT acc = eT(0);
   
    u32 h_i = h_n_elem - 1;
    
    for(u32 x_i = out_i - h_n_elem + 1; x_i <= out_i; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
      
    out_mem[out_i] = acc;
    }
  
  
  for(u32 out_i = out_n_elem - (h_n_elem-1); out_i < out_n_elem; ++out_i)
    {
    eT acc = eT(0);
    
    u32 h_i = h_n_elem - 1;
    
    for(u32 x_i = out_i - h_n_elem + 1; x_i < x_n_elem; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
    
    out_mem[out_i] = acc;
    }
  
  
  }



//! @}
