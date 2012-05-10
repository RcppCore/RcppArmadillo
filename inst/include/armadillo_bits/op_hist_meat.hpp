// Copyright (C) 2012 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_hist
//! @{



template<typename T1>
inline
void
op_hist::apply(Mat<uword>& out, const mtOp<uword, T1, op_hist>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword n_bins = X.aux_uword_a;
  
  const unwrap_check_mixed<T1> tmp(X.m, out);
  const Mat<eT>& A           = tmp.M;
  
  
        uword A_n_elem = A.n_elem;
  const eT*   A_mem    = A.memptr();
  
  eT min_val = priv::most_pos<eT>();
  eT max_val = priv::most_neg<eT>();
  
  uword i,j;
  for(i=0, j=1; j < A_n_elem; i+=2, j+=2)
    {
    const eT val_i = A_mem[i];
    const eT val_j = A_mem[j];
    
    if(min_val > val_i) { min_val = val_i; }
    if(min_val > val_j) { min_val = val_j; }
      
    if(max_val < val_i) { max_val = val_i; }
    if(max_val < val_j) { max_val = val_j; }
    }
  
  if(i < A_n_elem)
    {
    const eT val_i = A_mem[i];
    
    if(min_val > val_i) { min_val = val_i; }
    if(max_val < val_i) { max_val = val_i; }
    }
  
  if(arma_isfinite(min_val) == false) { min_val = priv::most_neg<eT>(); }
  if(arma_isfinite(max_val) == false) { max_val = priv::most_pos<eT>(); }
  
  if(n_bins >= 1)
    {
    Col<eT> c(n_bins);
    eT* c_mem = c.memptr();
    
    for(uword ii=0; ii < n_bins; ++ii)
      {
      c_mem[ii] = (0.5 + ii) / double(n_bins);   // TODO: may need to be modified for integer matrices
      }
    
    c = ((max_val - min_val) * c) + min_val;
    
    out = hist(A, c);
    }
  else
    {
    out.reset();
    }
  }



//! @}
