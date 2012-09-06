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


//! \addtogroup spop_misc
//! @{



template<typename T1>
arma_hot
inline
void
spop_scalar_times::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  out = in.m;
  
  out *= in.aux;
  }



template<typename T1>
arma_hot
inline
void
spop_square::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_square>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out = in.m;
  
  eT* values = access::rwp(out.values);
  
  const uword n_nz = out.n_nonzero;
  
  uword i,j;
  for(i=0, j=1; j < n_nz; i+=2, j+=2)
    {
    const eT tmp_i = values[i];
    const eT tmp_j = values[j];
    
    values[i] = tmp_i*tmp_i;
    values[j] = tmp_j*tmp_j;
    }
  
  if(i < n_nz)
    {
    const eT tmp_i = values[i];
    
    values[i] = tmp_i*tmp_i;
    }
  }



//! @}
