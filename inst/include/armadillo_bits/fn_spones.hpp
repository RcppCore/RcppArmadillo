// Copyright (C) 2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_spones
//! @{



//! Generate a sparse matrix with the non-zero values in the same locations as in the given sparse matrix X,
//! with the non-zero values set to one
template<typename T1>
inline
SpMat<typename T1::elem_type>
spones(const SpBase<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  SpMat<eT> out( X.get_ref() );
  
  const uword nnz = out.n_nonzero;
  
  eT* values = access::rwp(out.values);
  
  for(uword i=0; i < nnz; ++i)
    {
    values[i] = eT(1);
    }
  
  return out;
  }



//! @}
