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


//! \addtogroup fn_sprandu
//! @{



//! Generate a sparse matrix with a randomly selected subset of the elements
//! set to random values in the [0,1] interval (uniform distribution)
template<typename obj_type>
inline
obj_type
sprandu
  (
  const uword  n_rows,
  const uword  n_cols,
  const double density,
  const typename arma_SpMat_SpCol_SpRow_only<obj_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(is_SpCol<obj_type>::value == true)
    {
    arma_debug_check( (n_cols != 1), "sprandu(): incompatible size" );
    }
  else
  if(is_SpRow<obj_type>::value == true)
    {
    arma_debug_check( (n_rows != 1), "sprandu(): incompatible size" );
    }
  
  obj_type out;
  
  out.sprandu(n_rows, n_cols, density);
  
  return out;
  }



inline
sp_mat
sprandu(const uword n_rows, const uword n_cols, const double density)
  {
  arma_extra_debug_sigprint();
  
  sp_mat out;
  
  out.sprandu(n_rows, n_cols, density);
  
  return out;
  }



//! Generate a sparse matrix with the non-zero values in the same locations as in the given sparse matrix X,
//! with the non-zero values set to random values in the [0,1] interval (uniform distribution)
template<typename T1>
inline
SpMat<typename T1::elem_type>
sprandu(const SpBase<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  SpMat<eT> out( X.get_ref() );
  
  eop_aux_randu<eT>::fill( access::rwp(out.values), out.n_nonzero );
  
  return out;
  }



//! @}
