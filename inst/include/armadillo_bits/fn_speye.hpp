// Copyright (C) 2012 Conrad Sanderson
// Copyright (C) 2012 Ryan Curtin
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_speye
//! @{



//! Generate a sparse matrix with the values along the main diagonal set to one
template<typename obj_type>
inline
obj_type
speye(const uword n_rows, const uword n_cols, const typename arma_SpMat_SpCol_SpRow_only<obj_type>::result* junk = NULL)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  obj_type out;
  
  out.eye(n_rows, n_cols);
  
  return out;
  }



// Convenience shortcut method (no template parameter necessary)
inline
sp_mat
speye(const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  sp_mat out;
  
  out.eye(n_rows, n_cols);
  
  return out;
  }



//! @}
