// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_eye
//! @{



arma_inline
const Gen<mat::elem_type, gen_ones_diag>
eye(const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  return Gen<mat::elem_type, gen_ones_diag>(n_rows, n_cols);
  }



template<typename mat_type>
arma_inline
const Gen<typename mat_type::elem_type, gen_ones_diag>
eye(const uword n_rows, const uword n_cols, const typename arma_Mat_Col_Row_only<mat_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Gen<typename mat_type::elem_type, gen_ones_diag>(n_rows, n_cols);
  }



//! @}
