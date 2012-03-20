// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
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
const Gen<mat, gen_ones_diag>
eye(const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  return Gen<mat, gen_ones_diag>(n_rows, n_cols);
  }



template<typename obj_type>
arma_inline
const Gen<obj_type, gen_ones_diag>
eye(const uword n_rows, const uword n_cols, const typename arma_Mat_Col_Row_only<obj_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  if(is_Col<obj_type>::value == true)
    {
    arma_debug_check( (n_cols != 1), "eye(): incompatible size" );
    }
  else
  if(is_Row<obj_type>::value == true)
    {
    arma_debug_check( (n_rows != 1), "eye(): incompatible size" );
    }
  
  return Gen<obj_type, gen_ones_diag>(n_rows, n_cols);
  }



//! @}
