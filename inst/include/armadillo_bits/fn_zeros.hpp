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


//! \addtogroup fn_zeros
//! @{


//! Generate a vector with all elements set to zero
arma_inline
const Gen<vec::elem_type, gen_zeros>
zeros(const uword n_elem)
  {
  arma_extra_debug_sigprint();
  
  return Gen<vec::elem_type, gen_zeros>(n_elem, 1);
  }



template<typename vec_type>
arma_inline
const Gen<typename vec_type::elem_type, gen_zeros>
zeros(const uword n_elem, const arma_empty_class junk1 = arma_empty_class(), const typename arma_Mat_Col_Row_only<vec_type>::result* junk2 = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk1);
  arma_ignore(junk2);
  
  if(is_Row<vec_type>::value == true)
    {
    return Gen<typename vec_type::elem_type, gen_zeros>(1, n_elem);
    }
  else
    {
    return Gen<typename vec_type::elem_type, gen_zeros>(n_elem, 1);
    }
  }



//! Generate a dense matrix with all elements set to zero
arma_inline
const Gen<mat::elem_type, gen_zeros>
zeros(const uword n_rows, const uword n_cols)
  {
  arma_extra_debug_sigprint();
  
  return Gen<mat::elem_type, gen_zeros>(n_rows, n_cols);
  }



template<typename mat_type>
arma_inline
const Gen<typename mat_type::elem_type, gen_zeros>
zeros(const uword n_rows, const uword n_cols, const typename arma_Mat_Col_Row_only<mat_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return Gen<typename mat_type::elem_type, gen_zeros>(n_rows, n_cols);
  }



arma_inline
const GenCube<cube::elem_type, gen_zeros>
zeros(const uword n_rows, const uword n_cols, const uword n_slices)
  {
  arma_extra_debug_sigprint();
  
  return GenCube<cube::elem_type, gen_zeros>(n_rows, n_cols, n_slices);
  }



template<typename cube_type>
arma_inline
const GenCube<typename cube_type::elem_type, gen_zeros>
zeros(const uword n_rows, const uword n_cols, const uword n_slices, const typename arma_Cube_only<cube_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return GenCube<typename cube_type::elem_type, gen_zeros>(n_rows, n_cols, n_slices);
  }



//! @}
