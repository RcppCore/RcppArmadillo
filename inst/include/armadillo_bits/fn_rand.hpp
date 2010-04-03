// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_rand
//! @{



//! Generate a dense matrix with all elements set to random values in the [0,1] interval (uniform distribution)
arma_inline
const eOp<mat, eop_rand>
rand(const u32 n_rows, const u32 n_cols)
  {
  arma_extra_debug_sigprint();
  
  return eOp<mat, eop_rand>(n_rows, n_cols);
  }



arma_inline
const eOpCube<cube, eop_cube_rand>
rand(const u32 n_rows, const u32 n_cols, const u32 n_slices)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<cube, eop_cube_rand>(n_rows, n_cols, n_slices);
  }



template<typename mat_type>
arma_inline
const eOp<mat_type, eop_rand>
rand(const u32 n_rows, const u32 n_cols)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_Mat<mat_type>::value == false>::apply();
  
  return eOp<mat_type, eop_rand>(n_rows, n_cols);
  }



template<typename cube_type>
arma_inline
const eOpCube<cube_type, eop_cube_rand>
rand(const u32 n_rows, const u32 n_cols, const u32 n_slices)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_Cube<cube_type>::value == false>::apply();
  
  return eOpCube<cube_type, eop_cube_rand>(n_rows, n_cols, n_slices);
  }



//! Generate a vector with all elements set to random values in the [0,1] interval (uniform distribution)
arma_inline
const eOp<colvec, eop_rand>
rand(const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  return eOp<colvec, eop_rand>(n_elem, 1);
  }



template<typename vec_type>
arma_inline
const eOp<vec_type, eop_rand>
rand(const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< (is_Col<vec_type>::value == false) && (is_Row<vec_type>::value == false) >::apply();
  
  if(is_Row<vec_type>::value == true)
    {
    return eOp<vec_type, eop_rand>(1, n_elem);
    }
  else
    {
    return eOp<vec_type, eop_rand>(n_elem, 1);
    }
  }



//! @}
