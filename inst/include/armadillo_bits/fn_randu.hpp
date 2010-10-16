// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_randu
//! @{


inline
double
randu()
  {
  return double(eop_aux_randu<double>());
  }


template<typename eT>
inline
typename arma_scalar_only<eT>::result
randu()
  {
  return eT(eop_aux_randu<eT>());
  }



//! Generate a vector with all elements set to random values in the [0,1] interval (uniform distribution)
arma_inline
const eOp<colvec, eop_randu>
randu(const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  return eOp<colvec, eop_randu>(n_elem, 1);
  }



template<typename vec_type>
arma_inline
const eOp<vec_type, eop_randu>
randu(const u32 n_elem, const typename arma_Mat_Col_Row_only<vec_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  if(is_Row<vec_type>::value == true)
    {
    return eOp<vec_type, eop_randu>(1, n_elem);
    }
  else
    {
    return eOp<vec_type, eop_randu>(n_elem, 1);
    }
  }



//! Generate a dense matrix with all elements set to random values in the [0,1] interval (uniform distribution)
arma_inline
const eOp<mat, eop_randu>
randu(const u32 n_rows, const u32 n_cols)
  {
  arma_extra_debug_sigprint();
  
  return eOp<mat, eop_randu>(n_rows, n_cols);
  }



template<typename mat_type>
arma_inline
const eOp<mat_type, eop_randu>
randu(const u32 n_rows, const u32 n_cols, const typename arma_Mat_Col_Row_only<mat_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOp<mat_type, eop_randu>(n_rows, n_cols);
  }



arma_inline
const eOpCube<cube, eop_randu>
randu(const u32 n_rows, const u32 n_cols, const u32 n_slices)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<cube, eop_randu>(n_rows, n_cols, n_slices);
  }



template<typename cube_type>
arma_inline
const eOpCube<cube_type, eop_randu>
randu(const u32 n_rows, const u32 n_cols, const u32 n_slices, const typename arma_Cube_only<cube_type>::result* junk = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  return eOpCube<cube_type, eop_randu>(n_rows, n_cols, n_slices);
  }



//! @}
