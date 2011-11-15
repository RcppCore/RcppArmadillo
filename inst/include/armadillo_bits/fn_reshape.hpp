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


//! \addtogroup fn_reshape
//! @{



template<typename T1>
inline
const Op<T1, op_reshape>
reshape(const Base<typename T1::elem_type,T1>& X, const uword in_n_rows, const uword in_n_cols, const uword dim = 0)
  {
  arma_extra_debug_sigprint();

  arma_debug_check( (dim > 1), "reshape(): dim must be 0 or 1");

  typedef typename T1::elem_type eT;
  
  return Op<T1, op_reshape>(X.get_ref(), in_n_rows, in_n_cols, dim, 'j');
  }



template<typename T1>
inline
const OpCube<T1, op_reshape>
reshape(const BaseCube<typename T1::elem_type,T1>& X, const uword in_n_rows, const uword in_n_cols, const uword in_n_slices, const uword dim = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (dim > 1), "reshape(): dim must be 0 or 1");

  typedef typename T1::elem_type eT;
  
  return OpCube<T1, op_reshape>(X.get_ref(), in_n_rows, in_n_cols, in_n_slices, dim, 'j');
  }



//! @}
