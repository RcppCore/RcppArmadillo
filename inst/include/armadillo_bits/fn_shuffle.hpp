// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup fn_shuffle
//! @{

//! \brief
//! Shuffle the rows or the columns of a matrix or vector in random fashion.
//! If dim = 0, shuffle the columns (default operation).
//! If dim = 1, shuffle the rows.

template<typename T1>
inline
const Op<T1, op_shuffle>
shuffle(const Base<typename T1::elem_type,T1>& X, const uword dim = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (dim > 1), "shuffle(): dim must be 0 or 1");
  
  return Op<T1, op_shuffle>(X.get_ref(), dim, 0);
  }



//! @}
