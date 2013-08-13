// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



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
