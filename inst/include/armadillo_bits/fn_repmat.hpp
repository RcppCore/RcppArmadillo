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



//! \addtogroup fn_repmat
//! @{


//! \brief
//! delayed 'repeat matrix' construction of a matrix
template<typename T1>
arma_inline
const Op<T1, op_repmat>
repmat(const Base<typename T1::elem_type,T1>& A, const uword r, const uword c)
  {
  arma_extra_debug_sigprint();

  return Op<T1, op_repmat>(A.get_ref(), r, c);
  }



//! @}
