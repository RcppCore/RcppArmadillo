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


//! \addtogroup fn_diagmat
//! @{


//! interpret a matrix or a vector as a diagonal matrix (i.e. off-diagonal entries are zero)
template<typename T1>
arma_inline
const Op<T1, op_diagmat>
diagmat(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_diagmat>(X.get_ref());
  }



// TODO:
// create "op_diagmat2", to allow placement of vector onto a sub- or super- diagonal.
// op_diagmat2 is required, as other code assumes that op_diagmat indicates only the main diagonal)



//! @}
