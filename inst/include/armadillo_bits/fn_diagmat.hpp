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


//! \addtogroup fn_diagmat
//! @{


//! interpret a matrix or a vector as a diagonal matrix (i.e. off-diagonal entries are zero)
template<typename T1>
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value,
  const Op<T1, op_diagmat>
  >::result
diagmat(const T1& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_diagmat>(X);
  }



// TODO:
// create "op_diagmat2", to allow placement of vector onto a sub- or super- diagonal.
// op_diagmat2 is required, as other code assumes that op_diagmat indicates only the main diagonal)



//! @}
