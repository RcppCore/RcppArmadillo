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


//! \addtogroup fn_diagvec
//! @{


//! extract a diagonal from a matrix
template<typename T1>
arma_inline
const Op<T1, op_diagvec>
diagvec(const Base<typename T1::elem_type,T1>& X, const sword diag_id = 0)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_diagvec>(X.get_ref(), ((diag_id < 0) ? -diag_id : diag_id), ((diag_id < 0) ? 1 : 0) );
  }



//! @}
