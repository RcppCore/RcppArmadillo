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


//! \addtogroup op_diagvec
//! @{



template<typename T1>
inline
void
op_diagvec::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagvec>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const sword id = (X.aux_uword_b > 0) ? -sword(X.aux_uword_a) : sword(X.aux_uword_a);
  
  const unwrap_check<T1> tmp(X.m, out);
  const Mat<eT>& A     = tmp.M;

  out = A.diag(id);
  }



//! @}
