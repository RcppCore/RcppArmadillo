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


//! \addtogroup op_chol
//! @{



template<typename T1>
inline
void
op_chol::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_chol>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(X.m, out);
  const Mat<eT>&     A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "chol(): given matrix is not square");
  
  const bool ok = auxlib::chol(out, A);
  
  if(ok == false)
    {
    out.reset();
    arma_print("chol(): matrix factorisation failed");
    }
  }



//! @}
