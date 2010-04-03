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


//! \addtogroup fn_chol
//! @{



template<typename T1>
inline
bool
chol(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(X.get_ref(), out);
  const Mat<eT>&     A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "chol(): given matrix is not square");
  
  return auxlib::chol(out, A);
  }



template<typename T1>
inline
Mat<typename T1::elem_type>
chol(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> out;
  
  const bool ok = chol(out, X);
  
  if(ok == false)
    {
    arma_print("chol(): matrix factorisation failed");
    }
  
  return out;
  }



//! @}
