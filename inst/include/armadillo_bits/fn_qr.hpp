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


//! \addtogroup fn_qr
//! @{


//! QR decomposition
template<typename T1>
inline
void
qr(Mat<typename T1::elem_type>& Q, Mat<typename T1::elem_type>& R, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;

  arma_debug_check( (&Q == &R), "qr(): Q and R are the same object");
  
  const unwrap_check<T1> tmp1(X.get_ref(), Q);
  const Mat<eT>&     A = tmp1.M;
  
  const unwrap_check< Mat<eT> > tmp2(A, R);
  const Mat<eT>&            B = tmp2.M;
  
  const bool ok = auxlib::qr(Q, R, B);
  
  if(ok == false)
    {
    arma_print("qr(): matrix factorisation failed");
    }
  
  }


//! @}
