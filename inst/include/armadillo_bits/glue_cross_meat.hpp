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



//! \addtogroup glue_cross
//! @{



template<typename T1, typename T2>
inline
void
glue_cross::apply(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_cross>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> A_tmp(X.A);
  const unwrap<T2> B_tmp(X.B);
  
  const Mat<eT>& A = A_tmp.M;
  const Mat<eT>& B = B_tmp.M;
  
  arma_debug_check( ((A.n_elem != 3) || (B.n_elem != 3)), "cross(): input vectors must have 3 elements" );
  
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_ptr = out.memptr();
  const eT* A_ptr   = A.memptr();
  const eT* B_ptr   = B.memptr();
  
  // out_ptr[0] = A_ptr[1]*B_ptr[2] - A_ptr[2]*B_ptr[1];
  // out_ptr[1] = A_ptr[2]*B_ptr[0] - A_ptr[0]*B_ptr[2];
  // out_ptr[2] = A_ptr[0]*B_ptr[1] - A_ptr[1]*B_ptr[0];
  
  const eT ax = A_ptr[0];
  const eT ay = A_ptr[1];
  const eT az = A_ptr[2];
  
  const eT bx = B_ptr[0];
  const eT by = B_ptr[1];
  const eT bz = B_ptr[2];
  
  out_ptr[0] = ay*bz - az*by;
  out_ptr[1] = az*bx - ax*bz;
  out_ptr[2] = ax*by - ay*bx;
  }



//! @}
