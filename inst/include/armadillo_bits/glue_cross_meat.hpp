// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
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
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  const Proxy<T1> A(X.A);
  const Proxy<T2> B(X.B);
  
  arma_debug_check( ((A.get_n_elem() != 3) || (B.get_n_elem() != 3)), "cross(): input vectors must have 3 elements" );
  
  out.set_size(A.get_n_rows(), A.get_n_cols());
  
  eT*      out_mem = out.memptr();
  ea_type1 PA      = A.get_ea();
  ea_type2 PB      = B.get_ea();
  
  const eT ax = PA[0];
  const eT ay = PA[1];
  const eT az = PA[2];
  
  const eT bx = PB[0];
  const eT by = PB[1];
  const eT bz = PB[2];
  
  out_mem[0] = ay*bz - az*by;
  out_mem[1] = az*bx - ax*bz;
  out_mem[2] = ax*by - ay*bx;
  }



//! @}
