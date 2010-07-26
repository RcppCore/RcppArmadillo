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


//! \addtogroup fn_lu
//! @{


//! immediate lower upper decomposition
template<typename T1>
inline
void
lu
  (
        Mat<typename T1::elem_type>&     L,
        Mat<typename T1::elem_type>&     U,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( (&L == &U), "lu(): L and U are the same object");
  
  const unwrap_check<T1> tmp1(X.get_ref(), L);
  const Mat<eT>&     A = tmp1.M;
  
  const unwrap_check< Mat<eT> > tmp2(A, U);
  const Mat<eT>&            B = tmp2.M;
  
  auxlib::lu(L, U, B);
  }



//! immediate lower upper decomposition, also providing the permutation matrix
template<typename T1>
inline
void
lu
  (
        Mat<typename T1::elem_type>&     L,
        Mat<typename T1::elem_type>&     U, 
        Mat<typename T1::elem_type>&     P,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_debug_check( ( (&L == &U) || (&L == &P) || (&U == &P) ), "lu(): two or more output objects are the same object");

  const unwrap_check<T1> tmp1(X.get_ref(), L);
  const Mat<eT>&     A = tmp1.M;
  
  const unwrap_check< Mat<eT> > tmp2(A, U);
  const Mat<eT>&            B = tmp2.M;
  
  const unwrap_check< Mat<eT> > tmp3(B, P);
  const Mat<eT>&            C = tmp3.M;
  
  auxlib::lu(L, U, P, C);
  }


//! @}
