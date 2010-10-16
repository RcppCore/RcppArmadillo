// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_kron
//! @{



//! \brief
//! kronecker product of two matrices,
//! with the matrices having the same element type
template<typename T1, typename T2>
arma_inline
const Glue<T1,T2,glue_kron>
kron(const Base<typename T1::elem_type,T1>& A, const Base<typename T1::elem_type,T2>& B)
  {
  arma_extra_debug_sigprint();

  return Glue<T1, T2, glue_kron>(A.get_ref(), B.get_ref());
  }



//! \brief
//! kronecker product of two matrices,
//! with the matrices having different element types
template<typename T, typename T1, typename T2>
inline
Mat<typename eT_promoter<T1,T2>::eT>
kron(const Base<std::complex<T>,T1>& X, const Base<T,T2>& Y)
  {
  arma_extra_debug_sigprint();

  typedef typename std::complex<T> eT1;

  promote_type<eT1,T>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT1>& A = tmp1.M;
  const Mat<T  >& B = tmp2.M;

  Mat<eT1> out;
  
  glue_kron::direct_kron(out, A, B);
  
  return out;
  }



//! \brief
//! kronecker product of two matrices,
//! with the matrices having different element types
template<typename T, typename T1, typename T2>
inline
Mat<typename eT_promoter<T1,T2>::eT>
kron(const Base<T,T1>& X, const Base<std::complex<T>,T2>& Y)
  {
  arma_extra_debug_sigprint();

  typedef typename std::complex<T> eT2;  

  promote_type<T,eT2>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<T  >& A = tmp1.M;
  const Mat<eT2>& B = tmp2.M;

  Mat<eT2> out;
  
  glue_kron::direct_kron(out, A, B);
  
  return out;
  }



//! @}
