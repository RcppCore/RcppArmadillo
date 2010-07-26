// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
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


//! \addtogroup fn_princomp
//! @{



//! \brief
//! principal component analysis -- 4 arguments version
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
//! latent_out   -> eigenvalues of principal vectors
//! tsquared_out -> Hotelling's T^2 statistic
template<typename T1>
inline
void
princomp
  (
         Mat<typename T1::elem_type>&    coeff_out,
         Mat<typename T1::elem_type>&    score_out,
         Col<typename T1::pod_type>&     latent_out,
         Col<typename T1::elem_type>&    tsquared_out,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  op_princomp::direct_princomp(coeff_out, score_out, latent_out, tsquared_out, A);
  }



//! \brief
//! principal component analysis -- 3 arguments version
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
//! latent_out   -> eigenvalues of principal vectors
template<typename T1>
inline
void
princomp
  (
         Mat<typename T1::elem_type>&    coeff_out,
         Mat<typename T1::elem_type>&    score_out,
         Col<typename T1::pod_type>&     latent_out,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  op_princomp::direct_princomp(coeff_out, score_out, latent_out, A); 
  }



//! \brief
//! principal component analysis -- 2 arguments version
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
template<typename T1>
inline
void
princomp
  (
         Mat<typename T1::elem_type>&    coeff_out,
         Mat<typename T1::elem_type>&    score_out,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  op_princomp::direct_princomp(coeff_out, score_out, A); 
  }



//! \brief
//! principal component analysis -- 1 argument version
//! coeff_out    -> principal component coefficients
template<typename T1>
inline
const Op<T1, op_princomp>
princomp
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();

  return Op<T1, op_princomp>(X.get_ref());
  }



//! @}
