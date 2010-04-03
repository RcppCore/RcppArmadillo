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


//! \addtogroup fn_princomp_cov
//! @{



//! \brief
//! principal component analysis of a covariance matrix -- 3 arguments version
//! coeff_out     -> principal component coefficients
//! latent_out    -> principal component variances
//! explained_out -> percentage of the total variance explained by each principal component.
template<typename T1>
inline
void
princomp_cov
  (
         Mat<typename T1::elem_type>&    coeff_out,
         Col<typename T1::pod_type>&     latent_out,
         Col<typename T1::pod_type>&     explained_out,
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;

  op_princomp_cov::direct_princomp_cov(coeff_out, latent_out, explained_out, A); 
  }



//! \brief
//! principal component analysis of a covariance matrix -- 2 arguments version
//! coeff_out  -> principal component coefficients
//! latent_out -> principal component variances
template<typename T1>
inline
void
princomp_cov
  (
         Mat<typename T1::elem_type>&    coeff_out,
         Col<typename T1::pod_type>&     latent_out,
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();
 
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
 
  op_princomp_cov::direct_princomp_cov(coeff_out, latent_out, A);
  }



//! \brief
//! principal component analysis of a covariance matrix -- 1 argument version
//! coeff_out -> principal component coefficients
template<typename T1>
inline
const Op<T1, op_princomp_cov>
princomp_cov(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();

  return Op<T1, op_princomp_cov>(X.get_ref());
  }



//! @}
