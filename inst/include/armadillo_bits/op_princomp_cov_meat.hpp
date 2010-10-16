// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// Copyright (C) 2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_princomp_cov
//! @{



//! \brief
//! principal component analysis of a covariance matrix -- 3 arguments version
//! computation is done via singular value decomposition
//! coeff_out     -> principal component coefficients
//! latent_out    -> principal component variances
//! explained_out -> percentage of the total variance explained by each principal component
template<typename eT>
inline
void
op_princomp_cov::direct_princomp_cov
  (
        Mat<eT>& coeff_out,
        Col<eT>& latent_out,
        Col<eT>& explained_out,
  const Mat<eT>& in
  )
  {
  arma_extra_debug_sigprint();
  
  // computation of the covariance matrix
  const Mat<eT> in_cov = cov(in);
  
  // singular value decomposition
  Mat<eT> U;

  const bool svd_ok = svd(U, latent_out, coeff_out, in_cov);
    
  if(svd_ok == false)
    {
    arma_print("princomp_cov(): singular value decomposition failed");
      
    coeff_out.reset();
    latent_out.reset();
    explained_out.reset();
      
    return;
    }
  
  explained_out =  (eT(100) * latent_out) / sum(latent_out);
  }



//! \brief
//! principal component analysis of a covariance matrix -- 2 arguments version
//! computation is done via singular value decomposition
//! coeff_out     -> principal component coefficients
//! latent_out    -> principal component variances
template<typename eT>
inline
void
op_princomp_cov::direct_princomp_cov
  (
        Mat<eT>& coeff_out,
        Col<eT>& latent_out,
  const Mat<eT>& in
  )
  {
  arma_extra_debug_sigprint();
  
  // computation of the covariance matrix
  const Mat<eT> in_cov = cov(in);
  
  // singular value decomposition
  Mat<eT> U;

  const bool svd_ok = svd(U, latent_out, coeff_out, in_cov);
    
  if(svd_ok == false)
    {
    arma_print("princomp_cov(): singular value decomposition failed");
      
    coeff_out.reset();
    latent_out.reset();
      
    return;
    }

  }



//! \brief
//! principal component analysis of a covariance matrix -- 1 argument version
//! computation is done via singular value decomposition
//! coeff_out     -> principal component coefficients
template<typename eT>
inline
void
op_princomp_cov::direct_princomp_cov
  (
        Mat<eT>& coeff_out,
  const Mat<eT>& in
  )
  {
  arma_extra_debug_sigprint();
  
  // computation of the covariance matrix
  const Mat<eT> in_cov = cov(in);
  
  // singular value decomposition
  Mat<eT> U;
  Col<eT> s;

  const bool svd_ok = svd(U, s, coeff_out, in_cov);
    
  if(svd_ok == false)
    {
    arma_print("princomp_cov(): singular value decomposition failed");
      
    coeff_out.reset();
     
    return;
    }
  
  }



//! \brief
//! principal component analysis of a covariance matrix -- 3 arguments complex version
//! computation is done via singular value decomposition
//! coeff_out     -> principal component coefficients
//! latent_out    -> principal component variances
//! explained_out -> percentage of the total variance explained by each principal component
template<typename T>
inline
void
op_princomp_cov::direct_princomp_cov
  (
        Mat< std::complex<T> >& coeff_out,
                        Col<T>& latent_out,
                        Col<T>& explained_out,
  const Mat< std::complex<T> >& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  // computation of the covariance matrix
  const Mat<eT> in_cov = cov(in);
  
  // singular value decomposition
  Mat<eT> U;

  const bool svd_ok = svd(U, latent_out, coeff_out, in_cov);
    
  if(svd_ok == false)
    {
    arma_print("princomp_cov(): singular value decomposition failed");
      
    coeff_out.reset();
    latent_out.reset();
    explained_out.reset();
      
    return;
    }
  
  explained_out =  (T(100) * latent_out) / sum(latent_out);
  }



//! \brief
//! principal component analysis of a covariance matrix -- 2 arguments complex version
//! computation is done via singular value decomposition
//! coeff_out     -> principal component coefficients
//! latent_out    -> principal component variances
template<typename T>
inline
void
op_princomp_cov::direct_princomp_cov
  (
        Mat< std::complex<T> >& coeff_out,
                        Col<T>& latent_out,
  const Mat< std::complex<T> >& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  // computation of the covariance matrix
  const Mat<eT> in_cov = cov(in);
  
  // singular value decomposition
  Mat<eT> U;

  const bool svd_ok = svd(U, latent_out, coeff_out, in_cov);
    
  if(svd_ok == false)
    {
    arma_print("princomp_cov(): singular value decomposition failed");
      
    coeff_out.reset();
    latent_out.reset();
      
    return;
    }

  }



//! \brief
//! principal component analysis of a covariance matrix -- 1 argument complex version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
template<typename T>
inline
void
op_princomp_cov::direct_princomp_cov
  (
        Mat< std::complex<T> >& coeff_out,
  const Mat< std::complex<T> >& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  // computation of the covariance matrix
  const Mat<eT> in_cov = cov(in);
  
  // singular value decomposition
  Mat<eT> U;
  Col<T> s;

  const bool svd_ok = svd(U, s, coeff_out, in_cov);
    
  if(svd_ok == false)
    {
    arma_print("princomp_cov(): singular value decomposition failed");
      
    coeff_out.reset();
      
    return;
    }
  
  }



template<typename T1>
inline
void
op_princomp_cov::apply
  (
        Mat<typename T1::elem_type>& out,
  const Op<T1,op_princomp_cov>&      in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& A     = tmp.M;
  
  op_princomp_cov::direct_princomp_cov(out, A);  
  }



//! @}
