// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// Copyright (C) 2010 Dimitrios Bouzas
// Copyright (C) 2011 Stanislav Funiak
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_princomp
//! @{



//! \brief
//! principal component analysis -- 4 arguments version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
//! latent_out   -> eigenvalues of principal vectors
//! tsquared_out -> Hotelling's T^2 statistic
template<typename eT>
inline
bool
op_princomp::direct_princomp
  (
        Mat<eT>& coeff_out,
        Mat<eT>& score_out,
        Col<eT>& latent_out, 
        Col<eT>& tsquared_out,
  const Mat<eT>& in
  )
  {
  arma_extra_debug_sigprint();

  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows > 1) // more than one sample
    {
    // subtract the mean - use score_out as temporary matrix
    score_out = in - repmat(mean(in), n_rows, 1);
 	  
    // singular value decomposition
    Mat<eT> U;
    Col<eT> s;
    
    const bool svd_ok = svd(U,s,coeff_out,score_out);
    
    if(svd_ok == false)
      {
      return false;
      }
    
    
    //U.reset();  // TODO: do we need this ?  U will get automatically deleted anyway
    
    // normalize the eigenvalues
    s /= std::sqrt( double(n_rows - 1) );
    
    // project the samples to the principals
    score_out *= coeff_out;
    
    if(n_rows <= n_cols) // number of samples is less than their dimensionality
      {
      score_out.cols(n_rows-1,n_cols-1).zeros();
      
      //Col<eT> s_tmp = zeros< Col<eT> >(n_cols);
      Col<eT> s_tmp(n_cols);
      s_tmp.zeros();
      
      s_tmp.rows(0,n_rows-2) = s.rows(0,n_rows-2);
      s = s_tmp;
          
      // compute the Hotelling's T-squared
      s_tmp.rows(0,n_rows-2) = eT(1) / s_tmp.rows(0,n_rows-2);
      
      const Mat<eT> S = score_out * diagmat(Col<eT>(s_tmp));   
      tsquared_out = sum(S%S,1); 
      }
    else
      {
      // compute the Hotelling's T-squared   
      const Mat<eT> S = score_out * diagmat(Col<eT>( eT(1) / s));
      tsquared_out = sum(S%S,1);
      }
            
    // compute the eigenvalues of the principal vectors
    latent_out = s%s;
    }
  else // 0 or 1 samples
    {
    coeff_out.eye(n_cols, n_cols);
    
    score_out.copy_size(in);
    score_out.zeros();
    
    latent_out.set_size(n_cols);
    latent_out.zeros();
    
    tsquared_out.set_size(n_rows);
    tsquared_out.zeros();
    }
  
  return true;
  }



//! \brief
//! principal component analysis -- 3 arguments version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
//! latent_out   -> eigenvalues of principal vectors
template<typename eT>
inline
bool
op_princomp::direct_princomp
  (
        Mat<eT>& coeff_out,
        Mat<eT>& score_out,
        Col<eT>& latent_out,
  const Mat<eT>& in
  )
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows > 1) // more than one sample
    {
    // subtract the mean - use score_out as temporary matrix
    score_out = in - repmat(mean(in), n_rows, 1);
 	  
    // singular value decomposition
    Mat<eT> U;
    Col<eT> s;
    
    const bool svd_ok = svd(U,s,coeff_out,score_out);
    
    if(svd_ok == false)
      {
      return false;
      }
    
    
    // U.reset();
    
    // normalize the eigenvalues
    s /= std::sqrt( double(n_rows - 1) );
    
    // project the samples to the principals
    score_out *= coeff_out;
    
    if(n_rows <= n_cols) // number of samples is less than their dimensionality
      {
      score_out.cols(n_rows-1,n_cols-1).zeros();
      
      Col<eT> s_tmp = zeros< Col<eT> >(n_cols);
      s_tmp.rows(0,n_rows-2) = s.rows(0,n_rows-2);
      s = s_tmp;
      }
      
    // compute the eigenvalues of the principal vectors
    latent_out = s%s;
    
    }
  else // 0 or 1 samples
    {
    coeff_out.eye(n_cols, n_cols);
    
    score_out.copy_size(in);
    score_out.zeros();
    
    latent_out.set_size(n_cols);
    latent_out.zeros(); 
    }
  
  return true;
  }



//! \brief
//! principal component analysis -- 2 arguments version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
template<typename eT>
inline
bool
op_princomp::direct_princomp
  (
        Mat<eT>& coeff_out,
        Mat<eT>& score_out,
  const Mat<eT>& in
  )
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows > 1) // more than one sample
    {
    // subtract the mean - use score_out as temporary matrix
    score_out = in - repmat(mean(in), n_rows, 1);
 	  
    // singular value decomposition
    Mat<eT> U;
    Col<eT> s;
    
    const bool svd_ok = svd(U,s,coeff_out,score_out);
    
    if(svd_ok == false)
      {
      return false;
      }
    
    // U.reset();
    
    // normalize the eigenvalues
    s /= std::sqrt( double(n_rows - 1) );
    
    // project the samples to the principals
    score_out *= coeff_out;
    
    if(n_rows <= n_cols) // number of samples is less than their dimensionality
      {
      score_out.cols(n_rows-1,n_cols-1).zeros();
      
      Col<eT> s_tmp = zeros< Col<eT> >(n_cols);
      s_tmp.rows(0,n_rows-2) = s.rows(0,n_rows-2);
      s = s_tmp;
      }
    }
  else // 0 or 1 samples
    {
    coeff_out.eye(n_cols, n_cols);
    score_out.copy_size(in);
    score_out.zeros();
    }
  
  return true;
  }



//! \brief
//! principal component analysis -- 1 argument version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
template<typename eT>
inline
bool
op_princomp::direct_princomp
  (
        Mat<eT>& coeff_out,
  const Mat<eT>& in
  )
  {
  arma_extra_debug_sigprint();
  
  if(in.n_elem != 0)
    {
    // singular value decomposition
    Mat<eT> U;
    Col<eT> s;
    
    const Mat<eT> tmp = in - repmat(mean(in), in.n_rows, 1);
    
    const bool svd_ok = svd(U,s,coeff_out, tmp);
    
    if(svd_ok == false)
      {
      return false;
      }
    }
  else
    {
    coeff_out.eye(in.n_cols, in.n_cols);
    }
  
  return true;
  }



//! \brief
//! principal component analysis -- 4 arguments complex version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
//! latent_out   -> eigenvalues of principal vectors
//! tsquared_out -> Hotelling's T^2 statistic
template<typename T>
inline
bool
op_princomp::direct_princomp
  (
        Mat< std::complex<T> >& coeff_out,
        Mat< std::complex<T> >& score_out,
        Col<T>&                 latent_out,
        Col< std::complex<T> >& tsquared_out,
  const Mat< std::complex<T> >& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows > 1) // more than one sample
    {
    // subtract the mean - use score_out as temporary matrix
    score_out = in - repmat(mean(in), n_rows, 1);
 	  
    // singular value decomposition
    Mat<eT> U;
    Col<T> s;
    
    const bool svd_ok = svd(U,s,coeff_out,score_out); 
    
    if(svd_ok == false)
      {
      return false;
      }
    
    
    //U.reset();
    
    // normalize the eigenvalues
    s /= std::sqrt( double(n_rows - 1) );
    
    // project the samples to the principals
    score_out *= coeff_out;
    
    if(n_rows <= n_cols) // number of samples is less than their dimensionality
      {
      score_out.cols(n_rows-1,n_cols-1).zeros();
      
      Col<T> s_tmp = zeros< Col<T> >(n_cols);
      s_tmp.rows(0,n_rows-2) = s.rows(0,n_rows-2);
      s = s_tmp;
          
      // compute the Hotelling's T-squared   
      s_tmp.rows(0,n_rows-2) = 1.0 / s_tmp.rows(0,n_rows-2);
      const Mat<eT> S = score_out * diagmat(Col<T>(s_tmp));                     
      tsquared_out = sum(S%S,1); 
      }
    else
      {
      // compute the Hotelling's T-squared   
      const Mat<eT> S = score_out * diagmat(Col<T>(T(1) / s));                     
      tsquared_out = sum(S%S,1);
      }
    
    // compute the eigenvalues of the principal vectors
    latent_out = s%s;
    
    }
  else // 0 or 1 samples
    {
    coeff_out.eye(n_cols, n_cols);
    
    score_out.copy_size(in);
    score_out.zeros();
      
    latent_out.set_size(n_cols);
    latent_out.zeros();
      
    tsquared_out.set_size(n_rows);
    tsquared_out.zeros();
    }
  
  return true;
  }



//! \brief
//! principal component analysis -- 3 arguments complex version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
//! latent_out   -> eigenvalues of principal vectors
template<typename T>
inline
bool
op_princomp::direct_princomp
  (
        Mat< std::complex<T> >& coeff_out,
        Mat< std::complex<T> >& score_out,
        Col<T>&                 latent_out,
  const Mat< std::complex<T> >& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows > 1) // more than one sample
    {
    // subtract the mean - use score_out as temporary matrix
    score_out = in - repmat(mean(in), n_rows, 1);
 	  
    // singular value decomposition
    Mat<eT> U;
    Col< T> s;
    
    const bool svd_ok = svd(U,s,coeff_out,score_out);
    
    if(svd_ok == false)
      {
      return false;
      }
    
    
    // U.reset();
    
    // normalize the eigenvalues
    s /= std::sqrt( double(n_rows - 1) );
    
    // project the samples to the principals
    score_out *= coeff_out;
    
    if(n_rows <= n_cols) // number of samples is less than their dimensionality
      {
      score_out.cols(n_rows-1,n_cols-1).zeros();
      
      Col<T> s_tmp = zeros< Col<T> >(n_cols);
      s_tmp.rows(0,n_rows-2) = s.rows(0,n_rows-2);
      s = s_tmp;
      }
      
    // compute the eigenvalues of the principal vectors
    latent_out = s%s;

    }
  else // 0 or 1 samples
    {
    coeff_out.eye(n_cols, n_cols);

    score_out.copy_size(in);
    score_out.zeros();

    latent_out.set_size(n_cols);
    latent_out.zeros();
    }
  
  return true;
  }



//! \brief
//! principal component analysis -- 2 arguments complex version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
//! score_out    -> projected samples
template<typename T>
inline
bool
op_princomp::direct_princomp
  (
        Mat< std::complex<T> >& coeff_out,
        Mat< std::complex<T> >& score_out,
  const Mat< std::complex<T> >& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;
  
  const uword n_rows = in.n_rows;
  const uword n_cols = in.n_cols;
  
  if(n_rows > 1) // more than one sample
    {
    // subtract the mean - use score_out as temporary matrix
    score_out = in - repmat(mean(in), n_rows, 1);
 	  
    // singular value decomposition
    Mat<eT> U;
    Col< T> s;
    
    const bool svd_ok = svd(U,s,coeff_out,score_out);
    
    if(svd_ok == false)
      {
      return false;
      }
    
    // U.reset();
    
    // normalize the eigenvalues
    s /= std::sqrt( double(n_rows - 1) );

    // project the samples to the principals
    score_out *= coeff_out;

    if(n_rows <= n_cols) // number of samples is less than their dimensionality
      {
      score_out.cols(n_rows-1,n_cols-1).zeros();
      }

    }
  else // 0 or 1 samples
    {
    coeff_out.eye(n_cols, n_cols);
    
    score_out.copy_size(in);
    score_out.zeros();
    }
  
  return true;
  }



//! \brief
//! principal component analysis -- 1 argument complex version
//! computation is done via singular value decomposition
//! coeff_out    -> principal component coefficients
template<typename T>
inline
bool
op_princomp::direct_princomp
  (
        Mat< std::complex<T> >& coeff_out,
  const Mat< std::complex<T> >& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  if(in.n_elem != 0)
    {
 	  // singular value decomposition
 	  Mat<eT> U;
    Col< T> s;
    
    const Mat<eT> tmp = in - repmat(mean(in), in.n_rows, 1);
    
    const bool svd_ok = svd(U,s,coeff_out, tmp);
    
    if(svd_ok == false)
      {
      return false;
      }
    }
  else
    {
    coeff_out.eye(in.n_cols, in.n_cols);
    }
  
  return true;
  }



template<typename T1>
inline
void
op_princomp::apply
  (
        Mat<typename T1::elem_type>& out,
  const Op<T1,op_princomp>&          in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& A     = tmp.M;
  
  const bool status = op_princomp::direct_princomp(out, A);
  
  if(status == false)
    {
    out.reset();
    
    arma_bad("princomp(): failed to converge");
    }
  }



//! @}
