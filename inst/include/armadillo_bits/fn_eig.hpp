// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Edmund Highcock (edmund dot highcock at merton dot ox dot ac dot uk)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_eig
//! @{


//
// symmetric/hermitian matrices
//


//! Eigenvalues of real/complex symmetric/hermitian matrix X
template<typename T1>
inline
void
eig_sym(Col<typename T1::pod_type>& eigval, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // unwrap_check not used as T1::elem_type and T1::pod_type may not be the same.
  // furthermore, it doesn't matter if A is an alias of S, as auxlib::eig() makes a copy of A

  const unwrap<T1> tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;

  auxlib::eig_sym(eigval, A);
  }



//! Eigenvalues of real/complex symmetric/hermitian matrix X
template<typename T1>
inline
Col<typename T1::pod_type>
eig_sym(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Col<typename T1::pod_type> out;
  eig_sym(out, X);
  
  return out;
  }


//! Eigenvalues and eigenvectors of real/complex symmetric/hermitian matrix X
template<typename T1> 
inline
void
eig_sym
  (
  Col<typename T1::pod_type>& eigval,
  Mat<typename T1::elem_type>& eigvec,
  const Base<typename T1::elem_type,T1>& X
  )
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  auxlib::eig_sym(eigval, eigvec, A);
  }



//
// general matrices
//



//! Eigenvalues and eigenvectors (both left and right) of general real/complex square matrix X
template<typename T1>
inline
void
eig_gen
  (
  Col< std::complex<typename T1::pod_type> >& eigval, 
  Mat<typename T1::elem_type>&                l_eigvec,
  Mat<typename T1::elem_type>&                r_eigvec,
  const Base<typename T1::elem_type,T1>&      X
  )
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;

  auxlib::eig_gen(eigval, l_eigvec, r_eigvec, A, 'b');
  }



//! Eigenvalues and eigenvectors of general real square matrix X.
//! Optional argument 'side' specifies which eigenvectors should be computed:
//! 'r' for right (default) and 'l' for left.
template<typename eT, typename T1>
inline
void
eig_gen
  (
  Col< std::complex<eT> >& eigval, 
  Mat< std::complex<eT> >& eigvec,
  const Base<eT, T1>& X, 
  const char side = 'r'
  )
  {
  arma_extra_debug_sigprint();

  //std::cout << "real" << std::endl;

  const unwrap<T1> tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;

  Mat<eT> dummy_eigvec;
  Mat<eT> tmp_eigvec;
  
  switch(side)
    {
    case 'r':
      auxlib::eig_gen(eigval, dummy_eigvec, tmp_eigvec, A, side);
      break;

    case 'l':
      auxlib::eig_gen(eigval, tmp_eigvec, dummy_eigvec, A, side);
      break;
      
    default:
      arma_stop("eig_gen(): parameter 'side' is invalid");
    }


  const u32 n = A.n_rows;

  if(n > 0)
    {
    eigvec.set_size(n,n);

    for(u32 j=0; j<n; ++j)
      {
      if( (j < n-1) && (eigval[j] == std::conj(eigval[j+1])) )
        {
        // eigvec.col(j)   = Mat< std::complex<eT> >( tmp_eigvec.col(j),  tmp_eigvec.col(j+1) );
        // eigvec.col(j+1) = Mat< std::complex<eT> >( tmp_eigvec.col(j), -tmp_eigvec.col(j+1) );

        for(u32 i=0; i<n; ++i)
          {
          eigvec.at(i,j)   = std::complex<eT>( tmp_eigvec.at(i,j),  tmp_eigvec.at(i,j+1) );
          eigvec.at(i,j+1) = std::complex<eT>( tmp_eigvec.at(i,j), -tmp_eigvec.at(i,j+1) );
          }

        ++j;
        }
      else
        {
        // eigvec.col(i) = tmp_eigvec.col(i);

        for(u32 i=0; i<n; ++i)
          {
          eigvec.at(i,j) = std::complex<eT>(tmp_eigvec.at(i,j), eT(0));
          }

        }
      }
    }

  }



//! Eigenvalues and eigenvectors of general complex square matrix X
//! Optional argument 'side' specifies which eigenvectors should be computed:
//! 'r' for right (default) and 'l' for left.
template<typename T, typename T1>
inline
void
eig_gen
  (
  Col< std::complex<T> >& eigval, 
  Mat< std::complex<T> >& eigvec,
  const Base<std::complex<T>, T1>& X, 
  const char side = 'r'
  )
  {
  arma_extra_debug_sigprint();
  //std::cout << "complex" << std::endl;

  typedef typename std::complex<T> eT;

  const unwrap<T1> tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;

  Mat<eT> dummy_eigvec;
  
  switch(side)
    {
    case 'r':
      auxlib::eig_gen(eigval, dummy_eigvec, eigvec, A, side);
      break;

    case 'l':
      auxlib::eig_gen(eigval, eigvec, dummy_eigvec, A, side);
      break;
      
    default:
      arma_stop("eig_gen(): parameter 'side' is invalid");
    }
  }



//! @}

