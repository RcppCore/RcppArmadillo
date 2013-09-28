// Copyright (C) 2008-2013 Conrad Sanderson
// Copyright (C) 2008-2013 NICTA (www.nicta.com.au)
// Copyright (C) 2009 Edmund Highcock
// Copyright (C) 2011 Stanislav Funiak
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_eig
//! @{


//
// symmetric/hermitian matrices
//


//! Eigenvalues of real/complex symmetric/hermitian matrix X
template<typename T1>
inline
bool
eig_sym
  (
         Col<typename T1::pod_type>&     eigval,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  // unwrap_check not used as T1::elem_type and T1::pod_type may not be the same.
  // furthermore, it doesn't matter if X is an alias of eigval, as auxlib::eig_sym() makes a copy of X
  
  const bool status = auxlib::eig_sym(eigval, X);
  
  if(status == false)
    {
    eigval.reset();
    arma_bad("eig_sym(): failed to converge", false);
    }
  
  return status;
  }



//! Eigenvalues of real/complex symmetric/hermitian matrix X
template<typename T1>
inline
Col<typename T1::pod_type>
eig_sym
  (
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  Col<typename T1::pod_type> out;
  const bool status = auxlib::eig_sym(out, X);

  if(status == false)
    {
    out.reset();
    arma_bad("eig_sym(): failed to converge");
    }
  
  return out;
  }



//! Eigenvalues and eigenvectors of real/complex symmetric/hermitian matrix X
template<typename T1> 
inline
bool
eig_sym
  (
         Col<typename T1::pod_type>&     eigval,
         Mat<typename T1::elem_type>&    eigvec,
  const Base<typename T1::elem_type,T1>& X,
  const char* method =                   "",
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check( void_ptr(&eigval) == void_ptr(&eigvec), "eig_sym(): eigval is an alias of eigvec" );
  
  bool use_divide_and_conquer = false;
  
  const char sig = method[0];
  
  switch(sig)
    {
    case '\0':
    case 's':
      break;
      
    case 'd':
      use_divide_and_conquer = true;
      break;
    
    default:
      {
      arma_stop("eig_sym(): unknown method specified");
      return false;
      }
    }
  
  const bool status = (use_divide_and_conquer == false) ? auxlib::eig_sym(eigval, eigvec, X) : auxlib::eig_sym_dc(eigval, eigvec, X);
  
  if(status == false)
    {
    eigval.reset();
    eigvec.reset();
    arma_bad("eig_sym(): failed to converge", false);
    }
  
  return status;
  }



//
// general matrices
//



//! Eigenvalues and eigenvectors (both left and right) of general real/complex square matrix X
template<typename T1>
inline
bool
eig_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigval, 
         Mat<typename T1::elem_type>&                l_eigvec,
         Mat<typename T1::elem_type>&                r_eigvec,
  const Base<typename T1::elem_type,T1>&             X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check
    (
    ((&l_eigvec) == (&r_eigvec)),
    "eig_gen(): l_eigvec is an alias of r_eigvec"
    );
  
  arma_debug_check
    (
      (
      (((void*)(&eigval)) == ((void*)(&l_eigvec)))
      ||
      (((void*)(&eigval)) == ((void*)(&r_eigvec)))
      ),
    "eig_gen(): eigval is an alias of l_eigvec or r_eigvec"
    );
  
  const bool status = auxlib::eig_gen(eigval, l_eigvec, r_eigvec, X, 'b');
  
  if(status == false)
    {
    eigval.reset();
    l_eigvec.reset();
    r_eigvec.reset();
    arma_bad("eig_gen(): failed to converge", false);
    }
  
  return status;
  }



//! Eigenvalues and eigenvectors of general real square matrix X.
//! Optional argument 'side' specifies which eigenvectors should be computed:
//! 'r' for right (default) and 'l' for left.
template<typename eT, typename T1>
inline
bool
eig_gen
  (
        Col< std::complex<eT> >& eigval, 
        Mat< std::complex<eT> >& eigvec,
  const Base<eT, T1>&            X, 
  const char                     side,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  //std::cout << "real" << std::endl;
  
  arma_debug_check( ( ((void*)(&eigval)) == ((void*)(&eigvec)) ), "eig_gen(): eigval is an alias of eigvec" );
  
  Mat<eT> dummy_eigvec;
  Mat<eT> tmp_eigvec;
  
  bool status;
  
  switch(side)
    {
    case 'r':
      status = auxlib::eig_gen(eigval, dummy_eigvec, tmp_eigvec, X, side);
      break;
    
    case 'l':
      status = auxlib::eig_gen(eigval, tmp_eigvec, dummy_eigvec, X, side);
      break;
      
    default:
      arma_stop("eig_gen(): parameter 'side' is invalid");
      status = false;
    }
  
  if(status == false)
    {
    eigval.reset();
    eigvec.reset();
    arma_bad("eig_gen(): failed to converge", false);
    }
  else
    {
    const uword n = eigval.n_elem;
    
    if(n > 0)
      {
      eigvec.set_size(n,n);
      
      for(uword j=0; j<n; ++j)
        {
        if( (j < n-1) && (eigval[j] == std::conj(eigval[j+1])) )
          {
          // eigvec.col(j)   = Mat< std::complex<eT> >( tmp_eigvec.col(j),  tmp_eigvec.col(j+1) );
          // eigvec.col(j+1) = Mat< std::complex<eT> >( tmp_eigvec.col(j), -tmp_eigvec.col(j+1) );
          
          for(uword i=0; i<n; ++i)
            {
            eigvec.at(i,j)   = std::complex<eT>( tmp_eigvec.at(i,j),  tmp_eigvec.at(i,j+1) );
            eigvec.at(i,j+1) = std::complex<eT>( tmp_eigvec.at(i,j), -tmp_eigvec.at(i,j+1) );
            }
          
          ++j;
          }
        else
          {
          // eigvec.col(i) = tmp_eigvec.col(i);
          
          for(uword i=0; i<n; ++i)
            {
            eigvec.at(i,j) = std::complex<eT>(tmp_eigvec.at(i,j), eT(0));
            }
          
          }
        }
      }
    }
  
  return status;
  }



//! Eigenvalues and eigenvectors of general complex square matrix X
//! Optional argument 'side' specifies which eigenvectors should be computed:
//! 'r' for right (default) and 'l' for left.
template<typename T, typename T1>
inline
bool
eig_gen
  (
         Col<std::complex<T> >&    eigval, 
         Mat<std::complex<T> >&    eigvec,
  const Base<std::complex<T>, T1>& X, 
  const char                       side,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  //std::cout << "complex" << std::endl;
  
  arma_debug_check( ( ((void*)(&eigval)) == ((void*)(&eigvec)) ), "eig_gen(): eigval is an alias of eigvec" );
  
  Mat< std::complex<T> > dummy_eigvec;
  
  bool status;
  
  switch(side)
    {
    case 'r':
      status = auxlib::eig_gen(eigval, dummy_eigvec, eigvec, X, side);
      break;
    
    case 'l':
      status = auxlib::eig_gen(eigval, eigvec, dummy_eigvec, X, side);
      break;
      
    default:
      arma_stop("eig_gen(): parameter 'side' is invalid");
      status = false;
    }
  
  if(status == false)
    {
    eigval.reset();
    eigvec.reset();
    arma_bad("eig_gen(): failed to converge", false);
    }
  
  return status;
  }



template<typename eT, typename T1>
inline
bool
eig_gen
  (
        Col< std::complex<eT> >& eigval, 
        Mat< std::complex<eT> >& eigvec,
  const Base<eT, T1>&            X, 
  const char*                    side = "right",
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return eig_gen(eigval, eigvec, X, side[0]);
  }



template<typename T, typename T1>
inline
bool
eig_gen
  (
         Col<std::complex<T> >&    eigval, 
         Mat<std::complex<T> >&    eigvec,
  const Base<std::complex<T>, T1>& X, 
  const char*                      side = "right",
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  return eig_gen(eigval, eigvec, X, side[0]);
  }



//! @}
