// Copyright (C) 2015 National ICT Australia (NICTA)
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_eig_gen
//! @{


template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, Col< std::complex<typename T1::pod_type> > >::result
eig_gen(const Base<typename T1::pod_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat<T> eigvec_l;
  Mat<T> eigvec_r;
  
  Col< std::complex<T> > eigval;
  
  const bool status = auxlib::eig_gen(eigval, eigvec_l, eigvec_r, X.get_ref(), 'n');
  
  if(status == false)
    {
    eigval.reset();
    arma_bad("eig_gen(): decomposition failed");
    }
  
  return eigval;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, Col< std::complex<typename T1::pod_type> > >::result
eig_gen(const Base< std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type     T;
  typedef typename std::complex<T> eT;
  
  Mat<eT> eigvec_l;
  Mat<eT> eigvec_r;
  
  Col<eT> eigval;
  
  const bool status = auxlib::eig_gen(eigval, eigvec_l, eigvec_r, X.get_ref(), 'n');
  
  if(status == false)
    {
    eigval.reset();
    arma_bad("eig_gen(): decomposition failed");
    }
  
  return eigval;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen(Col< std::complex<typename T1::pod_type> >& eigval, const Base<typename T1::pod_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat<T> eigvec_l;
  Mat<T> eigvec_r;
  
  const bool status = auxlib::eig_gen(eigval, eigvec_l, eigvec_r, X.get_ref(), 'n');
  
  if(status == false)
    {
    eigval.reset();
    arma_bad("eig_gen(): decomposition failed", false);
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen(Col< std::complex<typename T1::pod_type> >& eigval, const Base< std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Mat< std::complex<T> > eigvec_l;
  Mat< std::complex<T> > eigvec_r;
  
  const bool status = auxlib::eig_gen(eigval, eigvec_l, eigvec_r, X.get_ref(), 'n');
  
  if(status == false)
    {
    eigval.reset();
    arma_bad("eig_gen(): decomposition failed", false);
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen
  (
        Col< std::complex<typename T1::pod_type> >& eigval,
        Mat< std::complex<typename T1::pod_type> >& eigvec,
  const Base<typename T1::pod_type, T1>&            X,
  const char                                        mode = 'r'
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  arma_debug_check( ( void_ptr(&eigval) == void_ptr(&eigvec) ), "eig_gen(): eigval is an alias of eigvec" );
  
  arma_debug_check( ((mode != 'l') && (mode != 'r')), "eig_gen(): parameter 'mode' must be either 'l' or 'r'" );
  
  Mat<T> tmp;
  Mat<T> junk;
  
  bool status = false;
  
  switch(mode)
    {
    case 'l':  status = auxlib::eig_gen(eigval, tmp,  junk, X.get_ref(), mode);  break;
    case 'r':  status = auxlib::eig_gen(eigval, junk, tmp,  X.get_ref(), mode);  break;
    default:   ;
    }
  
  if(status == false)
    {
    eigval.reset();
    eigvec.reset();
    arma_bad("eig_gen(): decomposition failed", false);
    return false;
    }
  
  const uword N = eigval.n_elem;
  
  eigvec.set_size(N,N);
  
  if(N == 0)  { return true; }
  
  // LAPACK masochism
  
  for(uword j=0; j<N; ++j)
    {
    if( (j < N-1) && (eigval[j] == std::conj(eigval[j+1])) )
      {
      for(uword i=0; i<N; ++i)
        {
        eigvec.at(i,j)   = std::complex<T>( tmp.at(i,j),  tmp.at(i,j+1) );
        eigvec.at(i,j+1) = std::complex<T>( tmp.at(i,j), -tmp.at(i,j+1) );
        }
      
      ++j;
      }
    else
      {
      for(uword i=0; i<N; ++i)
        {
        eigvec.at(i,j) = std::complex<T>(tmp.at(i,j), T(0));
        }
      }
    }
  
  return true;
  }



template<typename T1>
inline
typename enable_if2< is_supported_blas_type<typename T1::pod_type>::value, bool >::result
eig_gen
  (
         Col< std::complex<typename T1::pod_type> >&    eigval, 
         Mat< std::complex<typename T1::pod_type> >&    eigvec,
  const Base< std::complex<typename T1::pod_type>, T1>& X, 
  const char                                            mode = 'r'
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type     T;
  typedef typename std::complex<T> eT;
  
  arma_debug_check( ( void_ptr(&eigval) == void_ptr(&eigvec) ), "eig_gen(): eigval is an alias of eigvec" );
  
  arma_debug_check( ((mode != 'l') && (mode != 'r')), "eig_gen(): parameter 'mode' must be either 'l' or 'r'" );
  
  Mat<eT> junk;
  
  bool status = false;
  
  switch(mode)
    {
    case 'l':  status = auxlib::eig_gen(eigval, eigvec, junk,   X.get_ref(), mode);  break;
    case 'r':  status = auxlib::eig_gen(eigval, junk,   eigvec, X.get_ref(), mode);  break;
    default:   ;
    }
  
  if(status == false)
    {
    eigval.reset();
    eigvec.reset();
    arma_bad("eig_gen(): decomposition failed", false);
    }
  
  return status;
  }



//! NOTE: this form is deprecated -- don't use it
template<typename T1>
arma_deprecated
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
eig_gen
  (
         Col< std::complex<typename T1::pod_type> >& eigval, 
         Mat<typename T1::elem_type>&                eigvec_l,
         Mat<typename T1::elem_type>&                eigvec_r,
  const Base<typename T1::elem_type,T1>&             X
  )
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (void_ptr(&eigvec_l) == void_ptr(&eigvec_r)), "eig_gen(): parameter 'eigvec_l' is an alias of parameter 'eigvec_r'" );
  arma_debug_check( (void_ptr(&eigvec_l) == void_ptr(&eigval)  ), "eig_gen(): parameter 'eigvec_l' is an alias of parameter 'eigval'"   );
  arma_debug_check( (void_ptr(&eigvec_r) == void_ptr(&eigval)  ), "eig_gen(): parameter 'eigvec_r' is an alias of parameter 'eigval'"   );
  
  const bool status = auxlib::eig_gen(eigval, eigvec_l, eigvec_r, X.get_ref(), 'b');
  
  if(status == false)
    {
    eigval.reset();
    eigvec_l.reset();
    eigvec_r.reset();
    arma_bad("eig_gen(): decomposition failed", false);
    }
  
  return status;
  }


//! @}
