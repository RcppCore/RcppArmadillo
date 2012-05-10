// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_dot
//! @{




//! for two arrays, generic version
template<typename eT>
arma_hot
arma_pure
arma_inline
eT
op_dot::direct_dot_arma(const uword n_elem, const eT* const A, const eT* const B)
  {
  arma_extra_debug_sigprint();
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  uword i, j;
  
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    val1 += A[i] * B[i];
    val2 += A[j] * B[j];
    }
  
  if(i < n_elem)
    {
    val1 += A[i] * B[i];
    }
  
  return val1 + val2;
  }



//! for two arrays, float and double version
template<typename eT>
arma_hot
arma_pure
inline
typename arma_float_only<eT>::result
op_dot::direct_dot(const uword n_elem, const eT* const A, const eT* const B)
  {
  arma_extra_debug_sigprint();
  
  if( n_elem <= (128/sizeof(eT)) )
    {
    return op_dot::direct_dot_arma(n_elem, A, B);
    }
  else
    {
    #if defined(ARMA_USE_ATLAS)
      {
      arma_extra_debug_print("atlas::cblas_dot()");
      
      return atlas::cblas_dot(n_elem, A, B);
      }
    #elif defined(ARMA_USE_BLAS)
      {
      arma_extra_debug_print("blas::dot()");
      
      return blas::dot(n_elem, A, B);
      }
    #else
      {
      return op_dot::direct_dot_arma(n_elem, A, B);
      }
    #endif
    }
  }



//! for two arrays, complex version
template<typename eT>
inline
arma_hot
arma_pure
typename arma_cx_only<eT>::result
op_dot::direct_dot(const uword n_elem, const eT* const A, const eT* const B)
  {
  #if defined(ARMA_USE_ATLAS)
    {
    arma_extra_debug_print("atlas::cx_cblas_dot()");
    
    return atlas::cx_cblas_dot(n_elem, A, B);
    }
  #elif defined(ARMA_USE_BLAS)
    {
    // TODO: work out the mess with zdotu() and zdotu_sub() in BLAS
    return op_dot::direct_dot_arma(n_elem, A, B);
    }
  #else
    {
    return op_dot::direct_dot_arma(n_elem, A, B);
    }
  #endif
  }



//! for two arrays, integral version
template<typename eT>
arma_hot
arma_pure
inline
typename arma_integral_only<eT>::result
op_dot::direct_dot(const uword n_elem, const eT* const A, const eT* const B)
  {
  return op_dot::direct_dot_arma(n_elem, A, B);
  }




//! for three arrays
template<typename eT>
arma_hot
arma_pure
inline
eT
op_dot::direct_dot(const uword n_elem, const eT* const A, const eT* const B, const eT* C)
  {
  arma_extra_debug_sigprint();
  
  eT val = eT(0);
  
  for(uword i=0; i<n_elem; ++i)
    {
    val += A[i] * B[i] * C[i];
    }

  return val;
  }



template<typename T1, typename T2>
arma_hot
arma_inline
typename T1::elem_type
op_dot::apply(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  if( (is_Mat<T1>::value == true) && (is_Mat<T2>::value == true) )
    {
    return op_dot::apply_unwrap(X,Y);
    }
  else
    {
    return op_dot::apply_proxy(X,Y);
    }
  }



template<typename T1, typename T2>
arma_hot
arma_inline
typename T1::elem_type
op_dot::apply_unwrap(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (A.n_elem != B.n_elem), "dot(): objects must have the same number of elements" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  }



template<typename T1, typename T2>
arma_hot
inline
typename T1::elem_type
op_dot::apply_proxy(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
  
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor) && (Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    arma_debug_check( (A.get_n_elem() != B.get_n_elem()), "dot(): objects must have the same number of elements" );
  
    const uword    N  = A.get_n_elem();
          ea_type1 PA = A.get_ea();
          ea_type2 PB = B.get_ea();
    
    eT val1 = eT(0);
    eT val2 = eT(0);
    
    uword i,j;
    
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      val1 += PA[i] * PB[i];
      val2 += PA[j] * PB[j];
      }
    
    if(i < N)
      {
      val1 += PA[i] * PB[i];
      }
    
    return val1 + val2;
    }
  else
    {
    return op_dot::apply_unwrap(A.Q, B.Q);
    }
  }



//
// op_norm_dot



template<typename T1, typename T2>
arma_hot
inline
typename T1::elem_type
op_norm_dot::apply(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor) && (Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const Proxy<T1> A(X.get_ref());
    const Proxy<T2> B(Y.get_ref());
    
    arma_debug_check( (A.get_n_elem() != B.get_n_elem()), "norm_dot(): objects must have the same number of elements" );
    
    const uword    N  = A.get_n_elem();
          ea_type1 PA = A.get_ea();
          ea_type2 PB = B.get_ea();
    
    eT acc1 = eT(0);
    eT acc2 = eT(0);
    eT acc3 = eT(0);
    
    for(uword i=0; i<N; ++i)
      {
      const eT tmpA = PA[i];
      const eT tmpB = PB[i];
      
      acc1 += tmpA * tmpA;
      acc2 += tmpB * tmpB;
      acc3 += tmpA * tmpB;
      }
      
    return acc3 / ( std::sqrt(acc1 * acc2) );
    }
  else
    {
    return op_norm_dot::apply_unwrap(X, Y);
    }
  }



template<typename T1, typename T2>
arma_hot
inline
typename T1::elem_type
op_norm_dot::apply_unwrap(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  
  arma_debug_check( (A.n_elem != B.n_elem), "norm_dot(): objects must have the same number of elements" );
  
  const uword N = A.n_elem;
  
  const eT* A_mem = A.memptr();
  const eT* B_mem = B.memptr();
  
  eT acc1 = eT(0);
  eT acc2 = eT(0);
  eT acc3 = eT(0);
  
  for(uword i=0; i<N; ++i)
    {
    const eT tmpA = A_mem[i];
    const eT tmpB = B_mem[i];
    
    acc1 += tmpA * tmpA;
    acc2 += tmpB * tmpB;
    acc3 += tmpA * tmpB;
    }
    
  return acc3 / ( std::sqrt(acc1 * acc2) );
  }



//
// op_cdot



template<typename T1, typename T2>
arma_hot
arma_inline
typename T1::elem_type
op_cdot::apply(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  const Proxy<T1> A(X.get_ref());
  const Proxy<T2> B(Y.get_ref());
  
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor) || (Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    arma_debug_check( (A.get_n_elem() != B.get_n_elem()), "cdot(): objects must have the same number of elements" );
    
    const uword    N  = A.get_n_elem();
          ea_type1 PA = A.get_ea();
          ea_type2 PB = B.get_ea();
    
    eT val1 = eT(0);
    eT val2 = eT(0);
    
    uword i,j;
    for(i=0, j=1; j<N; i+=2, j+=2)
      {
      val1 += std::conj(PA[i]) * PB[i];
      val2 += std::conj(PA[j]) * PB[j];
      }
    
    if(i < N)
      {
      val1 += std::conj(PA[i]) * PB[i];
      }
    
    return val1 + val2;
    }
  else
    {
    const unwrap< typename Proxy<T1>::stored_type > tmp_A(A.Q);
    const unwrap< typename Proxy<T2>::stored_type > tmp_B(B.Q);
    
    return op_cdot::apply( tmp_A.M, tmp_B.M );
    }
  }



//! @}
