// Copyright (C) 2012 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup unwrap_spmat
//! @{



template<typename T1>
struct unwrap_spmat
  {
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_spmat(const T1& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<eT> M;
  };



template<typename eT>
struct unwrap_spmat< SpMat<eT> >
  {
  inline
  unwrap_spmat(const SpMat<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<eT>& M;
  };



template<typename eT>
struct unwrap_spmat< SpRow<eT> >
  {
  inline
  unwrap_spmat(const SpRow<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpRow<eT>& M;
  };



template<typename eT>
struct unwrap_spmat< SpCol<eT> >
  {
  inline
  unwrap_spmat(const SpCol<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpCol<eT>& M;
  };



template<typename T1, typename spop_type>
struct unwrap_spmat< SpOp<T1, spop_type> >
  {
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_spmat(const SpOp<T1, spop_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<eT> M;
  };



template<typename T1, typename T2, typename spglue_type>
struct unwrap_spmat< SpGlue<T1, T2, spglue_type> >
  {
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_spmat(const SpGlue<T1, T2, spglue_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<eT> M;
  };



template<typename out_eT, typename T1, typename spop_type>
struct unwrap_spmat< mtSpOp<out_eT, T1, spop_type> >
  {
  inline
  unwrap_spmat(const mtSpOp<out_eT, T1, spop_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const SpMat<out_eT> M;
  };



//! @}
