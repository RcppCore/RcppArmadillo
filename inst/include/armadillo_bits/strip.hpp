// Copyright (C) 2010-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup strip
//! @{



template<typename T1>
struct strip_diagmat
  {
  typedef T1 stored_type;
  
  inline strip_diagmat(const Base<typename T1::elem_type, T1>& X)
    : M(X.get_ref())
    {
    arma_extra_debug_sigprint();
    }
  
  static const bool do_diagmat = false;
  
  const T1& M;
  };



template<typename T1>
struct strip_diagmat< Op<T1, op_diagmat> >
  {
  typedef T1 stored_type;
  
  inline strip_diagmat(const Op<T1, op_diagmat>& X)
    : M(X.m)
    {
    arma_extra_debug_sigprint();
    }
  
  static const bool do_diagmat = true;
  
  const T1& M;
  };



template<typename T1>
struct strip_inv
  {
  typedef T1 stored_type;
  
  inline strip_inv(const T1& X)
    : M(X)
    {
    arma_extra_debug_sigprint();
    }
  
  static const bool do_inv = false;
  
  const T1& M;
  };



template<typename T1>
struct strip_inv< Op<T1, op_inv> >
  {
  typedef T1 stored_type;
  
  inline strip_inv(const Op<T1, op_inv>& X)
    : M(X.m)
    {
    arma_extra_debug_sigprint();
    }
  
  static const bool do_inv = true;
  
  const T1& M;
  };



//! @}
