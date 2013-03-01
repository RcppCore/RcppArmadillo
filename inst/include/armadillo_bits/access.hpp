// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup access
//! @{


class access
  {
  public:
  
  //! internal function to allow modification of data declared as read-only (use with caution)
  template<typename T1> arma_inline static T1&  rw (const T1& x)        { return const_cast<T1& >(x); }
  template<typename T1> arma_inline static T1*& rwp(const T1* const& x) { return const_cast<T1*&>(x); }
  
  //! internal function to obtain the real part of either a plain number or a complex number
  template<typename eT> arma_inline static const eT& tmp_real(const eT&              X) { return X;        }
  template<typename  T> arma_inline static const   T tmp_real(const std::complex<T>& X) { return X.real(); }
  };


//! @}
