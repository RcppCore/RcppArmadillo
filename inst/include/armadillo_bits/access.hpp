// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup access
//! @{


class access
  {
  public:
  
  //! internal function to allow modification of data declared as read-only
  template<typename T1> arma_inline static T1& rw(const T1& x) { return const_cast<T1&>(x); }

  //! internal function to obtain the real part of either a plain number or a complex number
  template<typename eT> arma_inline static const eT& tmp_real(const eT&              X) { return X;        }
  template<typename  T> arma_inline static const   T tmp_real(const std::complex<T>& X) { return X.real(); }
  };


//! @}
