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


//! \addtogroup cond_rel
//! @{


//
// for preventing pedantic compiler warnings

template<const bool do_eval>
class cond_rel
  {
  public:
  
  template<typename eT> arma_inline static bool lt(const eT A, const eT B);
  template<typename eT> arma_inline static bool gt(const eT A, const eT B);

  template<typename eT> arma_inline static bool leq(const eT A, const eT B);
  template<typename eT> arma_inline static bool geq(const eT A, const eT B);
  };



//! @}
