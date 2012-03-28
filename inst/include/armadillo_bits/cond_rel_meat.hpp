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



template<>
template<typename eT>
arma_inline
bool
cond_rel<true>::lt(const eT A, const eT B)
  {
  return (A < B);
  }
  


template<>
template<typename eT>
arma_inline
bool
cond_rel<false>::lt(const eT, const eT)
  {
  return false;
  }
  


template<>
template<typename eT>
arma_inline
bool
cond_rel<true>::gt(const eT A, const eT B)
  {
  return (A > B);
  }
  


template<>
template<typename eT>
arma_inline
bool
cond_rel<false>::gt(const eT, const eT)
  {
  return false;
  }
  


template<>
template<typename eT>
arma_inline
bool
cond_rel<true>::leq(const eT A, const eT B)
  {
  return (A <= B);
  }
  


template<>
template<typename eT>
arma_inline
bool
cond_rel<false>::leq(const eT, const eT)
  {
  return false;
  }
  


template<>
template<typename eT>
arma_inline
bool
cond_rel<true>::geq(const eT A, const eT B)
  {
  return (A >= B);
  }
  


template<>
template<typename eT>
arma_inline
bool
cond_rel<false>::geq(const eT, const eT)
  {
  return false;
  }
  


//! @}
