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



//! \addtogroup fn_trunc_exp
//! @{



template<typename eT>
inline
static
typename arma_float_only<eT>::result
trunc_exp(const eT x)
  {
  if(std::numeric_limits<eT>::is_iec559 && (x >= Math<eT>::log_max() ))
    {
    return std::numeric_limits<eT>::max();
    }
  else
    {
    return std::exp(x);
    }
  }



template<typename eT>
inline
static
typename arma_integral_only<eT>::result
trunc_exp(const eT x)
  {
  return eT( trunc_exp( double(x) ) );
  }



template<typename T>
arma_inline
static
std::complex<T>
trunc_exp(const std::complex<T>& x)
  {
  return std::exp(x);
  }
  
  
  
template<typename T1>
arma_inline
const eOp<T1, eop_trunc_exp>
trunc_exp(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_trunc_exp>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_trunc_exp>
trunc_exp(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_trunc_exp>(A.get_ref());
  }



//! @}
