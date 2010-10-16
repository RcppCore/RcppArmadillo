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


//! \addtogroup fn_trunc_log
//! @{



template<typename eT>
inline 
static
typename arma_float_only<eT>::result
trunc_log(const eT x)
  {
  if(std::numeric_limits<eT>::is_iec559)
    {
    if(x == std::numeric_limits<eT>::infinity())
      {
      return Math<eT>::log_max();
      }
    else
      {
      return (x <= eT(0)) ? Math<eT>::log_min() : std::log(x);
      }
    }
  else
    {
    return std::log(x);
    }
  }



template<typename eT>
inline 
static
typename arma_integral_only<eT>::result
trunc_log(const eT x)
  {
  return eT( trunc_log( double(x) ) );
  }



template<typename T>
inline 
static
std::complex<T>
trunc_log(const std::complex<T>& x)
  {
  return std::log(x);
  }



template<typename T1>
arma_inline
const eOp<T1, eop_trunc_log>
trunc_log(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_trunc_log>(A.get_ref());
  }



template<typename T1>
arma_inline
const eOpCube<T1, eop_trunc_log>
trunc_log(const BaseCube<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return eOpCube<T1, eop_trunc_log>(A.get_ref());
  }



//! @}
