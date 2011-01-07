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


//! \addtogroup fn_misc
//! @{



//! \brief
//! Generate a vector with 'num' elements.
//! The values of the elements linearly increase from 'start' upto (and including) 'end'.

template<typename vec_type>
inline
vec_type
linspace
  (
  const typename vec_type::pod_type start,
  const typename vec_type::pod_type end,
  const u32 num,
  const typename arma_Mat_Col_Row_only<vec_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  arma_debug_check( (num < 2), "linspace(): num must be >= 2");
  
  typedef typename vec_type::elem_type eT;
  typedef typename vec_type::pod_type   T;
  
  const u32 n_rows = (is_Row<vec_type>::value == true) ? 1   : num;
  const u32 n_cols = (is_Row<vec_type>::value == true) ? num : 1;
  
  Mat<eT> x(n_rows, n_cols);
  eT* x_mem = x.memptr();
  
  const u32 num_m1 = num - 1;
  
  if(is_non_integral<T>::value == true)
    {
    const T delta = (end-start)/T(num_m1);
    
    for(u32 i=0; i<num_m1; ++i)
      {
      x_mem[i] = eT(start + i*delta);
      }
    
    x_mem[num_m1] = eT(end);
    }
  else
    {
    const double delta = (end >= start) ? double(end-start)/double(num_m1) : -double(start-end)/double(num_m1);
    
    for(u32 i=0; i<num_m1; ++i)
      {
      x_mem[i] = eT(double(start) + i*delta);
      }
    
    x_mem[num_m1] = eT(end);
    }
  
  return x;
  }



inline
mat
linspace(const double start, const double end, const u32 num)
  {
  arma_extra_debug_sigprint();
  return linspace<mat>(start, end, num);
  }



//
// log_add

template<typename eT>
inline
typename arma_float_only<eT>::result
log_add(eT log_a, eT log_b)
  {
  if(log_a < log_b)
    {
    std::swap(log_a, log_b);
    }
  
  const eT negdelta = log_b - log_a;
  
  if( (negdelta < Math<eT>::log_min()) || (arma_isfinite(negdelta) == false) )
    {
    return log_a;
    }
  else
    {
    #if defined(ARMA_HAVE_LOG1P)
      return (log_a + log1p(std::exp(negdelta)));
    #else
      return (log_a + std::log(1.0 + std::exp(negdelta)));
    #endif
    }
  }



//! @}
