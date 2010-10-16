// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_rank
//! @{



template<typename T1>
inline
arma_warn_unused
u32
rank
  (
  const Base<typename T1::elem_type,T1>& X,
        typename T1::pod_type            tol = 0.0,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  u32    X_n_rows;
  u32    X_n_cols;
  Col<T> s;
  
  const bool status = auxlib::svd(s, X, X_n_rows, X_n_cols);
  
  if(status == true)
    {
    if(tol == T(0))
      {
      tol = (std::max)(X_n_rows, X_n_cols) * eop_aux::direct_eps(max(s));
      }
    
    // count non zero valued elements in s
    
    const T*  s_mem  = s.memptr();
    const u32 n_elem = s.n_elem;
          u32 count  = 0;
    
    for(u32 i=0; i<n_elem; ++i)
      {
      if(s_mem[i] > tol)
        {
        ++count;
        }
      }
    
    return count;
    }
  else
    {
    arma_print("rank(): failed to converge");
    return u32(0);
    }
  }



//! @}
