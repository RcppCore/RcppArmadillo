// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
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
rank(const Base<typename T1::elem_type,T1>& X, typename T1::pod_type tol = 0.0)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  Col<T> s;
  const bool status = auxlib::svd(s, A);
  
  if(status == true)
    {
    if(tol == T(0))
      {
      tol = (std::max)(A.n_rows, A.n_cols) * eop_aux::direct_eps(max(s));
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
    arma_print("rank(): singular value decomposition failed");
    return u32(0);
    }
   
  }



//! @}
