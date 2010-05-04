// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_find
//! @{



template<typename eT>
arma_hot
inline
u32
find_helper(u32* out_mem, const eT* in_mem, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  u32 n_nz = 0;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    if(in_mem[i] != eT(0))
      {
      out_mem[n_nz] = i;
      ++n_nz;
      }
    }
   
  return n_nz;
  }


// TODO:
// once the relational operators are changed into delayed operations,
// find() will need special handling of relational operators to gain speed
template<typename T1>
inline
Mat<u32>
find(const Base<typename T1::elem_type,T1>& X, const u32 k = 0, const char* direction = "first")
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  const u32 n_elem = A.n_elem;
  
  Mat<u32> indices(n_elem, 1);
  
  const u32 n_nz = find_helper(indices.memptr(), A.memptr(), n_elem);
  
  const char sig = direction[0];
  
  // return the indices of all n_nz elements, otherwise return an empty matrix
  if(n_nz > 0)
    {
    if(sig == 'f' || sig == 'F')   // "first"
      {
      return ( (k > 0 && k <= n_nz) ? indices.rows(0,      k-1   ) : indices.rows(0, n_nz-1) );
      }
    else
    if(sig == 'l' || sig == 'L')   // "last"
      {
      return ( (k > 0 && k <= n_nz) ? indices.rows(n_nz-k, n_nz-1) : indices.rows(0, n_nz-1) );
      }
    else
      {
      arma_stop("find(): 3rd input argument must be \"first\" or \"last\"");
      
      return Mat<u32>();
      }
    }
  else
    {
    return Mat<u32>();
    }
  }



//! @}
