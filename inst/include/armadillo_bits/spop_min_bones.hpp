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


//! \addtogroup spop_min
//! @{


class spop_min
  {
  public:
  
  template<typename T1>
  inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_min>& in);
  
  template<typename T1>
  inline static void apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& p, const uword dim, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0);
  
  template<typename T1>
  inline static void apply_noalias(SpMat<typename T1::elem_type>& out, const SpProxy<T1>& p, const uword dim, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0);
  
  
  
  template<typename T1>
  inline static typename T1::elem_type vector_min(const T1& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0);
  
  template<typename T1>
  inline static typename T1::elem_type vector_min(const T1& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0);
  };


//! @}
