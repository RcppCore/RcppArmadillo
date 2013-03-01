// Copyright (C) 2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


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
