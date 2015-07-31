// Copyright (C) 2008-2015 Conrad Sanderson
// Copyright (C) 2008-2015 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_sum
//! @{


class op_sum
  {
  public:
  
  template<typename T1>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1, op_sum>& in);
  
  template<typename T1>
  arma_hot inline static void apply_noalias(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword dim);
  
  template<typename T1>
  arma_hot inline static void apply_noalias_unwrap(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword dim);
  
  template<typename T1>
  arma_hot inline static void apply_noalias_proxy(Mat<typename T1::elem_type>& out, const Proxy<T1>& P, const uword dim);
  };


//! @}
