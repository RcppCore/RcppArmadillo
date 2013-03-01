// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_inv
//! @{



//! 'invert matrix' operation (general matrices)
class op_inv
  {
  public:
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Mat<eT>& A, const bool slow = false);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv>& in);
  
  template<typename T1>
  inline static void apply_diag(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type, T1>& X);
  };



//! 'invert matrix' operation (triangular matrices)
class op_inv_tr
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_tr>& in);
  };



//! 'invert matrix' operation (symmetric positive definite matrices)
class op_inv_sympd
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_sympd>& in);
  };



//! @}
