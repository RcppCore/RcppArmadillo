// Copyright (C) 2013 Conrad Sanderson
// Copyright (C) 2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_vectorise
//! @{



class op_vectorise_col
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const Op<T1,op_vectorise_col>& in);
  
  template<typename T1> inline static void apply_expr( Mat<typename T1::elem_type>& out, const T1& expr);
  };



class op_vectorise_all
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const Op<T1,op_vectorise_all>& in);
  };



//! @}
