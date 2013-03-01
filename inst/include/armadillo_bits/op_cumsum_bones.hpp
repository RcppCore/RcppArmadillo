// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_cumsum
//! @{



class op_cumsum_mat
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cumsum_mat>& in);
  };



class op_cumsum_vec
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cumsum_vec>& in);
  };

//! @}
