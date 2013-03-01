// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_repmat
//! @{



class op_repmat
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_repmat>& in);
  };



//! @}
