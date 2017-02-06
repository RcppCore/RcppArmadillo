// Copyright (C) 2016 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au



//! \addtogroup glue_polyfit
//! @{



class glue_polyfit
  {
  public:
  
  template<typename eT> inline static bool apply_noalias(Mat<eT>& out, const Col<eT>& X, const Col<eT>& Y, const uword N);
  
  template<typename T1, typename T2> inline static bool apply_direct(Mat<typename T1::elem_type>& out, const Base<typename T1::elem_type,T1>& X_expr, const Base<typename T1::elem_type, T2>& Y_expr, const uword N);
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_polyfit>& expr);
  };



//! @}

