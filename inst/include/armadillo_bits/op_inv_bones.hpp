// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


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



//! 'invert matrix' operation (positive definite symmetric matrices)
class op_inv_sym
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv_sym>& in);
  };



//! @}
