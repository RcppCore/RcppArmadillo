// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_trans
//! @{


//! 'matrix transpose' operation

class op_trans
  {
  public:
  
  template<typename eT>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trans>& in);
  
  
  // inline static void apply_inplace(mat &out);
  
  };



class op_trans2
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trans2>& in);
  };

  
  
//! @}
