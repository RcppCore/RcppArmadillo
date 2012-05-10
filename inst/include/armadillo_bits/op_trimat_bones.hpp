// Copyright (C) 2010-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2012 Conrad Sanderson
// Copyright (C) 2011 Ryan Curtin
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_trimat
//! @{



class op_trimat
  {
  public:
  
  template<typename eT>
  inline static void fill_zeros(Mat<eT>& A, const bool upper);
  
  //
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trimat>& in);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Op<T1,op_htrans>, op_trimat>& in);
  
  //
  
  template<typename eT>
  inline static void apply_htrans(Mat<eT>& out, const Mat<eT>& A, const bool upper, const typename arma_not_cx<eT>::result* junk = 0);
  
  template<typename eT>
  inline static void apply_htrans(Mat<eT>& out, const Mat<eT>& A, const bool upper, const typename arma_cx_only<eT>::result* junk = 0);
  };



//! @}
