// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_htrans
//! @{


//! 'hermitian transpose' operation

class op_htrans
  {
  public:
  
  template<typename eT>
  arma_hot arma_inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk = 0);
  
  template<typename eT>
  arma_hot inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk = 0);
  
  //
  
  template<typename eT>
  arma_hot arma_inline static void apply(Mat<eT>& out, const Mat<eT>& A, const typename arma_not_cx<eT>::result* junk = 0);
  
  template<typename eT>
  arma_hot inline static void apply(Mat<eT>& out, const Mat<eT>& A, const typename arma_cx_only<eT>::result* junk = 0);
  
  //
  
  template<typename T1>
  arma_hot inline static void apply_proxy(Mat<typename T1::elem_type>& out, const T1& X);
  
  //
  
  template<typename T1>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0);
  
  template<typename T1>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0);
  
  //
  
  template<typename T1>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Op< Op<T1, op_trimat>, op_htrans>& in);
  };



class op_htrans2
  {
  public:
  
  template<typename eT>
  arma_hot inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& A, const eT val);
  
  template<typename eT>
  arma_hot inline static void apply(Mat<eT>& out, const Mat<eT>& A, const eT val);
  
  //
  
  template<typename T1>
  arma_hot inline static void apply_proxy(Mat<typename T1::elem_type>& out, const T1& X, const typename T1::elem_type val);
  
  //
  
  template<typename T1>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0);
  
  template<typename T1>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_htrans2>& in, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0);
  };



//! @}
