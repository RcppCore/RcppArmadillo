// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup op_strans
//! @{


//! 'matrix transpose' operation

class op_strans
  {
  public:
  
  template<const bool do_flip, const uword row, const uword col>
  struct pos
    {
    static const uword n2 = (do_flip == false) ? (row + col*2) : (col + row*2);
    static const uword n3 = (do_flip == false) ? (row + col*3) : (col + row*3);
    static const uword n4 = (do_flip == false) ? (row + col*4) : (col + row*4);
    };
  
  template<typename eT>
  arma_hot inline static void apply_noalias_tinysq(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename eT>
  arma_hot inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename eT>
  arma_hot inline static void apply(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename T1>
  arma_hot inline static void apply_proxy(Mat<typename T1::elem_type>& out, const T1& X);
  
  template<typename T1>
  arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_strans>& in);
  };



class op_strans2
  {
  public:
  
  template<const bool do_flip, const uword row, const uword col>
  struct pos
    {
    static const uword n2 = (do_flip == false) ? (row + col*2) : (col + row*2);
    static const uword n3 = (do_flip == false) ? (row + col*3) : (col + row*3);
    static const uword n4 = (do_flip == false) ? (row + col*4) : (col + row*4);
    };
  
  template<typename eT>
  arma_hot inline static void apply_noalias_tinysq(Mat<eT>& out, const Mat<eT>& A, const eT val);
  
  template<typename eT>
  arma_hot inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& A, const eT val);
  
  template<typename eT>
  arma_hot inline static void apply(Mat<eT>& out, const Mat<eT>& A, const eT val);
  
  template<typename T1>
  arma_hot inline static void apply_proxy(Mat<typename T1::elem_type>& out, const T1& X, const typename T1::elem_type val);
  };



//! @}
