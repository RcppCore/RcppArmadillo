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
  inline static void apply_noalias_tinysq(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename eT>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Mat<eT>& A);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_strans>& in);
  
  
  // inline static void apply_inplace(mat &out);
  
  };



//! @}
