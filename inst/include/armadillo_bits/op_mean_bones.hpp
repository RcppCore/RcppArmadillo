// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_mean
//! @{


//! Class for finding mean values of a matrix
class op_mean
  {
  public:
  
  template<typename eT>
  inline static eT direct_mean(const eT* const X, const uword N);
  
  template<typename eT>
  inline static eT direct_mean(const Mat<eT>& X, const uword row);
  
  template<typename eT>
  inline static eT direct_mean(const subview<eT>& X);
  
  template<typename eT>
  inline static eT direct_mean(const diagview<eT>& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_mean>& in);
  
  
  template<typename eT>
  inline static eT direct_mean_robust(const eT* const X, const uword N);
  
  template<typename eT>
  inline static eT direct_mean_robust(const Mat<eT>& X, const uword row);
  
  template<typename eT>
  inline static eT direct_mean_robust(const subview<eT>& X);
  
  template<typename eT>
  inline static eT direct_mean_robust(const diagview<eT>& X);
  };



//! @}
