// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_max
//! @{



//! Class for finding maximum values in a matrix
class op_max
  {
  public:
  
  template<typename eT>
  inline static eT direct_max(const eT* const X, const u32 N);
  
  template<typename eT>
  inline static eT direct_max(const subview<eT>& X);
  
  template<typename eT>
  inline static eT direct_max(const diagview<eT>& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_max>& in);
  
  
  //
  // for complex numbers

  template<typename T>
  inline static std::complex<T> direct_max(const std::complex<T>* const X, const u32 n_elem);
  
  template<typename T>
  inline static std::complex<T> direct_max(const subview< std::complex<T> >&X);
  
  template<typename T>
  inline static std::complex<T> direct_max(const diagview< std::complex<T> >&X);
  
  template<typename T, typename T1>
  inline static void apply(Mat< std::complex<T> >& out, const Op<T1,op_max>& in);
  
  };



//! @}
