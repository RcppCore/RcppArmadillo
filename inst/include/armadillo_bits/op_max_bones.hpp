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


//! \addtogroup op_max
//! @{



//! Class for finding maximum values in a matrix
class op_max
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_max>& in);
  
  
  //
  // for non-complex numbers
  
  template<typename eT>
  inline static eT direct_max(const eT* const X, const uword N);
  
  template<typename eT>
  inline static eT direct_max(const eT* const X, const uword N, uword& index_of_max_val);
  
  template<typename eT>
  inline static eT direct_max(const Mat<eT>& X, const uword row);
  
  template<typename eT>
  inline static eT max(const subview<eT>& X);
  
  template<typename T1>
  inline static typename arma_not_cx<typename T1::elem_type>::result max(const Base<typename T1::elem_type, T1>& X);
  
  
  //
  // for complex numbers
  
  template<typename T>
  inline static std::complex<T> direct_max(const std::complex<T>* const X, const uword n_elem);
  
  template<typename T>
  inline static std::complex<T> direct_max(const std::complex<T>* const X, const uword n_elem, uword& index_of_max_val);
  
  template<typename T>
  inline static std::complex<T> direct_max(const Mat< std::complex<T> >& X, const uword row);
  
  template<typename T>
  inline static std::complex<T> max(const subview< std::complex<T> >& X);
  
  template<typename T1>
  inline static typename arma_cx_only<typename T1::elem_type>::result max(const Base<typename T1::elem_type, T1>& X);
  };



//! @}
