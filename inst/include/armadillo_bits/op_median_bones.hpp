// Copyright (C) 2009-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_median
//! @{


template<typename T>
struct arma_cx_median_packet
  {
  T   val;
  uword index;
  };



template<typename T>
arma_inline
bool
operator< (const arma_cx_median_packet<T>& A, const arma_cx_median_packet<T>& B)
  {
  return A.val < B.val;
  }



//! Class for finding median values of a matrix
class op_median
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_median>& in);
  
  template<typename T, typename T1>
  inline static void apply(Mat< std::complex<T> >& out, const Op<T1,op_median>& in);
  
  //
  //
  
  template<typename T1>
  inline static typename T1::elem_type median_vec(const T1& X, const typename arma_not_cx<typename T1::elem_type>::result* junk = 0);
  
  template<typename T1>
  inline static typename T1::elem_type median_vec(const T1& X, const typename arma_cx_only<typename T1::elem_type>::result* junk = 0);
  
  //
  //
  
  template<typename eT>
  inline static eT direct_median(std::vector<eT>& X);
  
  template<typename T>
  inline static void direct_cx_median_index(uword& out_index1, uword& out_index2, std::vector< arma_cx_median_packet<T> >& X);
  };



//! @}
