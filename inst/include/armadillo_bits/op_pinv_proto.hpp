// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_pinv
//! @{



class op_pinv
  {
  public:
  
  template<typename eT> inline static void direct_pinv(Mat<eT>&                 out, const Mat<eT>& X,                 eT tol);
  template<typename eT> inline static void direct_pinv(Mat< std::complex<eT> >& out, const Mat< std::complex<eT> >& X, eT tol);

  template<typename T1> inline static void apply(Mat<typename T1::pod_type>&                 out, const Op<T1,op_pinv>& in);
  template<typename T1> inline static void apply(Mat< std::complex<typename T1::pod_type> >& out, const Op<T1,op_pinv>& in);
  };



//! @}
