// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_cov
//! @{



class op_cov
  {
  public:
  
  template<typename eT> inline static void direct_cov(Mat<eT>&                out, const Mat<eT>& X,                const uword norm_type);
  template<typename  T> inline static void direct_cov(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& X, const uword norm_type);
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cov>& in);
  };



//! @}
