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



//! \addtogroup glue_kron
//! @{



class glue_kron
  {
  public:

  template<typename eT> inline static void direct_kron(Mat<eT>&                out, const Mat<eT>&                A, const Mat<eT>&                B);
  template<typename T>  inline static void direct_kron(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A, const Mat<T>&                 B);
  template<typename T>  inline static void direct_kron(Mat< std::complex<T> >& out, const Mat<T>&                 A, const Mat< std::complex<T> >& B);
  
  template<typename T1, typename T2>   inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_kron>& X);
  };



//! @}

