// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup spop_misc
//! @{


class spop_scalar_times
  {
  public:
  
  template<typename T1>
  arma_hot inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_scalar_times>& in);
  };



class spop_square
  {
  public:
  
  template<typename T1>
  arma_hot inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_square>& in);
  };



class spop_sqrt
  {
  public:
  
  template<typename T1>
  arma_hot inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_sqrt>& in);
  };



class spop_abs
  {
  public:
  
  template<typename T1>
  arma_hot inline static void apply(SpMat<typename T1::elem_type>& out, const SpOp<T1, spop_abs>& in);
  };



class spop_cx_abs
  {
  public:
  
  template<typename T1>
  arma_hot inline static void apply(SpMat<typename T1::pod_type>& out, const mtSpOp<typename T1::pod_type, T1, spop_cx_abs>& in);
  };



//! @}
