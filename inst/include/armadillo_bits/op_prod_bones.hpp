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


//! \addtogroup op_prod
//! @{

//! Class for finding products of values in a matrix  (e.g. along rows or columns)
class op_prod
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1, op_prod>& in);
  
  template<typename eT>
  inline static eT prod(const subview<eT>& S);
  
  template<typename T1>
  inline static typename T1::elem_type prod(const Base<typename T1::elem_type,T1>& X);
  };


//! @}
