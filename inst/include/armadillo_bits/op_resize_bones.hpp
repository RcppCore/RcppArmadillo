// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_resize
//! @{



class op_resize
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_resize>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_resize>& in);
  };



//! @}
