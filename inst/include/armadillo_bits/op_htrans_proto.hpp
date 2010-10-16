// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_htrans
//! @{


//! 'hermitian transpose' operation (only valid for complex number matrices)

class op_htrans
  {
  public:
  
  template<typename T>
  inline static void apply_noalias(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A);


  template<typename T>
  inline static void apply(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A);
  

  template<typename T, typename T1>
  inline static void apply(Mat< std::complex<T> >& out, const Op<T1,op_htrans>& in);
  
  };


//! @}
