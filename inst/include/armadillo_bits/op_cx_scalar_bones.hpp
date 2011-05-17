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


//! \addtogroup op_cx_scalar
//! @{



class op_cx_scalar_times
  {
  public:
  
  template<typename T1>
  inline static void
  apply
    (
          Mat< typename std::complex<typename T1::pod_type> >& out,
    const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>& X
    );
  
  template<typename T1>
  inline static void
  apply
    (
             Cube< typename std::complex<typename T1::pod_type> >& out,
    const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_times>& X
    );

  };



class op_cx_scalar_plus
  {
  public:
  
  template<typename T1>
  inline static void
  apply
    (
          Mat< typename std::complex<typename T1::pod_type> >& out,
    const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_plus>& X
    );
  
  template<typename T1>
  inline static void
  apply
    (
             Cube< typename std::complex<typename T1::pod_type> >& out,
    const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_plus>& X
    );

  };



class op_cx_scalar_minus_pre
  {
  public:
  
  template<typename T1>
  inline static void
  apply
    (
          Mat< typename std::complex<typename T1::pod_type> >& out,
    const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_pre>& X
    );
  
  template<typename T1>
  inline static void
  apply
    (
             Cube< typename std::complex<typename T1::pod_type> >& out,
    const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_pre>& X
    );

  };



class op_cx_scalar_minus_post
  {
  public:
  
  template<typename T1>
  inline static void
  apply
    (
          Mat< typename std::complex<typename T1::pod_type> >& out,
    const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_post>& X
    );
  
  template<typename T1>
  inline static void
  apply
    (
             Cube< typename std::complex<typename T1::pod_type> >& out,
    const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_minus_post>& X
    );

  };



class op_cx_scalar_div_pre
  {
  public:
  
  template<typename T1>
  inline static void
  apply
    (
          Mat< typename std::complex<typename T1::pod_type> >& out,
    const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_pre>& X
    );
  
  template<typename T1>
  inline static void
  apply
    (
             Cube< typename std::complex<typename T1::pod_type> >& out,
    const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_pre>& X
    );

  };



class op_cx_scalar_div_post
  {
  public:
  
  template<typename T1>
  inline static void
  apply
    (
          Mat< typename std::complex<typename T1::pod_type> >& out,
    const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_post>& X
    );
  
  template<typename T1>
  inline static void
  apply
    (
             Cube< typename std::complex<typename T1::pod_type> >& out,
    const mtOpCube<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_post>& X
    );

  };



//! @}
