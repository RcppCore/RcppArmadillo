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


//! \addtogroup Glue
//! @{



//! Class for storing data required for delayed binary operations,
//! such as the operands (e.g. two matrices) and the binary operator (e.g. addition).
//! The operands are stored as references (which can be optimised away),
//! while the operator is "stored" through the template definition (glue_type).
//! The operands can be 'Mat', 'Row', 'Col', 'Op', and 'Glue'.
//! Note that as 'Glue' can be one of the operands, more than two matrices can be stored.
//!
//! For example, we could have: Glue<Mat, Mat, glue_times>
//! 
//! Another example is: Glue< Op<Mat, op_htrans>, Op<Mat, op_inv>, glue_times >



template<typename T1, typename T2, typename glue_type>
class Glue : public Base<typename T1::elem_type, Glue<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  arma_inline  Glue(const T1& in_A, const T2& in_B);
  arma_inline  Glue(const T1& in_A, const T2& in_B, const uword in_aux_uword);
  arma_inline ~Glue();
  
  const T1&   A;          //!< first operand
  const T2&   B;          //!< second operand
        uword aux_uword;  //!< storage of auxiliary data, uword format
  };



//! @}
