// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup Op
//! @{



//! Class for storing data required for delayed unary operations,
//! such as the operand (e.g. the matrix to which the operation is to be applied) and the unary operator (e.g. inverse).
//! The operand is stored as a reference (which can be optimised away),
//! while the operator is "stored" through the template definition (op_type).
//! The operands can be 'Mat', 'Row', 'Col', 'Op', and 'Glue'.
//! Note that as 'Glue' can be one of the operands, more than one matrix can be stored.
//!
//! For example, we could have:
//! Op< Glue< Mat, Mat, glue_times >, op_htrans >

template<typename T1, typename op_type>
class Op : public Base<typename T1::elem_type, Op<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool is_row = ( T1::is_col && (is_same_type<op_type, op_strans>::value || is_same_type<op_type, op_htrans>::value || is_same_type<op_type, op_htrans2>::value) );
  static const bool is_col = ( T1::is_row && (is_same_type<op_type, op_strans>::value || is_same_type<op_type, op_htrans>::value || is_same_type<op_type, op_htrans2>::value) );
  
  inline explicit Op(const T1& in_m);
  inline          Op(const T1& in_m, const elem_type in_aux);
  inline          Op(const T1& in_m, const elem_type in_aux,         const uword in_aux_uword_a, const uword in_aux_uword_b);
  inline          Op(const T1& in_m, const uword     in_aux_uword_a, const uword in_aux_uword_b);
  inline          Op(const T1& in_m, const uword     in_aux_uword_a, const uword in_aux_uword_b, const uword in_aux_uword_c, const char junk);
  inline         ~Op();
    
  
  arma_aligned const T1&       m;            //!< storage of reference to the operand (eg. a matrix)
  arma_aligned       elem_type aux;          //!< storage of auxiliary data, user defined format
  arma_aligned       uword     aux_uword_a;  //!< storage of auxiliary data, uword format
  arma_aligned       uword     aux_uword_b;  //!< storage of auxiliary data, uword format
  arma_aligned       uword     aux_uword_c;  //!< storage of auxiliary data, uword format
  
  };



//! @}
