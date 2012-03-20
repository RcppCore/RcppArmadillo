// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup mtOp
//! @{



template<typename out_eT, typename T1, typename op_type>
class mtOp : public Base<out_eT, mtOp<out_eT, T1, op_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;

  typedef typename T1::elem_type                in_eT;

  static const bool is_row = T1::is_row && is_op_mixed_elem<op_type>::value;
  static const bool is_col = T1::is_col && is_op_mixed_elem<op_type>::value; 
  
  inline explicit mtOp(const T1& in_m);
  inline          mtOp(const T1& in_m, const in_eT in_aux);
  inline          mtOp(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b);
  inline          mtOp(const T1& in_m, const in_eT in_aux,         const uword in_aux_uword_a, const uword in_aux_uword_b);
  
  inline          mtOp(const char junk, const T1& in_m, const out_eT in_aux);
  
  inline         ~mtOp();
    
  
  arma_aligned const T1&    m;            //!< storage of reference to the operand (eg. a matrix)
  arma_aligned       in_eT  aux;          //!< storage of auxiliary data, using the element type as used by T1
  arma_aligned       out_eT aux_out_eT;   //!< storage of auxiliary data, using the element type as specified by the out_eT template parameter
  arma_aligned       uword  aux_uword_a;  //!< storage of auxiliary data, uword format
  arma_aligned       uword  aux_uword_b;  //!< storage of auxiliary data, uword format
  
  };



//! @}
