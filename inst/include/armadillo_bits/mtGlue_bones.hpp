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


//! \addtogroup mtGlue
//! @{



template<typename out_eT, typename T1, typename T2, typename glue_type>
class mtGlue : public Base<out_eT, mtGlue<out_eT, T1, T2, glue_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  
  static const bool is_row = ( is_glue_mixed_elem<glue_type>::value && (T1::is_row || T2::is_row) ) || ( is_glue_mixed_times<glue_type>::value && T1::is_row );
  static const bool is_col = ( is_glue_mixed_elem<glue_type>::value && (T1::is_col || T2::is_col) ) || ( is_glue_mixed_times<glue_type>::value && T2::is_col );
  
  arma_inline  mtGlue(const T1& in_A, const T2& in_B);
  arma_inline  mtGlue(const T1& in_A, const T2& in_B, const uword in_aux_uword);
  arma_inline ~mtGlue();
  
  arma_aligned const T1&   A;         //!< first operand
  arma_aligned const T2&   B;         //!< second operand
  arma_aligned       uword aux_uword; //!< storage of auxiliary data, uword format
  };



//! @}
