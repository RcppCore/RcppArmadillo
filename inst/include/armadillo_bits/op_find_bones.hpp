// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// Copyright (C) 2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_find
//! @{



class op_find
  {
  public:
  
  template<typename T1>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const Base<typename T1::elem_type, T1>& X
    );
  
  template<typename T1, typename op_type>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const mtOp<uword, T1, op_type>& X,
    const typename arma_op_rel_only<op_type>::result junk1 = 0,
    const typename arma_not_cx<typename T1::elem_type>::result junk2 = 0
    );
  
  template<typename T1, typename op_type>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const mtOp<uword, T1, op_type>& X,
    const typename arma_op_rel_only<op_type>::result junk1 = 0,
    const typename arma_cx_only<typename T1::elem_type>::result junk2 = 0
    );
  
  template<typename T1, typename T2, typename glue_type>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const mtGlue<uword, T1, T2, glue_type>& X,
    const typename arma_glue_rel_only<glue_type>::result junk1 = 0,
    const typename arma_not_cx<typename T1::elem_type>::result junk2 = 0,
    const typename arma_not_cx<typename T2::elem_type>::result junk3 = 0
    );
  
  template<typename T1, typename T2, typename glue_type>
  inline static uword
  helper
    (
    Mat<uword>& indices,
    const mtGlue<uword, T1, T2, glue_type>& X,
    const typename arma_glue_rel_only<glue_type>::result junk1 = 0,
    const typename arma_cx_only<typename T1::elem_type>::result junk2 = 0,
    const typename arma_cx_only<typename T2::elem_type>::result junk3 = 0
    );
  
  template<typename T1>
  inline static void apply(Mat<uword>& out, const mtOp<uword, T1, op_find>& X);
  };



//! @}
