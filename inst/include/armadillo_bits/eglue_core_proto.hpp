// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup eglue_core
//! @{



template<typename eglue_type>
struct eglue_core
  {
  
  // matrices
  
  template<typename T1, typename T2> arma_hot inline static void apply(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x);
  
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_plus (Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_minus(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_schur(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_div  (Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x);
  
  
  // cubes
  
  template<typename T1, typename T2> arma_hot inline static void apply(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_plus (Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_minus(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_schur(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_div  (Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  };



class eglue_plus;
class eglue_minus;
class eglue_div;
class eglue_schur;



//! @}
