// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup eglue_cube_core
//! @{



template<typename eglue_type>
struct eglue_cube_core
  {
  
  template<typename T1, typename T2> arma_inline static typename T1::elem_type get_elem(const eGlueCube<T1, T2, eglue_type>& x, const u32 i);
  template<typename T1, typename T2> arma_inline static typename T1::elem_type get_elem(const eGlueCube<T1, T2, eglue_type>& x, const u32 row, const u32 col, const u32 slice);
  
  template<typename T1, typename T2> arma_hot inline static void apply(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  
  template<typename T1, typename T2> arma_hot inline static void apply_proxy (Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_unwrap(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_plus (Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_minus(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_schur(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  template<typename T1, typename T2> arma_hot inline static void apply_inplace_div  (Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x);
  
  };



class eglue_cube_plus;
class eglue_cube_minus;
class eglue_cube_div;
class eglue_cube_schur;



//! @}
