// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup glue_mixed
//! @{



class glue_mixed_times
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_times>& X);
  };



class glue_mixed_plus
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_plus>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_plus>& X);
  };



class glue_mixed_minus
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_minus>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_minus>& X);
  };



class glue_mixed_div
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_div>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_div>& X);
  };



class glue_mixed_schur
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename eT_promoter<T1,T2>::eT>& out, const mtGlue<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_schur>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube<typename eT_promoter<T1,T2>::eT>& out, const mtGlueCube<typename eT_promoter<T1,T2>::eT, T1, T2, glue_mixed_schur>& X);
  };



//! @}
