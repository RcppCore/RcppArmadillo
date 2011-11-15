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


//! \addtogroup glue_relational
//! @{



class glue_rel_lt
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat <uword>& out, const mtGlue<uword, T1, T2, glue_rel_lt>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube <uword>& out, const mtGlueCube<uword, T1, T2, glue_rel_lt>& X);
  };



class glue_rel_gt
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat <uword>& out, const mtGlue<uword, T1, T2, glue_rel_gt>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube <uword>& out, const mtGlueCube<uword, T1, T2, glue_rel_gt>& X);
  };



class glue_rel_lteq
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat <uword>& out, const mtGlue<uword, T1, T2, glue_rel_lteq>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube <uword>& out, const mtGlueCube<uword, T1, T2, glue_rel_lteq>& X);
  };



class glue_rel_gteq
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat <uword>& out, const mtGlue<uword, T1, T2, glue_rel_gteq>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube <uword>& out, const mtGlueCube<uword, T1, T2, glue_rel_gteq>& X);
  };



class glue_rel_eq
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat <uword>& out, const mtGlue<uword, T1, T2, glue_rel_eq>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube <uword>& out, const mtGlueCube<uword, T1, T2, glue_rel_eq>& X);
  };



class glue_rel_noteq
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat <uword>& out, const mtGlue<uword, T1, T2, glue_rel_noteq>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Cube <uword>& out, const mtGlueCube<uword, T1, T2, glue_rel_noteq>& X);
  };



//! @}
