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


//! \addtogroup op_relational
//! @{



class op_rel_lt_pre
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_lt_pre>& X);
  };



class op_rel_lt_post
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_lt_post>& X);
  };



class op_rel_gt_pre
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_gt_pre>& X);
  };



class op_rel_gt_post
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_gt_post>& X);
  };



class op_rel_lteq_pre
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_lteq_pre>& X);
  };



class op_rel_lteq_post
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_lteq_post>& X);
  };



class op_rel_gteq_pre
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_gteq_pre>& X);
  };



class op_rel_gteq_post
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_gteq_post>& X);
  };



class op_rel_eq
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_eq>& X);
  };



class op_rel_noteq
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<u32>& out, const mtOp<u32, T1, op_rel_noteq>& X);
  };



//! @}
