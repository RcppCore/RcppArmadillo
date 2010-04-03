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


//! \addtogroup Col
//! @{

//! Class for column vectors (matrices with only column)

template<typename eT>
class Col : public Mat<eT>, public BaseVec< eT, Col<eT> >
  {
  public:
  
  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  
  
  inline                     Col();
  inline explicit            Col(const u32 n_elem);
  inline                     Col(const u32 in_rows, const u32 in_cols);
  
  inline                     Col(const char*        text);
  inline const Col&    operator=(const char*        text);
  inline                     Col(const std::string& text);
  inline const Col&    operator=(const std::string& text);
  
  inline                     Col(const Col& X);
  inline const Col&    operator=(const Col& X);
  
  //inline explicit            Col(const Mat<eT>& X);
  inline                     Col(const Mat<eT>& X);
  inline const Col&    operator=(const Mat<eT>& X);
  inline const Col&   operator*=(const Mat<eT>& X);
  
  inline Col(      eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const bool copy_aux_mem = true);
  inline Col(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols);
  
  inline Col(      eT* aux_mem, const u32 aux_length, const bool copy_aux_mem = true);
  inline Col(const eT* aux_mem, const u32 aux_length);
  
  template<typename T1, typename T2>
  inline explicit Col(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);
  
  inline                     Col(const subview<eT>& X);
  inline const Col&    operator=(const subview<eT>& X);
  inline const Col&   operator*=(const subview<eT>& X);
  
  inline                     Col(const subview_cube<eT>& X);
  inline const Col&    operator=(const subview_cube<eT>& X);
  inline const Col&   operator*=(const subview_cube<eT>& X);
  
  inline                     Col(const diagview<eT>& X);
  inline const Col&    operator=(const diagview<eT>& X);
  inline const Col&   operator*=(const diagview<eT>& X);
  
  arma_inline eT& row(const u32 row_num);
  arma_inline eT  row(const u32 row_num) const;
  
  arma_inline       subview_col<eT> rows(const u32 in_row1, const u32 in_row2);
  arma_inline const subview_col<eT> rows(const u32 in_row1, const u32 in_row2) const;
  
  template<typename T1, typename op_type> inline                   Col(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Col&  operator=(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Col& operator*=(const Op<T1, op_type>& X);
  
  template<typename T1, typename eop_type> inline                   Col(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Col&  operator=(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Col& operator*=(const eOp<T1, eop_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline                   Col(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Col&  operator=(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Col& operator*=(const Glue<T1, T2, glue_type>& X);
  
  template<typename T1, typename T2, typename eglue_type> inline                   Col(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Col&  operator=(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Col& operator*=(const eGlue<T1, T2, eglue_type>& X);
  
  inline void  set_size(const u32 n_elem);
  inline void  set_size(const u32 n_rows, const u32 n_cols);
  
  template<typename eT2>
  inline void copy_size(const Mat<eT2>& m);
  
  inline void zeros();
  inline void zeros(const u32 n_elem);
  inline void zeros(const u32 n_rows, const u32 n_cols);
  
  inline void ones();
  inline void ones(const u32 n_elem);
  inline void ones(const u32 n_rows, const u32 n_cols);
  
  inline void load(const std::string   name, const file_type type = auto_detect);
  inline void load(      std::istream& is,   const file_type type = auto_detect);
  };


//! @}
