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


//! \addtogroup Row
//! @{

//! Class for row vectors (matrices with only one row)

template<typename eT>
class Row : public Mat<eT>, public BaseVec< eT, Row<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  
  inline                     Row();
  inline explicit            Row(const u32 N);
  inline                     Row(const u32 in_rows, const u32 in_cols);
  
  inline                     Row(const char*        text);
  inline const Row&    operator=(const char*        text);
  inline                     Row(const std::string& text);
  inline const Row&    operator=(const std::string& text);
  
  inline                     Row(const Row& X);
  inline const Row&    operator=(const Row& X);
  inline const Row&   operator*=(const Row& X);
  
  //inline explicit            Row(const Mat<eT>& X);
  inline                     Row(const Mat<eT>& X);
  inline const Row&    operator=(const Mat<eT>& X);
  inline const Row&   operator*=(const Mat<eT>& X);

  inline Row(      eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const bool copy_aux_mem = true);
  inline Row(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols);
  
  inline Row(      eT* aux_mem, const u32 aux_length, const bool copy_aux_mem = true);
  inline Row(const eT* aux_mem, const u32 aux_length);

  template<typename T1, typename T2>
  inline explicit Row(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);

  inline                     Row(const subview<eT>& X);
  inline const Row&    operator=(const subview<eT>& X);
  inline const Row&   operator*=(const subview<eT>& X);

  inline                     Row(const subview_cube<eT>& X);
  inline const Row&    operator=(const subview_cube<eT>& X);
  inline const Row&   operator*=(const subview_cube<eT>& X);
  
  inline explicit            Row(const diagview<eT>& X);
  inline const Row&    operator=(const diagview<eT>& X);
  inline const Row&   operator*=(const diagview<eT>& X);

  inline mat_injector<Row> operator<<(const eT val);
  
  arma_inline eT& col(const u32 col_num);
  arma_inline eT  col(const u32 col_num) const;
  
  arma_inline       subview_row<eT> cols(const u32 in_col1, const u32 in_col2);
  arma_inline const subview_row<eT> cols(const u32 in_col1, const u32 in_col2) const;
  
  template<typename T1, typename op_type> inline                   Row(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Row&  operator=(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Row& operator*=(const Op<T1, op_type>& X);
  
  template<typename T1, typename eop_type> inline                   Row(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Row&  operator=(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Row& operator*=(const eOp<T1, eop_type>& X);
  
  template<typename T1, typename op_type> inline                   Row(const mtOp<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline const Row&  operator=(const mtOp<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline const Row& operator*=(const mtOp<eT, T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline                   Row(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Row&  operator=(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Row& operator*=(const Glue<T1, T2, glue_type>& X);
  
  template<typename T1, typename T2, typename eglue_type> inline                   Row(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Row&  operator=(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Row& operator*=(const eGlue<T1, T2, eglue_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline                   Row(const mtGlue<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Row&  operator=(const mtGlue<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Row& operator*=(const mtGlue<eT, T1, T2, glue_type>& X);
  
#ifdef ARMA_EXTRA_ROW_PROTO
#include ARMA_EXTRA_ROW_PROTO
#endif

  inline void  set_size(const u32 N);
  inline void  set_size(const u32 n_rows, const u32 n_cols);
  inline void   reshape(const u32 n_rows, const u32 n_cols, const u32 dim = 0);
  
  template<typename eT2>
  inline void copy_size(const Mat<eT2>& m);
  
  inline void zeros();
  inline void zeros(const u32 N);
  inline void zeros(const u32 n_rows, const u32 n_cols);
  
  inline void ones();
  inline void ones(const u32 N);
  inline void ones(const u32 n_rows, const u32 n_cols);
  
  
  inline void load(const std::string   name, const file_type type = auto_detect, const bool print_status = true);
  inline void load(      std::istream& is,   const file_type type = auto_detect, const bool print_status = true);
  
  inline void quiet_load(const std::string   name, const file_type type = auto_detect);
  inline void quiet_load(      std::istream& is,   const file_type type = auto_detect);
  
  
  typedef       eT*       row_iterator;
  typedef const eT* const_row_iterator;
  
  inline       row_iterator begin_row(const u32 row_num);
  inline const_row_iterator begin_row(const u32 row_num) const;
  
  inline       row_iterator end_row  (const u32 row_num);
  inline const_row_iterator end_row  (const u32 row_num) const;
  };


//! @}
