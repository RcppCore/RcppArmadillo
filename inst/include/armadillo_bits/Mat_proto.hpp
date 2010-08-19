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


//! \addtogroup Mat
//! @{

//! Dense matrix class

template<typename eT>
class Mat : public Base< eT, Mat<eT> >
  {
  public:
  
  typedef eT                                elem_type;  //!< the type of elements stored in the matrix
  typedef typename get_pod_type<eT>::result pod_type;   //!< if eT is non-complex, pod_type is same as eT. otherwise, pod_type is the underlying type used by std::complex
  
  const u32  n_rows;      //!< number of rows in the matrix (read-only)
  const u32  n_cols;      //!< number of columns in the matrix (read-only)
  const u32  n_elem;      //!< number of elements in the matrix (read-only)
  const bool use_aux_mem; //!< true if externally managed memory is being used (read-only)
  
  arma_aligned const eT* const mem;  //!< pointer to the memory used by the matrix (memory is read-only)
  
  protected:
  arma_aligned eT mem_local[ 16 ];
  
  
  public:
  
  inline ~Mat();
  inline  Mat();
  
  inline Mat(const u32 in_rows, const u32 in_cols);
  
  inline                  Mat(const char*        text);
  inline const Mat& operator=(const char*        text);
  inline                  Mat(const std::string& text);
  inline const Mat& operator=(const std::string& text);
  
  inline Mat(      eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const bool copy_aux_mem = true);
  inline Mat(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols);
  
  arma_inline const Mat&  operator=(const eT val);
  arma_inline const Mat& operator+=(const eT val);
  arma_inline const Mat& operator-=(const eT val);
  arma_inline const Mat& operator*=(const eT val);
  arma_inline const Mat& operator/=(const eT val);
  
  inline                   Mat(const Mat& m);
  inline const Mat&  operator=(const Mat& m);
  inline const Mat& operator+=(const Mat& m);
  inline const Mat& operator-=(const Mat& m);
  inline const Mat& operator*=(const Mat& m);
  inline const Mat& operator%=(const Mat& m);
  inline const Mat& operator/=(const Mat& m);

  template<typename T1, typename T2>
  inline explicit Mat(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);

  inline                   Mat(const subview<eT>& X);
  inline const Mat&  operator=(const subview<eT>& X);
  inline const Mat& operator+=(const subview<eT>& X);
  inline const Mat& operator-=(const subview<eT>& X);
  inline const Mat& operator*=(const subview<eT>& X);
  inline const Mat& operator%=(const subview<eT>& X);
  inline const Mat& operator/=(const subview<eT>& X);

  //inline explicit          Mat(const subview_cube<eT>& X);
  inline                   Mat(const subview_cube<eT>& X);
  inline const Mat&  operator=(const subview_cube<eT>& X);
  inline const Mat& operator+=(const subview_cube<eT>& X);
  inline const Mat& operator-=(const subview_cube<eT>& X);
  inline const Mat& operator*=(const subview_cube<eT>& X);
  inline const Mat& operator%=(const subview_cube<eT>& X);
  inline const Mat& operator/=(const subview_cube<eT>& X);

  
  //inline explicit          Mat(const diagview<eT>& X);
  inline                   Mat(const diagview<eT>& X);
  inline const Mat&  operator=(const diagview<eT>& X);
  inline const Mat& operator+=(const diagview<eT>& X);
  inline const Mat& operator-=(const diagview<eT>& X);
  inline const Mat& operator*=(const diagview<eT>& X);
  inline const Mat& operator%=(const diagview<eT>& X);
  inline const Mat& operator/=(const diagview<eT>& X);
  
  
  inline mat_injector<Mat> operator<<(const eT val);
  inline mat_injector<Mat> operator<<(const injector_helper x);
  
  
  arma_inline       subview_row<eT> row(const u32 row_num);
  arma_inline const subview_row<eT> row(const u32 row_num) const;
  
  arma_inline       subview_col<eT> col(const u32 col_num);
  arma_inline const subview_col<eT> col(const u32 col_num) const;
  
  arma_inline       subview<eT> rows(const u32 in_row1, const u32 in_row2);
  arma_inline const subview<eT> rows(const u32 in_row1, const u32 in_row2) const;
  
  arma_inline       subview<eT> cols(const u32 in_col1, const u32 in_col2);
  arma_inline const subview<eT> cols(const u32 in_col1, const u32 in_col2) const;
  
  arma_inline       subview<eT> submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2);
  arma_inline const subview<eT> submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2) const;

  arma_inline       diagview<eT> diag(const s32 in_id = 0);
  arma_inline const diagview<eT> diag(const s32 in_id = 0) const;
    
  inline void swap_rows(const u32 in_row1, const u32 in_row2);
  inline void swap_cols(const u32 in_col1, const u32 in_col2);
  
  template<typename T1, typename op_type> inline                   Mat(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat&  operator=(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator+=(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator-=(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator*=(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator%=(const Op<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator/=(const Op<T1, op_type>& X);
  
  template<typename T1, typename eop_type> inline                   Mat(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Mat&  operator=(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Mat& operator+=(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Mat& operator-=(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Mat& operator*=(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Mat& operator%=(const eOp<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Mat& operator/=(const eOp<T1, eop_type>& X);
  
  template<typename T1, typename op_type> inline                   Mat(const mtOp<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat&  operator=(const mtOp<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator+=(const mtOp<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator-=(const mtOp<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator*=(const mtOp<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator%=(const mtOp<eT, T1, op_type>& X);
  template<typename T1, typename op_type> inline const Mat& operator/=(const mtOp<eT, T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline                   Mat(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat&  operator=(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator+=(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator-=(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator*=(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator%=(const Glue<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator/=(const Glue<T1, T2, glue_type>& X);
  
  template<typename T1, typename T2>                     inline const Mat& operator+=(const Glue<T1, T2, glue_times>& X);
  template<typename T1, typename T2>                     inline const Mat& operator-=(const Glue<T1, T2, glue_times>& X);
  
  template<typename T1, typename T2, typename eglue_type> inline                   Mat(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Mat&  operator=(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Mat& operator+=(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Mat& operator-=(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Mat& operator*=(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Mat& operator%=(const eGlue<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Mat& operator/=(const eGlue<T1, T2, eglue_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline                   Mat(const mtGlue<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat&  operator=(const mtGlue<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator+=(const mtGlue<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator-=(const mtGlue<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator*=(const mtGlue<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator%=(const mtGlue<eT, T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Mat& operator/=(const mtGlue<eT, T1, T2, glue_type>& X);
  
#ifdef ARMA_EXTRA_MAT_PROTO
#include ARMA_EXTRA_MAT_PROTO
#endif
  
  arma_inline eT& operator[] (const u32 i);
  arma_inline eT  operator[] (const u32 i) const;
  arma_inline eT& operator() (const u32 i);
  arma_inline eT  operator() (const u32 i) const;
  
  arma_inline eT& at         (const u32 in_row, const u32 in_col);
  arma_inline eT  at         (const u32 in_row, const u32 in_col) const;
  arma_inline eT& operator() (const u32 in_row, const u32 in_col);
  arma_inline eT  operator() (const u32 in_row, const u32 in_col) const;
  
  arma_inline const Mat& operator++();
  arma_inline void       operator++(int);
  
  arma_inline const Mat& operator--();
  arma_inline void       operator--(int);
  
  arma_inline bool is_vec() const;
  arma_inline bool is_square() const;
  arma_inline bool is_finite() const;
  
  arma_inline       eT* colptr(const u32 in_col);
  arma_inline const eT* colptr(const u32 in_col) const;
  
  arma_inline       eT* memptr();
  arma_inline const eT* memptr() const;
  
  inline void print(const std::string extra_text = "") const;
  inline void print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void print_trans(const std::string extra_text = "") const;
  inline void print_trans(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void raw_print(const std::string extra_text = "") const;
  inline void raw_print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void raw_print_trans(const std::string extra_text = "") const;
  inline void raw_print_trans(std::ostream& user_stream, const std::string extra_text = "") const;
  
  template<typename eT2>
  inline void copy_size(const Mat<eT2>& m);
  
  inline void  set_size(const u32 in_rows, const u32 in_cols);
  inline void   reshape(const u32 in_rows, const u32 in_cols, const u32 dim = 0);
  
  arma_hot inline void fill(const eT val);
  
  inline void zeros();
  inline void zeros(const u32 in_rows, const u32 in_cols);
  
  inline void ones();
  inline void ones(const u32 in_rows, const u32 in_cols);
  
  inline void reset();
  
  
  inline bool save(const std::string   name, const file_type type = arma_binary, const bool print_status = true) const;
  inline bool save(      std::ostream& os,   const file_type type = arma_binary, const bool print_status = true) const;
  
  inline bool load(const std::string   name, const file_type type = auto_detect, const bool print_status = true);
  inline bool load(      std::istream& is,   const file_type type = auto_detect, const bool print_status = true);
  
  inline bool quiet_save(const std::string   name, const file_type type = arma_binary) const;
  inline bool quiet_save(      std::ostream& os,   const file_type type = arma_binary) const;
  
  inline bool quiet_load(const std::string   name, const file_type type = auto_detect);
  inline bool quiet_load(      std::istream& is,   const file_type type = auto_detect);
  
  
  // iterators
  
  typedef       eT*       iterator;
  typedef const eT* const_iterator;
  
  typedef       eT*       col_iterator;
  typedef const eT* const_col_iterator;
  
  class row_iterator
    {
    public:
    
    inline row_iterator(Mat<eT>& in_M, const u32 in_row);
    
    inline eT& operator* ();
    
    inline row_iterator& operator++();
    inline void          operator++(int);
    
    inline row_iterator& operator--();
    inline void          operator--(int);
    
    inline bool operator!=(const row_iterator& X) const;
    inline bool operator==(const row_iterator& X) const;
    
    arma_aligned Mat<eT>& M;
    arma_aligned u32      row;
    arma_aligned u32      col;
    };
  
  
  class const_row_iterator
    {
    public:
    
    const_row_iterator(const Mat<eT>& in_M, const u32 in_row);
    const_row_iterator(const row_iterator& X);
    
    inline eT operator*() const;
    
    inline const_row_iterator& operator++();
    inline void                operator++(int);
    
    inline const_row_iterator& operator--();
    inline void                operator--(int);
    
    inline bool operator!=(const const_row_iterator& X) const;
    inline bool operator==(const const_row_iterator& X) const;
    
    arma_aligned const Mat<eT>& M;
    arma_aligned       u32      row;
    arma_aligned       u32      col;
    };
  
  inline       iterator begin();
  inline const_iterator begin() const;
  
  inline       iterator end();
  inline const_iterator end()   const;
  
  inline       col_iterator begin_col(const u32 col_num);
  inline const_col_iterator begin_col(const u32 col_num) const;
  
  inline       col_iterator end_col  (const u32 col_num);
  inline const_col_iterator end_col  (const u32 col_num) const;
  
  inline       row_iterator begin_row(const u32 row_num);
  inline const_row_iterator begin_row(const u32 row_num) const;
  
  inline       row_iterator end_row  (const u32 row_num);
  inline const_row_iterator end_row  (const u32 row_num) const;
  
  
  protected:
  
  inline void init(const u32 in_rows, const u32 in_cols);
  inline void init(const std::string& text);
  inline void init(const Mat& x);
  
  inline Mat(const char junk, const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols);
  
  friend class Cube<eT>;
  friend class glue_join;
  };



class Mat_aux
  {
  public:

  template<typename eT> arma_inline static void prefix_pp(Mat<eT>& x);
  template<typename T>  arma_inline static void prefix_pp(Mat< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void postfix_pp(Mat<eT>& x);
  template<typename T>  arma_inline static void postfix_pp(Mat< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void prefix_mm(Mat<eT>& x);
  template<typename T>  arma_inline static void prefix_mm(Mat< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void postfix_mm(Mat<eT>& x);
  template<typename T>  arma_inline static void postfix_mm(Mat< std::complex<T> >& x);
  };



//! @}
