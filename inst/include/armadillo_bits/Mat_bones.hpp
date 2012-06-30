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


//! \addtogroup Mat
//! @{



//! Dense matrix class

template<typename eT>
class Mat : public Base< eT, Mat<eT> >
  {
  public:
  
  typedef eT                                elem_type;  //!< the type of elements stored in the matrix
  typedef typename get_pod_type<eT>::result pod_type;   //!< if eT is non-complex, pod_type is same as eT. otherwise, pod_type is the underlying type used by std::complex
  
  static const bool is_col = false;
  static const bool is_row = false;
  
  
  const uword  n_rows;    //!< number of rows in the matrix (read-only)
  const uword  n_cols;    //!< number of columns in the matrix (read-only)
  const uword  n_elem;    //!< number of elements in the matrix (read-only)
  const uhword vec_state; //!< 0: matrix layout; 1: column vector layout; 2: row vector layout
  const uhword mem_state;
  
  // mem_state = 0: normal matrix that can be resized; 
  // mem_state = 1: use auxiliary memory until change in the number of elements is requested;  
  // mem_state = 2: use auxiliary memory and don't allow the number of elements to be changed; 
  // mem_state = 3: fixed size (e.g. via template based size specification).
  
  arma_aligned const eT* const mem;  //!< pointer to the memory used by the matrix (memory is read-only)
  
  protected:
  arma_aligned eT mem_local[ arma_config::mat_prealloc ];
  
  
  public:
  
  inline ~Mat();
  inline  Mat();
  
  inline Mat(const uword in_rows, const uword in_cols);
  
  inline                  Mat(const char*        text);
  inline const Mat& operator=(const char*        text);
  
  inline                  Mat(const std::string& text);
  inline const Mat& operator=(const std::string& text);
  
  #if defined(ARMA_USE_CXX11)
  inline                  Mat(const std::initializer_list<eT>& list);
  inline const Mat& operator=(const std::initializer_list<eT>& list);
  #endif
  
  inline Mat(      eT* aux_mem, const uword aux_n_rows, const uword aux_n_cols, const bool copy_aux_mem = true, const bool strict = true);
  inline Mat(const eT* aux_mem, const uword aux_n_rows, const uword aux_n_cols);
  
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
  
  template<typename T1> inline                   Mat(const BaseCube<eT,T1>& X);
  template<typename T1> inline const Mat&  operator=(const BaseCube<eT,T1>& X);
  template<typename T1> inline const Mat& operator+=(const BaseCube<eT,T1>& X);
  template<typename T1> inline const Mat& operator-=(const BaseCube<eT,T1>& X);
  template<typename T1> inline const Mat& operator*=(const BaseCube<eT,T1>& X);
  template<typename T1> inline const Mat& operator%=(const BaseCube<eT,T1>& X);
  template<typename T1> inline const Mat& operator/=(const BaseCube<eT,T1>& X);
  
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
  
  template<typename T1> inline                   Mat(const subview_elem1<eT,T1>& X);
  template<typename T1> inline const Mat& operator= (const subview_elem1<eT,T1>& X);
  template<typename T1> inline const Mat& operator+=(const subview_elem1<eT,T1>& X);
  template<typename T1> inline const Mat& operator-=(const subview_elem1<eT,T1>& X);
  template<typename T1> inline const Mat& operator*=(const subview_elem1<eT,T1>& X);
  template<typename T1> inline const Mat& operator%=(const subview_elem1<eT,T1>& X);
  template<typename T1> inline const Mat& operator/=(const subview_elem1<eT,T1>& X);
  
  template<typename T1, typename T2> inline                   Mat(const subview_elem2<eT,T1,T2>& X);
  template<typename T1, typename T2> inline const Mat& operator= (const subview_elem2<eT,T1,T2>& X);
  template<typename T1, typename T2> inline const Mat& operator+=(const subview_elem2<eT,T1,T2>& X);
  template<typename T1, typename T2> inline const Mat& operator-=(const subview_elem2<eT,T1,T2>& X);
  template<typename T1, typename T2> inline const Mat& operator*=(const subview_elem2<eT,T1,T2>& X);
  template<typename T1, typename T2> inline const Mat& operator%=(const subview_elem2<eT,T1,T2>& X);
  template<typename T1, typename T2> inline const Mat& operator/=(const subview_elem2<eT,T1,T2>& X);
  
  
  inline mat_injector<Mat> operator<<(const eT val);
  inline mat_injector<Mat> operator<<(const injector_end_of_row<>& x);
  
  
  arma_inline       subview_row<eT> row(const uword row_num);
  arma_inline const subview_row<eT> row(const uword row_num) const;
  
  inline            subview_row<eT> operator()(const uword row_num, const span& col_span);
  inline      const subview_row<eT> operator()(const uword row_num, const span& col_span) const;
  
  
  arma_inline       subview_col<eT> col(const uword col_num);
  arma_inline const subview_col<eT> col(const uword col_num) const;
  
  inline            subview_col<eT> operator()(const span& row_span, const uword col_num);
  inline      const subview_col<eT> operator()(const span& row_span, const uword col_num) const;
  
  inline            Col<eT>  unsafe_col(const uword col_num);
  inline      const Col<eT>  unsafe_col(const uword col_num) const;
  
  
  arma_inline       subview<eT> rows(const uword in_row1, const uword in_row2);
  arma_inline const subview<eT> rows(const uword in_row1, const uword in_row2) const;
  
  arma_inline       subview<eT> cols(const uword in_col1, const uword in_col2);
  arma_inline const subview<eT> cols(const uword in_col1, const uword in_col2) const;
  
  arma_inline       subview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2);
  arma_inline const subview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const;
  
  
  inline            subview<eT> submat    (const span& row_span, const span& col_span);
  inline      const subview<eT> submat    (const span& row_span, const span& col_span) const;
  
  inline            subview<eT> operator()(const span& row_span, const span& col_span);
  inline      const subview<eT> operator()(const span& row_span, const span& col_span) const;
  
  
  template<typename T1> arma_inline       subview_elem1<eT,T1> elem(const Base<uword,T1>& a);
  template<typename T1> arma_inline const subview_elem1<eT,T1> elem(const Base<uword,T1>& a) const;
  
  
  template<typename T1, typename T2> arma_inline       subview_elem2<eT,T1,T2> elem(const Base<uword,T1>& ri, const Base<uword,T2>& ci);
  template<typename T1, typename T2> arma_inline const subview_elem2<eT,T1,T2> elem(const Base<uword,T1>& ri, const Base<uword,T2>& ci) const;
  
  template<typename T1, typename T2> arma_inline       subview_elem2<eT,T1,T2> submat(const Base<uword,T1>& ri, const Base<uword,T2>& ci);
  template<typename T1, typename T2> arma_inline const subview_elem2<eT,T1,T2> submat(const Base<uword,T1>& ri, const Base<uword,T2>& ci) const;
  
  template<typename T1, typename T2> arma_inline       subview_elem2<eT,T1,T2> operator()(const Base<uword,T1>& ri, const Base<uword,T2>& ci);
  template<typename T1, typename T2> arma_inline const subview_elem2<eT,T1,T2> operator()(const Base<uword,T1>& ri, const Base<uword,T2>& ci) const;
  
  template<typename T1> arma_inline       subview_elem2<eT,T1,T1> rows(const Base<uword,T1>& ri);
  template<typename T1> arma_inline const subview_elem2<eT,T1,T1> rows(const Base<uword,T1>& ri) const;
  
  template<typename T2> arma_inline       subview_elem2<eT,T2,T2> cols(const Base<uword,T2>& ci);
  template<typename T2> arma_inline const subview_elem2<eT,T2,T2> cols(const Base<uword,T2>& ci) const;
  
  
  
  arma_inline       diagview<eT> diag(const sword in_id = 0);
  arma_inline const diagview<eT> diag(const sword in_id = 0) const;
  
  
  inline void swap_rows(const uword in_row1, const uword in_row2);
  inline void swap_cols(const uword in_col1, const uword in_col2);
  
  inline void shed_row(const uword row_num);
  inline void shed_col(const uword col_num);
  
  inline void shed_rows(const uword in_row1, const uword in_row2);
  inline void shed_cols(const uword in_col1, const uword in_col2);
  
  inline void insert_rows(const uword row_num, const uword N, const bool set_to_zero = true);
  inline void insert_cols(const uword col_num, const uword N, const bool set_to_zero = true);
  
  template<typename T1> inline void insert_rows(const uword row_num, const Base<eT,T1>& X);
  template<typename T1> inline void insert_cols(const uword col_num, const Base<eT,T1>& X);
  
  
  template<typename T1, typename gen_type> inline                   Mat(const Gen<T1, gen_type>& X);
  template<typename T1, typename gen_type> inline const Mat&  operator=(const Gen<T1, gen_type>& X);
  template<typename T1, typename gen_type> inline const Mat& operator+=(const Gen<T1, gen_type>& X);
  template<typename T1, typename gen_type> inline const Mat& operator-=(const Gen<T1, gen_type>& X);
  template<typename T1, typename gen_type> inline const Mat& operator*=(const Gen<T1, gen_type>& X);
  template<typename T1, typename gen_type> inline const Mat& operator%=(const Gen<T1, gen_type>& X);
  template<typename T1, typename gen_type> inline const Mat& operator/=(const Gen<T1, gen_type>& X);
  
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
  
  
  arma_inline arma_warn_unused eT& operator[] (const uword ii);
  arma_inline arma_warn_unused eT  operator[] (const uword ii) const;
  arma_inline arma_warn_unused eT& at         (const uword ii);
  arma_inline arma_warn_unused eT  at         (const uword ii) const;
  arma_inline arma_warn_unused eT& operator() (const uword ii);
  arma_inline arma_warn_unused eT  operator() (const uword ii) const;
  
  arma_inline arma_warn_unused eT& at         (const uword in_row, const uword in_col);
  arma_inline arma_warn_unused eT  at         (const uword in_row, const uword in_col) const;
  arma_inline arma_warn_unused eT& operator() (const uword in_row, const uword in_col);
  arma_inline arma_warn_unused eT  operator() (const uword in_row, const uword in_col) const;
  
  arma_inline const Mat& operator++();
  arma_inline void       operator++(int);
  
  arma_inline const Mat& operator--();
  arma_inline void       operator--(int);
  
  arma_inline arma_warn_unused bool is_empty()  const;
  arma_inline arma_warn_unused bool is_vec()    const;
  arma_inline arma_warn_unused bool is_rowvec() const;
  arma_inline arma_warn_unused bool is_colvec() const;
  arma_inline arma_warn_unused bool is_square() const;
       inline arma_warn_unused bool is_finite() const;
  
  arma_inline arma_warn_unused bool in_range(const uword ii) const;
  arma_inline arma_warn_unused bool in_range(const span& x ) const;
  
  arma_inline arma_warn_unused bool in_range(const uword   in_row, const uword   in_col) const;
  arma_inline arma_warn_unused bool in_range(const span& row_span, const uword   in_col) const;
  arma_inline arma_warn_unused bool in_range(const uword   in_row, const span& col_span) const;
  arma_inline arma_warn_unused bool in_range(const span& row_span, const span& col_span) const;
  
  arma_inline arma_warn_unused       eT* colptr(const uword in_col);
  arma_inline arma_warn_unused const eT* colptr(const uword in_col) const;
  
  arma_inline arma_warn_unused       eT* memptr();
  arma_inline arma_warn_unused const eT* memptr() const;
  
  
  inline void impl_print(const std::string& extra_text) const;
  inline void impl_print(std::ostream& user_stream, const std::string& extra_text) const;
  
  inline void impl_raw_print(const std::string& extra_text) const;
  inline void impl_raw_print(std::ostream& user_stream, const std::string& extra_text) const;
  
  
  template<typename eT2, typename expr>
  inline void copy_size(const Base<eT2,expr>& X);
  
  inline void set_size(const uword in_elem);
  inline void set_size(const uword in_rows, const uword in_cols);
  
  inline void   resize(const uword in_elem);
  inline void   resize(const uword in_rows, const uword in_cols);
  inline void  reshape(const uword in_rows, const uword in_cols, const uword dim = 0);
  
  
  arma_hot inline const Mat& fill(const eT val);
  
  inline const Mat& zeros();
  inline const Mat& zeros(const uword in_elem);
  inline const Mat& zeros(const uword in_rows, const uword in_cols);
  
  inline const Mat& ones();
  inline const Mat& ones(const uword in_elem);
  inline const Mat& ones(const uword in_rows, const uword in_cols);
  
  inline const Mat& randu();
  inline const Mat& randu(const uword in_elem);
  inline const Mat& randu(const uword in_rows, const uword in_cols);
  
  inline const Mat& randn();
  inline const Mat& randn(const uword in_elem);
  inline const Mat& randn(const uword in_rows, const uword in_cols);
  
  inline const Mat& eye();
  inline const Mat& eye(const uword in_rows, const uword in_cols);
  
  inline void reset();
  
  
  template<typename T1> inline void set_real(const Base<pod_type,T1>& X);
  template<typename T1> inline void set_imag(const Base<pod_type,T1>& X);
  
  
  inline arma_warn_unused eT min() const;
  inline arma_warn_unused eT max() const;
  
  inline eT min(uword& index_of_min_val) const;
  inline eT max(uword& index_of_max_val) const;
  
  inline eT min(uword& row_of_min_val, uword& col_of_min_val) const;
  inline eT max(uword& row_of_max_val, uword& col_of_max_val) const;
  
  
  inline bool save(const std::string   name, const file_type type = arma_binary, const bool print_status = true) const;
  inline bool save(      std::ostream& os,   const file_type type = arma_binary, const bool print_status = true) const;
  
  inline bool load(const std::string   name, const file_type type = auto_detect, const bool print_status = true);
  inline bool load(      std::istream& is,   const file_type type = auto_detect, const bool print_status = true);
  
  inline bool quiet_save(const std::string   name, const file_type type = arma_binary) const;
  inline bool quiet_save(      std::ostream& os,   const file_type type = arma_binary) const;
  
  inline bool quiet_load(const std::string   name, const file_type type = auto_detect);
  inline bool quiet_load(      std::istream& is,   const file_type type = auto_detect);
  
  
  // for container-like functionality
  
  typedef eT    value_type;
  typedef uword size_type;
  
  typedef       eT*       iterator;
  typedef const eT* const_iterator;
  
  typedef       eT*       col_iterator;
  typedef const eT* const_col_iterator;
  
  class row_iterator
    {
    public:
    
    inline row_iterator(Mat<eT>& in_M, const uword in_row);
    
    inline eT& operator* ();
    
    inline row_iterator& operator++();
    inline void          operator++(int);
    
    inline row_iterator& operator--();
    inline void          operator--(int);
    
    inline bool operator!=(const row_iterator& X) const;
    inline bool operator==(const row_iterator& X) const;
    
    arma_aligned Mat<eT>& M;
    arma_aligned uword    row;
    arma_aligned uword    col;
    };
  
  
  class const_row_iterator
    {
    public:
    
    const_row_iterator(const Mat<eT>& in_M, const uword in_row);
    const_row_iterator(const row_iterator& X);
    
    inline eT operator*() const;
    
    inline const_row_iterator& operator++();
    inline void                operator++(int);
    
    inline const_row_iterator& operator--();
    inline void                operator--(int);
    
    inline bool operator!=(const const_row_iterator& X) const;
    inline bool operator==(const const_row_iterator& X) const;
    
    arma_aligned const Mat<eT>& M;
    arma_aligned       uword    row;
    arma_aligned       uword    col;
    };
  
  inline       iterator begin();
  inline const_iterator begin() const;
  
  inline       iterator end();
  inline const_iterator end()   const;
  
  inline       col_iterator begin_col(const uword col_num);
  inline const_col_iterator begin_col(const uword col_num) const;
  
  inline       col_iterator end_col  (const uword col_num);
  inline const_col_iterator end_col  (const uword col_num) const;
  
  inline       row_iterator begin_row(const uword row_num);
  inline const_row_iterator begin_row(const uword row_num) const;
  
  inline       row_iterator end_row  (const uword row_num);
  inline const_row_iterator end_row  (const uword row_num) const;
  
  inline void  clear();
  inline bool  empty() const;
  inline uword size()  const;
  
  inline void steal_mem(Mat& X);  //!< don't use this unless you're writing code internal to Armadillo
  
  template<uword fixed_n_rows, uword fixed_n_cols> class fixed;
  
  
  protected:
  
  inline void init_cold();
  inline void init_warm(uword in_rows, uword in_cols);
  
  inline void init(const std::string& text);
  
  #if defined(ARMA_USE_CXX11)
  inline void init(const std::initializer_list<eT>& list);
  #endif
  
  template<typename T1, typename T2>
  inline void init(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);
  
  inline Mat(const char junk, const eT* aux_mem, const uword aux_n_rows, const uword aux_n_cols);
  
  inline Mat(const arma_vec_indicator&, const uhword in_vec_state);
  inline Mat(const arma_vec_indicator&, const uword in_n_rows, const uword in_n_cols, const uhword in_vec_state);
  
  inline Mat(const arma_fixed_indicator&, const uword in_n_rows, const uword in_n_cols, const uhword in_vec_state, const eT* in_mem);
  
  
  friend class Cube<eT>;
  friend class glue_join;
  friend class op_strans;
  friend class op_htrans;
  friend class op_resize;
  
  public:
  
  #ifdef ARMA_EXTRA_MAT_PROTO
    #include ARMA_INCFILE_WRAP(ARMA_EXTRA_MAT_PROTO)
  #endif
  };



template<typename eT>
template<uword fixed_n_rows, uword fixed_n_cols> 
class Mat<eT>::fixed : public Mat<eT>
  {
  private:
  
  static const uword fixed_n_elem = fixed_n_rows * fixed_n_cols;
  static const bool  use_extra    = (fixed_n_elem > arma_config::mat_prealloc);
  
  arma_aligned eT mem_local_extra[ (use_extra) ? fixed_n_elem : 1 ];
  
  
  public:
  
  typedef fixed<fixed_n_rows, fixed_n_cols> Mat_fixed_type;
  
  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  
  static const bool is_col = (fixed_n_cols == 1) ? true : false;
  static const bool is_row = (fixed_n_rows == 1) ? true : false;

  static const uword n_rows = fixed_n_rows;
  static const uword n_cols = fixed_n_cols;
  static const uword n_elem = fixed_n_elem;
  
  arma_inline fixed();
  arma_inline fixed(const fixed<fixed_n_rows, fixed_n_cols>& X);
  
  template<typename T1>              inline fixed(const Base<eT,T1>& A);
  template<typename T1, typename T2> inline fixed(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);
  
  inline fixed(const eT* aux_mem);
  
  inline fixed(const char*        text);
  inline fixed(const std::string& text);
  
  using Mat<eT>::operator=;
  using Mat<eT>::operator();
  
  #if defined(ARMA_USE_CXX11)
  inline                fixed(const std::initializer_list<eT>& list);
  inline const Mat& operator=(const std::initializer_list<eT>& list);
  #endif
  
  
  arma_inline arma_warn_unused eT& operator[] (const uword i);
  arma_inline arma_warn_unused eT  operator[] (const uword i) const;
  arma_inline arma_warn_unused eT& at         (const uword i);
  arma_inline arma_warn_unused eT  at         (const uword i) const;
  arma_inline arma_warn_unused eT& operator() (const uword i);
  arma_inline arma_warn_unused eT  operator() (const uword i) const;
  
  arma_inline arma_warn_unused eT& at         (const uword in_row, const uword in_col);
  arma_inline arma_warn_unused eT  at         (const uword in_row, const uword in_col) const;
  arma_inline arma_warn_unused eT& operator() (const uword in_row, const uword in_col);
  arma_inline arma_warn_unused eT  operator() (const uword in_row, const uword in_col) const;
  
  arma_inline arma_warn_unused       eT* memptr();
  arma_inline arma_warn_unused const eT* memptr() const;
  
  arma_hot inline const Mat<eT>& fill(const eT val);
  arma_hot inline const Mat<eT>& zeros();
  arma_hot inline const Mat<eT>& ones();
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
  
  template<typename eT, typename T1> inline static void set_real(Mat<eT>&                out, const Base<eT,T1>& X);
  template<typename T,  typename T1> inline static void set_real(Mat< std::complex<T> >& out, const Base< T,T1>& X);
  
  template<typename eT, typename T1> inline static void set_imag(Mat<eT>&                out, const Base<eT,T1>& X);
  template<typename T,  typename T1> inline static void set_imag(Mat< std::complex<T> >& out, const Base< T,T1>& X);
  };



//! @}
