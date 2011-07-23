// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// Copyright (C)      2011 James Sanders
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup subview
//! @{


//! Class for storing data required to construct or apply operations to a submatrix
//! (i.e. where the submatrix starts and ends as well as a reference/pointer to the original matrix),
template<typename eT>
class subview : public Base<eT, subview<eT> >
  {
  public:    arma_aligned const Mat<eT>& m;
  protected: arma_aligned       Mat<eT>* m_ptr;
  
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  const u32 aux_row1;
  const u32 aux_col1;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  
  protected:
  
  arma_inline subview(const Mat<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_n_rows, const u32 in_n_cols);
  arma_inline subview(      Mat<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_n_rows, const u32 in_n_cols);
  
  
  public:
  
  inline ~subview();
  
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  // deliberately returning void
  template<typename T1> inline void operator=  (const Base<eT,T1>& x);
  template<typename T1> inline void operator+= (const Base<eT,T1>& x);
  template<typename T1> inline void operator-= (const Base<eT,T1>& x);
  template<typename T1> inline void operator%= (const Base<eT,T1>& x);
  template<typename T1> inline void operator/= (const Base<eT,T1>& x);
  
  inline void operator=  (const subview& x);
  inline void operator+= (const subview& x);
  inline void operator-= (const subview& x);
  inline void operator%= (const subview& x);
  inline void operator/= (const subview& x);
  
  inline static void extract(Mat<eT>& out, const subview& in);
  
  inline static void  plus_inplace(Mat<eT>& out, const subview& in);
  inline static void minus_inplace(Mat<eT>& out, const subview& in);
  inline static void schur_inplace(Mat<eT>& out, const subview& in);
  inline static void   div_inplace(Mat<eT>& out, const subview& in);
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  inline void eye();
  
  inline eT& operator[](const u32 i);
  inline eT  operator[](const u32 i) const;
  
  inline eT& operator()(const u32 i);
  inline eT  operator()(const u32 i) const;
  
  inline eT& operator()(const u32 in_row, const u32 in_col);
  inline eT  operator()(const u32 in_row, const u32 in_col) const;
  
  inline eT&         at(const u32 in_row, const u32 in_col);
  inline eT          at(const u32 in_row, const u32 in_col) const;
  
  arma_inline       eT* colptr(const u32 in_col);
  arma_inline const eT* colptr(const u32 in_col) const;
  
  inline bool check_overlap(const subview& x) const;
  
  inline bool is_vec() const;
  
  inline       subview_row<eT> row(const u32 row_num);
  inline const subview_row<eT> row(const u32 row_num) const;
  
  inline            subview_row<eT> operator()(const u32 row_num, const span& col_span);
  inline      const subview_row<eT> operator()(const u32 row_num, const span& col_span) const;
  
  inline       subview_col<eT> col(const u32 col_num);
  inline const subview_col<eT> col(const u32 col_num) const;
  
  inline            subview_col<eT> operator()(const span& row_span, const u32 col_num);
  inline      const subview_col<eT> operator()(const span& row_span, const u32 col_num) const;
  
  inline            Col<eT>  unsafe_col(const u32 col_num);
  inline      const Col<eT>  unsafe_col(const u32 col_num) const;
  
  inline       subview<eT> rows(const u32 in_row1, const u32 in_row2);
  inline const subview<eT> rows(const u32 in_row1, const u32 in_row2) const;
  
  inline       subview<eT> cols(const u32 in_col1, const u32 in_col2);
  inline const subview<eT> cols(const u32 in_col1, const u32 in_col2) const;
  
  inline       subview<eT> submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2);
  inline const subview<eT> submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2) const;
  
  inline            subview<eT> submat    (const span& row_span, const span& col_span);
  inline      const subview<eT> submat    (const span& row_span, const span& col_span) const;
  
  inline            subview<eT> operator()(const span& row_span, const span& col_span);
  inline      const subview<eT> operator()(const span& row_span, const span& col_span) const;
  
  inline       diagview<eT> diag(const s32 in_id = 0);
  inline const diagview<eT> diag(const s32 in_id = 0) const;
  
  inline void swap_rows(const u32 in_row1, const u32 in_row2);
  inline void swap_cols(const u32 in_col1, const u32 in_col2);
  
  
  inline void print(const std::string extra_text = "") const;
  inline void print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void print_trans(const std::string extra_text = "") const;
  inline void print_trans(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void raw_print(const std::string extra_text = "") const;
  inline void raw_print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void raw_print_trans(const std::string extra_text = "") const;
  inline void raw_print_trans(std::ostream& user_stream, const std::string extra_text = "") const;
  
  
  // // primitive forward iterator
  // class iter
  //   {
  //   public:
  //   
  //   inline iter(const subview<eT>& in_M);
  //   
  //   arma_inline eT operator* () const;
  //   
  //   inline void operator++();
  //   inline void operator++(int);
  //   
  //   
  //   private:
  //   
  //   arma_aligned const eT* mem;
  //   
  //   arma_aligned u32 n_rows;
  //   
  //   arma_aligned u32 row_start;
  //   arma_aligned u32 row_end_p1;
  //   
  //   arma_aligned u32 row;
  //   arma_aligned u32 col;
  //   arma_aligned u32 i;
  //   };
  
  
  private:
  
  friend class Mat<eT>;
  subview();
  };



template<typename eT>
class subview_col : public subview<eT>
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline void operator= (const subview<eT>& x);
  inline void operator= (const subview_col& x);
  
  template<typename T1>
  inline void operator= (const Base<eT,T1>& x);
  
  inline       subview_col<eT> rows(const u32 in_row1, const u32 in_row2);
  inline const subview_col<eT> rows(const u32 in_row1, const u32 in_row2) const;
  
  inline       subview_col<eT> subvec(const u32 in_row1, const u32 in_row2);
  inline const subview_col<eT> subvec(const u32 in_row1, const u32 in_row2) const;
  
  
  protected:
  
  inline subview_col(const Mat<eT>& in_m, const u32 in_col);
  inline subview_col(      Mat<eT>& in_m, const u32 in_col);
  
  inline subview_col(const Mat<eT>& in_m, const u32 in_col, const u32 in_row1, const u32 in_n_rows);
  inline subview_col(      Mat<eT>& in_m, const u32 in_col, const u32 in_row1, const u32 in_n_rows);
  
  
  private:
  
  friend class Mat<eT>;
  friend class Col<eT>;
  friend class subview<eT>;
  
  subview_col();
  };



template<typename eT>
class subview_row : public subview<eT>
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline void operator= (const subview<eT>& x);
  inline void operator= (const subview_row& x);
  
  template<typename T1>
  inline void operator= (const Base<eT,T1>& x);
  
  inline       subview_row<eT> cols(const u32 in_col1, const u32 in_col2);
  inline const subview_row<eT> cols(const u32 in_col1, const u32 in_col2) const;
  
  inline       subview_row<eT> subvec(const u32 in_col1, const u32 in_col2);
  inline const subview_row<eT> subvec(const u32 in_col1, const u32 in_col2) const;
  
  
  protected:
  
  inline subview_row(const Mat<eT>& in_m, const u32 in_row);
  inline subview_row(      Mat<eT>& in_m, const u32 in_row);
  
  inline subview_row(const Mat<eT>& in_m, const u32 in_row, const u32 in_col1, const u32 in_n_cols);
  inline subview_row(      Mat<eT>& in_m, const u32 in_row, const u32 in_col1, const u32 in_n_cols);
  
  
  private:
  
  friend class Mat<eT>;
  friend class Row<eT>;
  friend class subview<eT>;
  
  subview_row();
  };



//! @}
