// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// Copyright (C)      2011 James Sanders
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup subview
//! @{


//! Class for storing data required to construct or apply operations to a submatrix
//! (i.e. where the submatrix starts and ends as well as a reference/pointer to the original matrix),
template<typename eT>
class subview : public Base<eT, subview<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  arma_aligned const Mat<eT>& m;
  
  static const bool is_row = false;
  static const bool is_col = false;
  
  const uword aux_row1;
  const uword aux_col1;
  
  const uword n_rows;
  const uword n_cols;
  const uword n_elem;
  
  protected:
  
  arma_inline subview(const Mat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols);
  
  
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
  
  inline eT& operator[](const uword ii);
  inline eT  operator[](const uword ii) const;
  
  inline eT& operator()(const uword ii);
  inline eT  operator()(const uword ii) const;
  
  inline eT& operator()(const uword in_row, const uword in_col);
  inline eT  operator()(const uword in_row, const uword in_col) const;
  
  inline eT&         at(const uword in_row, const uword in_col);
  inline eT          at(const uword in_row, const uword in_col) const;
  
  arma_inline       eT* colptr(const uword in_col);
  arma_inline const eT* colptr(const uword in_col) const;
  
  inline bool check_overlap(const subview& x) const;
  
  inline bool is_vec() const;
  
  inline       subview_row<eT> row(const uword row_num);
  inline const subview_row<eT> row(const uword row_num) const;
  
  inline            subview_row<eT> operator()(const uword row_num, const span& col_span);
  inline      const subview_row<eT> operator()(const uword row_num, const span& col_span) const;
  
  inline       subview_col<eT> col(const uword col_num);
  inline const subview_col<eT> col(const uword col_num) const;
  
  inline            subview_col<eT> operator()(const span& row_span, const uword col_num);
  inline      const subview_col<eT> operator()(const span& row_span, const uword col_num) const;
  
  inline            Col<eT>  unsafe_col(const uword col_num);
  inline      const Col<eT>  unsafe_col(const uword col_num) const;
  
  inline       subview<eT> rows(const uword in_row1, const uword in_row2);
  inline const subview<eT> rows(const uword in_row1, const uword in_row2) const;
  
  inline       subview<eT> cols(const uword in_col1, const uword in_col2);
  inline const subview<eT> cols(const uword in_col1, const uword in_col2) const;
  
  inline       subview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2);
  inline const subview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const;
  
  inline            subview<eT> submat    (const span& row_span, const span& col_span);
  inline      const subview<eT> submat    (const span& row_span, const span& col_span) const;
  
  inline            subview<eT> operator()(const span& row_span, const span& col_span);
  inline      const subview<eT> operator()(const span& row_span, const span& col_span) const;
  
  inline subview_each1< subview<eT>, 0 > each_col();
  inline subview_each1< subview<eT>, 1 > each_row();
  
  template<typename T1> inline subview_each2< subview<eT>, 0, T1 > each_col(const Base<uword, T1>& indices);
  template<typename T1> inline subview_each2< subview<eT>, 1, T1 > each_row(const Base<uword, T1>& indices);
  
  inline       diagview<eT> diag(const sword in_id = 0);
  inline const diagview<eT> diag(const sword in_id = 0) const;
  
  inline void swap_rows(const uword in_row1, const uword in_row2);
  inline void swap_cols(const uword in_col1, const uword in_col2);
  
  
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
  //   arma_aligned uword n_rows;
  //   
  //   arma_aligned uword row_start;
  //   arma_aligned uword row_end_p1;
  //   
  //   arma_aligned uword row;
  //   arma_aligned uword col;
  //   arma_aligned uword i;
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
  
  static const bool is_row = false;
  static const bool is_col = true;
  
  const eT* colmem;
  
  inline void operator= (const subview<eT>& x);
  inline void operator= (const subview_col& x);
  
  template<typename T1>
  inline void operator= (const Base<eT,T1>& x);
  
  arma_inline const Op<subview_col<eT>,op_htrans>  t() const;
  arma_inline const Op<subview_col<eT>,op_htrans> ht() const;
  arma_inline const Op<subview_col<eT>,op_strans> st() const;
  
  inline eT& operator[](const uword i);
  inline eT  operator[](const uword i) const;
  
  inline eT& operator()(const uword i);
  inline eT  operator()(const uword i) const;
  
  inline eT& operator()(const uword in_row, const uword in_col);
  inline eT  operator()(const uword in_row, const uword in_col) const;
  
  inline eT&         at(const uword in_row, const uword in_col);
  inline eT          at(const uword in_row, const uword in_col) const;
  
  arma_inline       eT* colptr(const uword in_col);
  arma_inline const eT* colptr(const uword in_col) const;
  
  inline       subview_col<eT> rows(const uword in_row1, const uword in_row2);
  inline const subview_col<eT> rows(const uword in_row1, const uword in_row2) const;
  
  inline       subview_col<eT> subvec(const uword in_row1, const uword in_row2);
  inline const subview_col<eT> subvec(const uword in_row1, const uword in_row2) const;
  
  // TODO: add operator()(span)
  
  protected:
  
  inline subview_col(const Mat<eT>& in_m, const uword in_col);
  inline subview_col(const Mat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows);
  
  
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
  
  static const bool is_row = true;
  static const bool is_col = false;
  
  inline void operator= (const subview<eT>& x);
  inline void operator= (const subview_row& x);
  
  template<typename T1>
  inline void operator= (const Base<eT,T1>& x);
  
  arma_inline const Op<subview_row<eT>,op_htrans>  t() const;
  arma_inline const Op<subview_row<eT>,op_htrans> ht() const;
  arma_inline const Op<subview_row<eT>,op_strans> st() const;
  
  inline eT& operator[](const uword i);
  inline eT  operator[](const uword i) const;
  
  inline eT& operator()(const uword i);
  inline eT  operator()(const uword i) const;
  
  inline eT& operator()(const uword in_row, const uword in_col);
  inline eT  operator()(const uword in_row, const uword in_col) const;
  
  inline eT&         at(const uword in_row, const uword in_col);
  inline eT          at(const uword in_row, const uword in_col) const;
  
  inline       subview_row<eT> cols(const uword in_col1, const uword in_col2);
  inline const subview_row<eT> cols(const uword in_col1, const uword in_col2) const;
  
  inline       subview_row<eT> subvec(const uword in_col1, const uword in_col2);
  inline const subview_row<eT> subvec(const uword in_col1, const uword in_col2) const;
  
  // TODO: add operator()(span)
  
  protected:
  
  inline subview_row(const Mat<eT>& in_m, const uword in_row);
  inline subview_row(const Mat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols);
  
  
  private:
  
  friend class Mat<eT>;
  friend class Row<eT>;
  friend class subview<eT>;
  
  subview_row();
  };



//! @}
