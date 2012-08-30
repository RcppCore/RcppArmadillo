// Copyright (C) 2011-2012 Ryan Curtin <ryan@igglybob.com>
// Copyright (C) 2011 Matthew Amidon
// Copyright (C) 2012 Conrad Sanderson
//
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose.  You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)

//! \addtogroup SpSubview
//! @{

template<typename eT>
class SpSubview : public SpBase<eT, SpSubview<eT> >
  {
  public:
  
  const SpMat<eT>& m;
  
  typedef eT elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool is_row = false;
  static const bool is_col = false;
  
  const uword aux_row1;
  const uword aux_col1;
  const uword n_rows;
  const uword n_cols;
  const uword n_elem;
  const uword n_nonzero;

  // So that SpValProxy can call add_element() and delete_element().
  friend class SpValProxy<SpSubview<eT> >;

  protected:

  arma_inline SpSubview(const SpMat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols);
  arma_inline SpSubview(      SpMat<eT>& in_m, const uword in_row1, const uword in_col1, const uword in_n_rows, const uword in_n_cols);

  public:

  inline ~SpSubview();

  inline const SpSubview& operator+= (const eT val);
  inline const SpSubview& operator-= (const eT val);
  inline const SpSubview& operator*= (const eT val);
  inline const SpSubview& operator/= (const eT val);

  inline const SpSubview& operator=(const SpSubview& x);

  template<typename T1> inline const SpSubview& operator= (const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator+=(const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator-=(const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator*=(const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator%=(const Base<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator/=(const Base<eT, T1>& x);

  template<typename T1> inline const SpSubview& operator= (const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator+=(const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator-=(const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator*=(const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator%=(const SpBase<eT, T1>& x);
  template<typename T1> inline const SpSubview& operator/=(const SpBase<eT, T1>& x);

  /*
  inline static void extract(SpMat<eT>& out, const SpSubview& in);

  inline static void  plus_inplace(Mat<eT>& out, const subview& in);
  inline static void minus_inplace(Mat<eT>& out, const subview& in);
  inline static void schur_inplace(Mat<eT>& out, const subview& in);
  inline static void   div_inplace(Mat<eT>& out, const subview& in);
  */

  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  inline void eye();

  inline SpValProxy<SpSubview<eT> > operator[](const uword i);
  inline eT                         operator[](const uword i) const;

  inline SpValProxy<SpSubview<eT> > operator()(const uword i);
  inline eT                         operator()(const uword i) const;

  inline SpValProxy<SpSubview<eT> > operator()(const uword in_row, const uword in_col);
  inline eT                         operator()(const uword in_row, const uword in_col) const;

  inline SpValProxy<SpSubview<eT> > at(const uword i);
  inline eT                         at(const uword i) const;

  inline SpValProxy<SpSubview<eT> > at(const uword in_row, const uword in_col);
  inline eT                         at(const uword in_row, const uword in_col) const;

  inline bool check_overlap(const SpSubview& x) const;

  inline bool is_vec() const;

/* not yet
  inline       SpSubview_row<eT> row(const uword row_num);
  inline const SpSubview_row<eT> row(const uword row_num) const;

  inline            SpSubview_row<eT> operator()(const uword row_num, const span& col_span);
  inline      const SpSubview_row<eT> operator()(const uword row_num, const span& col_span) const;

  inline       SpSubview_col<eT> col(const uword col_num);
  inline const SpSubview_col<eT> col(const uword col_num) const;

  inline            SpSubview_col<eT> operator()(const span& row_span, const uword col_num);
  inline      const SpSubview_col<eT> operator()(const span& row_span, const uword col_num) const;

  inline            Col<eT>  unsafe_col(const uword col_num);
  inline      const Col<eT>  unsafe_col(const uword col_num) const;

  inline       SpSubview<eT> rows(const uword in_row1, const uword in_row2);
  inline const SpSubview<eT> rows(const uword in_row1, const uword in_row2) const;

  inline       SpSubview<eT> cols(const uword in_col1, const uword in_col2);
  inline const SpSubview<eT> cols(const uword in_col1, const uword in_col2) const;

  inline       SpSubview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2);
  inline const SpSubview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const;

  inline            SpSubview<eT> submat    (const span& row_span, const span& col_span);
  inline      const SpSubview<eT> submat    (const span& row_span, const span& col_span) const;

  inline            SpSubview<eT> operator()(const span& row_span, const span& col_span);
  inline      const SpSubview<eT> operator()(const span& row_span, const span& col_span) const;

  inline       diagview<eT> diag(const s32 in_id = 0);
  inline const diagview<eT> diag(const s32 in_id = 0) const;
*/

  inline void swap_rows(const uword in_row1, const uword in_row2);
  inline void swap_cols(const uword in_col1, const uword in_col2);

  // Similar to SpMat iterators but automatically iterate over values not in the subview.
  class iterator_base
    {
    public:

    inline iterator_base(const SpSubview& in_M);
    inline iterator_base(const SpSubview& in_M, const uword col, const uword pos, const uword skip_pos);

    inline eT operator*() const;

    inline bool operator!=(const iterator_base& rhs) const;
    inline bool operator==(const iterator_base& rhs) const;

    inline bool operator!=(const typename SpMat<eT>::iterator_base& rhs) const;
    inline bool operator==(const typename SpMat<eT>::iterator_base& rhs) const;

    // Don't hold location internally; call "dummy" methods to get that information.
    arma_inline uword row() const { return M.m.row_indices[internal_pos + skip_pos] - M.aux_row1; }
    arma_inline uword col() const { return internal_col;                                          }
    arma_inline uword pos() const { return internal_pos;                                          }

    arma_aligned const SpSubview& M;
    arma_aligned       uword      internal_col;
    arma_aligned       uword      internal_pos;
    arma_aligned       uword      skip_pos; // not used in row_iterator or const_row_iterator
    };

  class const_iterator : public iterator_base
    {
    public:

    inline const_iterator(const SpSubview& in_M, uword initial_pos = 0);
    inline const_iterator(const SpSubview& in_M, uword in_row, uword in_col);
    inline const_iterator(const const_iterator& other);

    inline const_iterator& operator++();
    inline void            operator++(int);

    inline const_iterator& operator--();
    inline void            operator--(int);
    };

  class iterator : public const_iterator
    {
    public:

    inline iterator(SpSubview& in_M, const uword initial_pos = 0) : const_iterator(in_M, initial_pos) { }
    inline iterator(SpSubview& in_M, const uword in_row, const uword in_col) : const_iterator(in_M, in_row, in_col) { }
    inline iterator(const iterator& other) : const_iterator(other) { }

    inline SpValProxy<SpSubview<eT> > operator*();

    // overloads needed for return type correctness
    inline iterator& operator++();
    inline void      operator++(int);

    inline iterator& operator--();
    inline void      operator--(int);
    };

  class const_row_iterator : public iterator_base
    {
    public:

    inline const_row_iterator(const SpSubview& in_M, uword initial_pos = 0);
    inline const_row_iterator(const SpSubview& in_M, uword in_row, uword in_col);
    inline const_row_iterator(const const_row_iterator& other);

    inline const_row_iterator& operator++();
    inline void                operator++(int);

    inline const_row_iterator& operator--();
    inline void                operator--(int);

    uword internal_row; // Hold row internally because we use internal_pos differently.
    uword actual_pos; // Actual position in subview's parent matrix.

    arma_inline eT operator*() const { return iterator_base::M.m.values[actual_pos]; }

    arma_inline uword row() const { return internal_row; }
    };

  class row_iterator : public const_row_iterator
    {
    public:

    inline row_iterator(SpSubview& in_M, uword initial_pos = 0) : const_row_iterator(in_M, initial_pos) { }
    inline row_iterator(SpSubview& in_M, uword in_row, uword in_col) : const_row_iterator(in_M, in_row, in_col) { }
    inline row_iterator(const row_iterator& other) : const_row_iterator(other) { }

    inline SpValProxy<SpSubview<eT> > operator*();

    // overloads needed for return type correctness
    inline row_iterator& operator++();
    inline void          operator++(int);

    inline row_iterator& operator--();
    inline void          operator--(int);
    };

  inline iterator           begin();
  inline const_iterator     begin() const;

  inline iterator           begin_col(const uword col_num);
  inline const_iterator     begin_col(const uword col_num) const;

  inline row_iterator       begin_row();
  inline const_row_iterator begin_row() const;

  inline iterator           end();
  inline const_iterator     end() const;

  inline row_iterator       end_row();
  inline const_row_iterator end_row() const;


  private:
  friend class SpMat<eT>;
  SpSubview();

  // For use by SpValProxy.  We just update n_nonzero and pass the call on to the matrix.
  inline arma_hot arma_warn_unused eT&  add_element(const uword in_row, const uword in_col, const eT in_val = 0.0);
  inline arma_hot                  void delete_element(const uword in_row, const uword in_col);

  };

/*
template<typename eT>
class SpSubview_col : public SpSubview<eT>
  {
  public:

  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;

  inline void operator= (const SpSubview<eT>& x);
  inline void operator= (const SpSubview_col& x);

  template<typename T1>
  inline void operator= (const Base<eT,T1>& x);

  inline       SpSubview_col<eT> rows(const uword in_row1, const uword in_row2);
  inline const SpSubview_col<eT> rows(const uword in_row1, const uword in_row2) const;

  inline       SpSubview_col<eT> subvec(const uword in_row1, const uword in_row2);
  inline const SpSubview_col<eT> subvec(const uword in_row1, const uword in_row2) const;


  protected:

  inline SpSubview_col(const Mat<eT>& in_m, const uword in_col);
  inline SpSubview_col(      Mat<eT>& in_m, const uword in_col);

  inline SpSubview_col(const Mat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows);
  inline SpSubview_col(      Mat<eT>& in_m, const uword in_col, const uword in_row1, const uword in_n_rows);


  private:

  friend class Mat<eT>;
  friend class Col<eT>;
  friend class SpSubview<eT>;

  SpSubview_col();
  };

template<typename eT>
class SpSubview_row : public SpSubview<eT>
  {
  public:

  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;

  inline void operator= (const SpSubview<eT>& x);
  inline void operator= (const SpSubview_row& x);

  template<typename T1>
  inline void operator= (const Base<eT,T1>& x);

  inline       SpSubview_row<eT> cols(const uword in_col1, const uword in_col2);
  inline const SpSubview_row<eT> cols(const uword in_col1, const uword in_col2) const;

  inline       SpSubview_row<eT> subvec(const uword in_col1, const uword in_col2);
  inline const SpSubview_row<eT> subvec(const uword in_col1, const uword in_col2) const;


  protected:

  inline SpSubview_row(const Mat<eT>& in_m, const uword in_row);
  inline SpSubview_row(      Mat<eT>& in_m, const uword in_row);

  inline SpSubview_row(const Mat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols);
  inline SpSubview_row(      Mat<eT>& in_m, const uword in_row, const uword in_col1, const uword in_n_cols);


  private:

  friend class Mat<eT>;
  friend class Row<eT>;
  friend class SpSubview<eT>;

  SpSubview_row();
  };
*/

//! @}
