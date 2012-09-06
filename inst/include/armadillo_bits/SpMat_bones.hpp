// Copyright (C) 2011-2012 Ryan Curtin
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

//! \addtogroup SpMat
//! @{

//! Sparse matrix class, with data stored in compressed sparse column (CSC) format

template<typename eT>
class SpMat : public SpBase< eT, SpMat<eT> >
  {
  public:
  
  typedef eT                                elem_type;  //!< the type of elements stored in the matrix
  typedef typename get_pod_type<eT>::result pod_type;   //!< if eT is non-complex, pod_type is the same as eT; otherwise, pod_type is the underlying type used by std::complex
  
  static const bool is_row = false;
  static const bool is_col = false;
  
  const uword n_rows;    //!< number of rows in the matrix (read-only)
  const uword n_cols;    //!< number of columns in the matrix (read-only)
  const uword n_elem;    //!< number of elements in the matrix (read-only)
  const uword n_nonzero; //!< number of nonzero elements in the matrix (read-only)
  const uword vec_state; //!< 0: matrix; 1: column vector; 2: row vector
  
  // So that SpValProxy can call add_element() and delete_element().
  friend class SpValProxy<SpMat<eT> >;
  friend class SpSubview<eT>;
  
  /**
   * The memory used to store the values of the matrix.
   * In accordance with the CSC format, this stores only the actual values.
   * The correct locations of the values are assembled from the row indices
   * and the column pointers.
   */
  const eT* const values;
  
  /**
   * The row indices of each value.  row_indices[i] is the row of values[i].
   */
  const uword* const row_indices;
  
  /**
   * The column pointers.  This stores the index of the first item in column i.
   * That is, values[col_ptrs[i]] is the first value in column i, and it is in
   * row row_indices[col_ptrs[i]].
   */
  const uword* const col_ptrs;
  
  inline  SpMat();  //! Size will be 0x0 (empty).
  inline ~SpMat();
  
  inline SpMat(const uword in_rows, const uword in_cols);
  
  inline                  SpMat(const char*        text);
  inline const SpMat& operator=(const char*        text);
  inline                  SpMat(const std::string& text);
  inline const SpMat& operator=(const std::string& text);
  inline                  SpMat(const SpMat<eT>&   x);
  
  
  inline const SpMat&  operator=(const eT val); //! Sets size to 1x1.
  inline const SpMat& operator*=(const eT val);
  inline const SpMat& operator/=(const eT val);
  // operator+=(val) and operator-=(val) are not defined as they don't make sense for sparse matrices
  
  /**
   * Operators on other sparse matrices.  These work as though you would expect.
   */
  inline const SpMat&  operator=(const SpMat& m);
  inline const SpMat& operator+=(const SpMat& m);
  inline const SpMat& operator-=(const SpMat& m);
  inline const SpMat& operator*=(const SpMat& m);
  inline const SpMat& operator%=(const SpMat& m);
  inline const SpMat& operator/=(const SpMat& m);
  
  /**
   * Operators on other regular matrices.  These work as though you would expect.
   */
  template<typename T1> inline explicit          SpMat(const Base<eT, T1>& m);
  template<typename T1> inline const SpMat&  operator=(const Base<eT, T1>& m);
  template<typename T1> inline const SpMat& operator*=(const Base<eT, T1>& m);
  template<typename T1> inline const SpMat& operator/=(const Base<eT, T1>& m);
  template<typename T1> inline const SpMat& operator%=(const Base<eT, T1>& m);
  
  
  //! construction of complex matrix out of two non-complex matrices;
  template<typename T1, typename T2>
  inline explicit SpMat(const SpBase<pod_type, T1>& A, const SpBase<pod_type, T2>& B);
  
  /**
   * Operations on sparse subviews.
   */
  inline                   SpMat(const SpSubview<eT>& X);
  inline const SpMat&  operator=(const SpSubview<eT>& X);
  inline const SpMat& operator+=(const SpSubview<eT>& X);
  inline const SpMat& operator-=(const SpSubview<eT>& X);
  inline const SpMat& operator*=(const SpSubview<eT>& X);
  inline const SpMat& operator%=(const SpSubview<eT>& X);
  inline const SpMat& operator/=(const SpSubview<eT>& X);
  
  /**
   * Operations on regular subviews.
   */
  inline                   SpMat(const subview<eT>& x);
  inline const SpMat&  operator=(const subview<eT>& x);
  inline const SpMat& operator+=(const subview<eT>& x);
  inline const SpMat& operator-=(const subview<eT>& x);
  inline const SpMat& operator*=(const subview<eT>& x);
  inline const SpMat& operator%=(const subview<eT>& x);
  inline const SpMat& operator/=(const subview<eT>& x);


  // delayed unary ops
  template<typename T1, typename spop_type> inline                   SpMat(const SpOp<T1, spop_type>& X);
  template<typename T1, typename spop_type> inline const SpMat&  operator=(const SpOp<T1, spop_type>& X);
  template<typename T1, typename spop_type> inline const SpMat& operator+=(const SpOp<T1, spop_type>& X);
  template<typename T1, typename spop_type> inline const SpMat& operator-=(const SpOp<T1, spop_type>& X);
  template<typename T1, typename spop_type> inline const SpMat& operator*=(const SpOp<T1, spop_type>& X);
  template<typename T1, typename spop_type> inline const SpMat& operator%=(const SpOp<T1, spop_type>& X);
  template<typename T1, typename spop_type> inline const SpMat& operator/=(const SpOp<T1, spop_type>& X);
  
  // delayed binary ops
  template<typename T1, typename T2, typename spglue_type> inline                   SpMat(const SpGlue<T1, T2, spglue_type>& X);
  template<typename T1, typename T2, typename spglue_type> inline const SpMat&  operator=(const SpGlue<T1, T2, spglue_type>& X);
  template<typename T1, typename T2, typename spglue_type> inline const SpMat& operator+=(const SpGlue<T1, T2, spglue_type>& X);
  template<typename T1, typename T2, typename spglue_type> inline const SpMat& operator-=(const SpGlue<T1, T2, spglue_type>& X);
  template<typename T1, typename T2, typename spglue_type> inline const SpMat& operator*=(const SpGlue<T1, T2, spglue_type>& X);
  template<typename T1, typename T2, typename spglue_type> inline const SpMat& operator%=(const SpGlue<T1, T2, spglue_type>& X);
  template<typename T1, typename T2, typename spglue_type> inline const SpMat& operator/=(const SpGlue<T1, T2, spglue_type>& X);
  
  /**
   * Submatrix methods.
   */
  arma_inline       SpSubview<eT> row(const uword row_num);
  arma_inline const SpSubview<eT> row(const uword row_num) const;
  
  inline            SpSubview<eT> operator()(const uword row_num, const span& col_span);
  inline      const SpSubview<eT> operator()(const uword row_num, const span& col_span) const;
  
  
  arma_inline       SpSubview<eT> col(const uword col_num);
  arma_inline const SpSubview<eT> col(const uword col_num) const;
  
  inline            SpSubview<eT> operator()(const span& row_span, const uword col_num);
  inline      const SpSubview<eT> operator()(const span& row_span, const uword col_num) const;
  
  /**
   * Row- and column-related functions.
   */
  inline void swap_rows(const uword in_row1, const uword in_row2);
  inline void swap_cols(const uword in_col1, const uword in_col2);
  
  inline void shed_row(const uword row_num);
  inline void shed_col(const uword col_num);
  
  inline void shed_rows(const uword in_row1, const uword in_row2);
  inline void shed_cols(const uword in_col1, const uword in_col2);
  
  arma_inline       SpSubview<eT> rows(const uword in_row1, const uword in_row2);
  arma_inline const SpSubview<eT> rows(const uword in_row1, const uword in_row2) const;
  
  arma_inline       SpSubview<eT> cols(const uword in_col1, const uword in_col2);
  arma_inline const SpSubview<eT> cols(const uword in_col1, const uword in_col2) const;
  
  arma_inline       SpSubview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2);
  arma_inline const SpSubview<eT> submat(const uword in_row1, const uword in_col1, const uword in_row2, const uword in_col2) const;
  
  
  inline            SpSubview<eT> submat    (const span& row_span, const span& col_span);
  inline      const SpSubview<eT> submat    (const span& row_span, const span& col_span) const;
  
  inline            SpSubview<eT> operator()(const span& row_span, const span& col_span);
  inline      const SpSubview<eT> operator()(const span& row_span, const span& col_span) const;
  
  /**
   * Element access; access the i'th element (works identically to the Mat accessors).
   * If there is nothing at element i, 0 is returned.
   *
   * @param i Element to access.
   */
  arma_inline arma_warn_unused SpValProxy<SpMat<eT> > operator[] (const uword i);
  arma_inline arma_warn_unused eT                     operator[] (const uword i) const;
  arma_inline arma_warn_unused SpValProxy<SpMat<eT> > at         (const uword i);
  arma_inline arma_warn_unused eT                     at         (const uword i) const;
  arma_inline arma_warn_unused SpValProxy<SpMat<eT> > operator() (const uword i);
  arma_inline arma_warn_unused eT                     operator() (const uword i) const;
  
  /**
   * Element access; access the element at row in_row and column in_col.
   * If there is nothing at that position, 0 is returned.
   */
  arma_inline arma_warn_unused SpValProxy<SpMat<eT> > at         (const uword in_row, const uword in_col);
  arma_inline arma_warn_unused eT                     at         (const uword in_row, const uword in_col) const;
  arma_inline arma_warn_unused SpValProxy<SpMat<eT> > operator() (const uword in_row, const uword in_col);
  arma_inline arma_warn_unused eT                     operator() (const uword in_row, const uword in_col) const;
  
  
  /**
   * Information boolean checks on matrices.
   */
  arma_inline arma_warn_unused bool is_empty()  const;
  arma_inline arma_warn_unused bool is_vec()    const;
  arma_inline arma_warn_unused bool is_rowvec() const;
  arma_inline arma_warn_unused bool is_colvec() const;
  arma_inline arma_warn_unused bool is_square() const;
       inline arma_warn_unused bool is_finite() const;
  
  arma_inline arma_warn_unused bool in_range(const uword i) const;
  arma_inline arma_warn_unused bool in_range(const span& x) const;
  
  arma_inline arma_warn_unused bool in_range(const uword   in_row, const uword   in_col) const;
  arma_inline arma_warn_unused bool in_range(const span& row_span, const uword   in_col) const;
  arma_inline arma_warn_unused bool in_range(const uword   in_row, const span& col_span) const;
  arma_inline arma_warn_unused bool in_range(const span& row_span, const span& col_span) const;
  
  /**
   * Printing the matrix.
   *
   * @param extra_text Text to prepend to output.
   */
  inline void impl_print(const std::string& extra_text) const;
  inline void impl_print(std::ostream& user_stream, const std::string& extra_text) const;

  inline void impl_raw_print(const std::string& extra_text) const;
  inline void impl_raw_print(std::ostream& user_stream, const std::string& extra_text) const;

  inline void impl_print_dense(const std::string& extra_text) const;
  inline void impl_print_dense(std::ostream& user_stream, const std::string& extra_text) const;
  
  inline void impl_raw_print_dense(const std::string& extra_text) const;
  inline void impl_raw_print_dense(std::ostream& user_stream, const std::string& extra_text) const;
  
  //! Copy the size of another matrix.
  template<typename eT2> inline void copy_size(const SpMat<eT2>& m);
  template<typename eT2> inline void copy_size(const Mat<eT2>& m);

  /**
   * Resize the matrix to a given size.  The matrix will be resized to be a column vector (i.e. in_elem columns, 1 row).
   *
   * @param in_elem Number of elements to allow.
   */
  inline void set_size(const uword in_elem);

  /**
   * Resize the matrix to a given size.
   *
   * @param in_rows Number of rows to allow.
   * @param in_cols Number of columns to allow.
   */
  inline void set_size(const uword in_rows, const uword in_cols);
  
  inline void  reshape(const uword in_rows, const uword in_cols, const uword dim = 0);
  
  inline const SpMat& zeros();
  inline const SpMat& zeros(const uword in_elem);
  inline const SpMat& zeros(const uword in_rows, const uword in_cols);
  
  inline const SpMat& eye();
  inline const SpMat& eye(const uword in_rows, const uword in_cols);
  
  inline const SpMat& speye();
  inline const SpMat& speye(const uword in_rows, const uword in_cols);
  
  inline const SpMat& sprandu(const uword in_rows, const uword in_cols, const double density);
  inline const SpMat& sprandn(const uword in_rows, const uword in_cols, const double density);
  
  inline void reset();
  
  /**
   * Get the minimum or maximum of the matrix.
   */
  inline arma_warn_unused eT min() const;
  inline                  eT min(uword& index_of_min_val) const;
  inline                  eT min(uword& row_of_min_val, uword& col_of_min_val) const;

  inline arma_warn_unused eT max() const;
  inline                  eT max(uword& index_of_max_val) const;
  inline                  eT max(uword& row_of_min_val, uword& col_of_min_val) const;

  // These forward declarations are necessary.
  class iterator_base;
  class iterator;
  class const_iterator;
  class row_iterator;
  class const_row_iterator;

  // Iterator base provides comparison operators but not the actual logic on how
  // to iterate.  The validity of the position is not checked (that is left up
  // to the user).
  class iterator_base
    {
    public:

    inline iterator_base(const SpMat& in_M);
    inline iterator_base(const SpMat& in_M, const uword col, const uword pos);

    inline eT operator*() const;

    inline bool operator==(const iterator_base& rhs) const;
    inline bool operator!=(const iterator_base& rhs) const;

    inline bool operator==(const typename SpSubview<eT>::iterator_base& rhs) const;
    inline bool operator!=(const typename SpSubview<eT>::iterator_base& rhs) const;
    
    // Don't hold location internally; call "dummy" methods to get that information.
    arma_inline uword row() const { return M.row_indices[internal_pos]; }
    arma_inline uword col() const { return internal_col;                }
    arma_inline uword pos() const { return internal_pos;                }

    arma_aligned const SpMat& M;
    arma_aligned       uword  internal_col;
    arma_aligned       uword  internal_pos;
    };

  class const_iterator : public iterator_base
    {
    public:
    
    inline const_iterator(const SpMat& in_M, uword initial_pos = 0); // Assumes initial_pos is valid.
    //! Once initialized, will be at the first nonzero value after the given position (using forward columnwise traversal).
    inline const_iterator(const SpMat& in_M, uword in_row, uword in_col);
    inline const_iterator(const const_iterator& other);

    inline const_iterator& operator++();
    inline void            operator++(int);
    
    inline const_iterator& operator--();
    inline void            operator--(int);
    };

  /**
   * So that we can iterate over nonzero values, we need an iterator
   * implementation.  This can't be as simple as Mat's, which is just a pointer
   * to an eT.  If a value is set to 0 using this iterator, the iterator is no
   * longer valid!
   */
  class iterator : public const_iterator
    {
    public:

    inline iterator(SpMat& in_M, uword initial_pos = 0) : const_iterator(in_M, initial_pos) { }
    inline iterator(SpMat& in_M, uword in_row, uword in_col) : const_iterator(in_M, in_row, in_col) { }
    inline iterator(const const_iterator& other) : const_iterator(other) { }

    inline SpValProxy<SpMat<eT> > operator*();

    // overloads needed for return type correctness
    inline iterator& operator++();
    inline void      operator++(int);

    inline iterator& operator--();
    inline void      operator--(int);
    };

  class const_row_iterator : public iterator_base
    {
    public:
    
    inline const_row_iterator(const SpMat& in_M, uword initial_pos = 0);
    //! Once initialized, will be at the first nonzero value after the given position (using forward row-wise traversal).
    inline const_row_iterator(const SpMat& in_M, uword in_row, uword in_col);
    inline const_row_iterator(const const_row_iterator& other);

    inline const_row_iterator& operator++();
    inline void                operator++(int);
    
    inline const_row_iterator& operator--();
    inline void                operator--(int);

    uword internal_row; // Hold row internally because we use internal_pos differently.
    uword actual_pos; // Actual position in matrix.

    arma_inline eT operator*() const { return iterator_base::M.values[actual_pos]; }

    arma_inline uword row() const { return internal_row; }
    };

  class row_iterator : public const_row_iterator
    {
    public:
    
    inline row_iterator(SpMat& in_M, uword initial_pos = 0) : const_row_iterator(in_M, initial_pos) { }
    //! Once initialized, will be at the first nonzero value after the given position (using forward row-wise traversal).
    inline row_iterator(SpMat& in_M, uword in_row, uword in_col) : const_row_iterator(in_M, in_row, in_col) { }
    inline row_iterator(const row_iterator& other) : const_row_iterator(other) { }
    
    inline SpValProxy<SpMat<eT> > operator*();

    // overloads required for return type correctness
    inline row_iterator& operator++();
    inline void          operator++(int);

    inline row_iterator& operator--();
    inline void          operator--(int);
    };
  
  inline       iterator     begin();
  inline const_iterator     begin() const;
  
  inline       iterator     end();
  inline const_iterator     end() const;
  
  inline       iterator     begin_col(const uword col_num);
  inline const_iterator     begin_col(const uword col_num) const;
  
  inline       iterator     end_col(const uword col_num);
  inline const_iterator     end_col(const uword col_num) const;
  
  inline       row_iterator begin_row(const uword row_num = 0);
  inline const_row_iterator begin_row(const uword row_num = 0) const;
  
  inline       row_iterator end_row();
  inline const_row_iterator end_row() const;
  
  inline       row_iterator end_row(const uword row_num);
  inline const_row_iterator end_row(const uword row_num) const;
  
  inline void  clear();
  inline bool  empty() const;
  inline uword size()  const;
  
  /**
   * Resize memory.  You are responsible for updating the column pointers and
   * filling the new memory (if the new size is larger).  If the new size is
   * smaller, the first new_n_nonzero elements will be copied.  n_nonzero is
   * updated.
   */
  inline void mem_resize(const uword new_n_nonzero);
  
  inline void steal_mem(SpMat& X); //!< don't use this unless you're writing internal Armadillo code
  
  
  protected:

  /**
   * Initialize the matrix to the specified size.  Data is not preserved, so the matrix is assumed to be entirely sparse (empty).
   */
  inline void init(uword in_rows, uword in_cols);

  /**
   * Initialize the matrix from text.  Data is (of course) not preserved, and
   * the size will be reset.
   */
  inline void init(const std::string& text);

  /**
   * Initialize from another matrix (copy).
   */
  inline void init(const SpMat& x);
  
  
  private:
  
  /**
   * Return the given element.
   */
  inline arma_hot arma_warn_unused SpValProxy<SpMat<eT> > get_value(const uword i);
  inline arma_hot arma_warn_unused eT                     get_value(const uword i) const;
  
  inline arma_hot arma_warn_unused SpValProxy<SpMat<eT> > get_value(const uword in_row, const uword in_col);
  inline arma_hot arma_warn_unused eT                     get_value(const uword in_row, const uword in_col) const;
  
  /**
   * Given the index representing which of the nonzero values this is, return
   * its actual location, either in row/col or just the index.
   */
  arma_inline arma_hot arma_warn_unused uword get_position(const uword i) const;
  arma_inline arma_hot                  void  get_position(const uword i, uword& row_of_i, uword& col_of_i) const;
  
  /**
   * Add an element at the given position, and return a reference to it.  The
   * element will be set to 0 (unless otherwise specified).  If the element
   * already exists, its value will be overwritten.
   *
   * @param in_row Row of new element.
   * @param in_col Column of new element.
   * @param in_val Value to set new element to (default 0.0).
   */
  inline arma_hot arma_warn_unused eT& add_element(const uword in_row, const uword in_col, const eT in_val = 0.0);
  
  /**
   * Delete an element at the given position.
   *
   * @param in_row Row of element to be deleted.
   * @param in_col Column of element to be deleted.
   */
  inline arma_hot void delete_element(const uword in_row, const uword in_col);
  
  
  public:
    
  #ifdef ARMA_EXTRA_SPMAT_PROTO
    #include ARMA_INCFILE_WRAP(ARMA_EXTRA_SPMAT_PROTO)
  #endif
  };



#define ARMA_HAS_SPMAT



//! @}
