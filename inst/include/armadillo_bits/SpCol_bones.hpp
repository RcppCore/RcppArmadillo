// Copyright (C) 2011-2012 Ryan Curtin
// Copyright (C) 2011 Matthew Amidon
//
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup SpCol
//! @{

//! Class for sparse column vectors (matrices with only one column)

template<typename eT>
class SpCol : public SpMat<eT>
  {
  public:

  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  
  static const bool is_row = false;
  static const bool is_col = true;
  
  
  inline          SpCol();
  inline explicit SpCol(const uword n_elem);
  inline          SpCol(const uword in_rows, const uword in_cols);

  inline                  SpCol(const char*        text);
  inline const SpCol& operator=(const char*        text);

  inline                  SpCol(const std::string& text);
  inline const SpCol& operator=(const std::string& text);

  inline const SpCol& operator=(const eT val);

  template<typename T1> inline                  SpCol(const Base<eT,T1>& X);
  template<typename T1> inline const SpCol& operator=(const Base<eT,T1>& X);

  template<typename T1> inline                  SpCol(const SpBase<eT,T1>& X);
  template<typename T1> inline const SpCol& operator=(const SpBase<eT,T1>& X);

  template<typename T1, typename T2>
  inline explicit SpCol(const SpBase<pod_type,T1>& A, const SpBase<pod_type,T2>& B);

  inline SpValProxy<SpMat<eT> > row(const uword row_num);
  inline eT                     row(const uword row_num) const;

//  arma_inline       subview_col<eT> rows(const uword in_row1, const uword in_row2);
//  arma_inline const subview_col<eT> rows(const uword in_row1, const uword in_row2) const;

//  arma_inline       subview_col<eT> subvec(const uword in_row1, const uword in_row2);
//  arma_inline const subview_col<eT> subvec(const uword in_row1, const uword in_row2) const;

//  arma_inline       subview_col<eT> subvec(const span& row_span);
//  arma_inline const subview_col<eT> subvec(const span& row_span) const;

  inline void shed_row (const uword row_num);
  inline void shed_rows(const uword in_row1, const uword in_row2);

//                         inline void insert_rows(const uword row_num, const uword N, const bool set_to_zero = true);
//   template<typename T1> inline void insert_rows(const uword row_num, const Base<eT,T1>& X);


  typedef typename SpMat<eT>::iterator       row_iterator;
  typedef typename SpMat<eT>::const_iterator const_row_iterator;

  inline       row_iterator begin_row(const uword row_num = 0);
  inline const_row_iterator begin_row(const uword row_num = 0) const;

  inline       row_iterator end_row  (const uword row_num = 0);
  inline const_row_iterator end_row  (const uword row_num = 0) const;
  
  
  #ifdef ARMA_EXTRA_SPCOL_PROTO
    #include ARMA_INCFILE_WRAP(ARMA_EXTRA_SPCOL_PROTO)
  #endif
  };
