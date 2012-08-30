// Copyright (C) 2011-2012 Ryan Curtin <ryan@igglybob.com>
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


//! \addtogroup SpRow
//! @{

//! Class for sparse row vectors (sparse matrices with only one row)

template<typename eT>
class SpRow : public SpMat<eT>
  {
  public:

  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  
  static const bool is_row = true;
  static const bool is_col = false;
  
  
  inline          SpRow();
  inline explicit SpRow(const uword N);
  inline          SpRow(const uword in_rows, const uword in_cols);

  inline                  SpRow(const char*        text);
  inline const SpRow& operator=(const char*        text);

  inline                  SpRow(const std::string& text);
  inline const SpRow& operator=(const std::string& text);

  inline const SpRow& operator=(const eT val);

  template<typename T1> inline                   SpRow(const Base<eT,T1>& X);
  template<typename T1> inline const SpRow&  operator=(const Base<eT,T1>& X);

  template<typename T1, typename T2>
  inline explicit SpRow(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);

  arma_inline SpValProxy<SpMat<eT> > col(const uword col_num);
  arma_inline eT                     col(const uword col_num) const;

//  arma_inline       subview_row<eT> cols(const uword in_col1, const uword in_col2);
//  arma_inline const subview_row<eT> cols(const uword in_col1, const uword in_col2) const;

//  arma_inline       subview_row<eT> subvec(const uword in_col1, const uword in_col2);
//  arma_inline const subview_row<eT> subvec(const uword in_col1, const uword in_col2) const;

//  arma_inline       subview_row<eT> subvec(const span& col_span);
//  arma_inline const subview_row<eT> subvec(const span& col_span) const;

//  arma_inline       subview_row<eT> operator()(const span& col_span);
//  arma_inline const subview_row<eT> operator()(const span& col_span) const;

  inline void shed_col (const uword col_num);
  inline void shed_cols(const uword in_col1, const uword in_col2);

                        inline void insert_cols(const uword col_num, const uword N, const bool set_to_zero = true);
  template<typename T1> inline void insert_cols(const uword col_num, const Base<eT,T1>& X);


  typedef typename SpMat<eT>::iterator       row_iterator;
  typedef typename SpMat<eT>::const_iterator const_row_iterator;

  inline       row_iterator begin_row();
  inline const_row_iterator begin_row() const;

  inline       row_iterator end_row();
  inline const_row_iterator end_row() const;
  
  #ifdef ARMA_EXTRA_ROW_PROTO
    #include ARMA_INCFILE_WRAP(ARMA_EXTRA_ROW_PROTO)
  #endif
  };



//! @}
