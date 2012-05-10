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


//! \addtogroup diagview
//! @{


//! Class for storing data required to extract and set the diagonals of a matrix
template<typename eT>
class diagview : public Base<eT, diagview<eT> >
  {
  public:
  
  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  
  arma_aligned const Mat<eT>& m;
  
  static const bool is_row = false;
  static const bool is_col = true;
  
  const uword row_offset;
  const uword col_offset;
  
  const uword n_rows;     // equal to n_elem
  const uword n_elem;
  
  static const uword n_cols = 1;
  
  
  protected:
  
  arma_inline diagview(const Mat<eT>& in_m, const uword in_row_offset, const uword in_col_offset, const uword len);
  
  
  public:
  
  inline ~diagview();
  
  inline void operator=(const diagview& x);
  
  inline void operator+=(const eT val);
  inline void operator-=(const eT val);
  inline void operator*=(const eT val);
  inline void operator/=(const eT val);
  
  template<typename T1> inline void operator= (const Base<eT,T1>& x);
  template<typename T1> inline void operator+=(const Base<eT,T1>& x);
  template<typename T1> inline void operator-=(const Base<eT,T1>& x);
  template<typename T1> inline void operator%=(const Base<eT,T1>& x);
  template<typename T1> inline void operator/=(const Base<eT,T1>& x);
  
  
  arma_inline eT& operator[](const uword ii);
  arma_inline eT  operator[](const uword ii) const;
  
  arma_inline eT&         at(const uword ii);
  arma_inline eT          at(const uword ii) const;
  
  arma_inline eT& operator()(const uword ii);
  arma_inline eT  operator()(const uword ii) const;
  
  arma_inline eT&         at(const uword in_n_row, const uword);
  arma_inline eT          at(const uword in_n_row, const uword) const;
   
  arma_inline eT& operator()(const uword in_n_row, const uword in_n_col);
  arma_inline eT  operator()(const uword in_n_row, const uword in_n_col) const;
  
  
  arma_inline const Op<diagview<eT>,op_htrans>  t() const;
  arma_inline const Op<diagview<eT>,op_htrans> ht() const;
  arma_inline const Op<diagview<eT>,op_strans> st() const;
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  
  inline static void extract(Mat<eT>& out, const diagview& in);
  
  inline static void  plus_inplace(Mat<eT>& out, const diagview& in);
  inline static void minus_inplace(Mat<eT>& out, const diagview& in);
  inline static void schur_inplace(Mat<eT>& out, const diagview& in);
  inline static void   div_inplace(Mat<eT>& out, const diagview& in);
  
  
  private:
  
  friend class Mat<eT>;
  friend class subview<eT>;
  
  diagview();
  //diagview(const diagview&);  // making this private causes an error under gcc 4.1/4.2, but not 4.3
  };


//! @}
