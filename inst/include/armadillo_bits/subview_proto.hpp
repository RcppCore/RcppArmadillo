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
  
  const u32 aux_row2;
  const u32 aux_col2;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  
  protected:
  
  arma_inline subview(const Mat<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_row2,  const u32 in_col2);
  arma_inline subview(      Mat<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_row2,  const u32 in_col2);
  
  
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
  
  arma_inline eT& operator[](const u32 i);
  arma_inline eT  operator[](const u32 i) const;
  
  arma_inline eT& operator()(const u32 i);
  arma_inline eT  operator()(const u32 i) const;
  
  arma_inline eT& operator()(const u32 in_row, const u32 in_col);
  arma_inline eT  operator()(const u32 in_row, const u32 in_col) const;
  
  arma_inline eT&         at(const u32 in_row, const u32 in_col);
  arma_inline eT          at(const u32 in_row, const u32 in_col) const;
  
  arma_inline       eT* colptr(const u32 in_col);
  arma_inline const eT* colptr(const u32 in_col) const;
  
  inline bool check_overlap(const subview& x) const;
  
  inline bool is_vec() const;
  
  
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
  
  
  protected:
  
  arma_inline subview_col(const Mat<eT>& in_m, const u32 in_col);
  arma_inline subview_col(      Mat<eT>& in_m, const u32 in_col);
  
  arma_inline subview_col(const Mat<eT>& in_m, const u32 in_col, const u32 in_row1, const u32 in_row2);
  arma_inline subview_col(      Mat<eT>& in_m, const u32 in_col, const u32 in_row1, const u32 in_row2);
  
  
  private:
  
  friend class Mat<eT>;
  friend class Col<eT>;
  
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
  
  
  protected:
  
  arma_inline subview_row(const Mat<eT>& in_m, const u32 in_row);
  arma_inline subview_row(      Mat<eT>& in_m, const u32 in_row);
  
  arma_inline subview_row(const Mat<eT>& in_m, const u32 in_row, const u32 in_col1, const u32 in_col2);
  arma_inline subview_row(      Mat<eT>& in_m, const u32 in_row, const u32 in_col1, const u32 in_col2);
  
  
  private:
  
  friend class Mat<eT>;
  friend class Row<eT>;
  
  subview_row();
  };



//! @}
