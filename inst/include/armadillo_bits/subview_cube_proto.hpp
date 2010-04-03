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


//! \addtogroup subview_cube
//! @{


//! Class for storing data required to construct or apply operations to a subcube
//! (i.e. where the subcube starts and ends as well as a reference/pointer to the original cube),
template<typename eT>
class subview_cube : public BaseCube<eT, subview_cube<eT> >
  {
  public:    arma_aligned const Cube<eT>& m;
  protected: arma_aligned       Cube<eT>* m_ptr;
  
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;

  const u32 aux_row1;
  const u32 aux_col1;
  const u32 aux_slice1;
  
  const u32 aux_row2;
  const u32 aux_col2;
  const u32 aux_slice2;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem_slice;
  const u32 n_slices;
  const u32 n_elem;
  
  
  protected:
  
  arma_inline subview_cube(const Cube<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_slice1, const u32 in_row2, const u32 in_col2, const u32 in_slice2);
  arma_inline subview_cube(      Cube<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_slice1, const u32 in_row2, const u32 in_col2, const u32 in_slice2);
  
  
  public:
  
  inline ~subview_cube();
  
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  // deliberately returning void
  template<typename T1> inline void operator=  (const BaseCube<eT,T1>& x);
  template<typename T1> inline void operator+= (const BaseCube<eT,T1>& x);
  template<typename T1> inline void operator-= (const BaseCube<eT,T1>& x);
  template<typename T1> inline void operator%= (const BaseCube<eT,T1>& x);
  template<typename T1> inline void operator/= (const BaseCube<eT,T1>& x);
  
  inline void operator=  (const subview_cube& x);
  inline void operator+= (const subview_cube& x);
  inline void operator-= (const subview_cube& x);
  inline void operator%= (const subview_cube& x);
  inline void operator/= (const subview_cube& x);
  
  template<typename T1> inline void operator=  (const Base<eT,T1>& x);
  template<typename T1> inline void operator+= (const Base<eT,T1>& x);
  template<typename T1> inline void operator-= (const Base<eT,T1>& x);
  template<typename T1> inline void operator%= (const Base<eT,T1>& x);
  template<typename T1> inline void operator/= (const Base<eT,T1>& x);

  inline static void extract(Cube<eT>& out, const subview_cube& in);
  inline static void extract(Mat<eT>&  out, const subview_cube& in);
  
  inline static void  plus_inplace(Cube<eT>& out, const subview_cube& in);
  inline static void minus_inplace(Cube<eT>& out, const subview_cube& in);
  inline static void schur_inplace(Cube<eT>& out, const subview_cube& in);
  inline static void   div_inplace(Cube<eT>& out, const subview_cube& in);
  
  inline static void  plus_inplace(Mat<eT>& out, const subview_cube& in);
  inline static void minus_inplace(Mat<eT>& out, const subview_cube& in);
  inline static void schur_inplace(Mat<eT>& out, const subview_cube& in);
  inline static void   div_inplace(Mat<eT>& out, const subview_cube& in);

  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  
  arma_inline eT& operator[](const u32 i);
  arma_inline eT  operator[](const u32 i) const;
  
  arma_inline eT& operator()(const u32 i);
  arma_inline eT  operator()(const u32 i) const;
  
  arma_inline eT& operator()(const u32 in_row, const u32 in_col, const u32 in_slice);
  arma_inline eT  operator()(const u32 in_row, const u32 in_col, const u32 in_slice) const;
  
  arma_inline eT&         at(const u32 in_row, const u32 in_col, const u32 in_slice);
  arma_inline eT          at(const u32 in_row, const u32 in_col, const u32 in_slice) const;
  
  arma_inline       eT* slice_colptr(const u32 in_slice, const u32 in_col);
  arma_inline const eT* slice_colptr(const u32 in_slice, const u32 in_col) const;
  
  inline bool check_overlap(const subview_cube& x) const;
  inline bool check_overlap(const Mat<eT>&      x) const;
  
  
  private:
  
  friend class  Mat<eT>;
  friend class Cube<eT>;
  
  subview_cube();
  };



//! @}
