// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup subview_field
//! @{


//! Class for storing data required to construct or apply operations to a subfield
//! (i.e. where the subfield starts and ends as well as a reference/pointer to the original field),
template<typename oT>
class subview_field
  {
  public:  const field<oT>& f;
  protected:     field<oT>* f_ptr;
  
  public:
  
  typedef oT object_type;
  
  const u32 aux_row1;
  const u32 aux_col1;
  
  const u32 aux_row2;
  const u32 aux_col2;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  
  protected:
  
  arma_inline subview_field(const field<oT>& in_f, const u32 in_row1, const u32 in_col1, const u32 in_row2,  const u32 in_col2);
  arma_inline subview_field(      field<oT>& in_f, const u32 in_row1, const u32 in_col1, const u32 in_row2,  const u32 in_col2);
  
  
  public:
  
  inline ~subview_field();
  
  inline void operator= (const field<oT>& x);
  inline void operator= (const subview_field& x);
  
  arma_inline       oT& operator[](const u32 i);
  arma_inline const oT& operator[](const u32 i) const;
  
  arma_inline       oT& operator()(const u32 i);
  arma_inline const oT& operator()(const u32 i) const;
  
  arma_inline       oT&         at(const u32 row, const u32 col);
  arma_inline const oT&         at(const u32 row, const u32 col) const;
  
  arma_inline       oT& operator()(const u32 row, const u32 col);
  arma_inline const oT& operator()(const u32 row, const u32 col) const;
  
  inline bool check_overlap(const subview_field& x) const;
  
  inline static void extract(field<oT>& out, const subview_field& in);
  
  
  private:
  
  friend class field<oT>;
  
  
  subview_field();
  //subview_field(const subview_field&);
  };


//! @}
