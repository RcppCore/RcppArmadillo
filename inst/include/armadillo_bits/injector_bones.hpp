// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup injector
//! @{



template<typename eT>
class mat_injector_row
  {
  public:
  
  inline      mat_injector_row();
  
  inline void insert(const eT val) const;
  
  mutable uword        n_cols;
  mutable podarray<eT> A;
  mutable podarray<eT> B;
  };



template<typename T1>
class mat_injector
  {
  public:
  
  typedef typename T1::elem_type elem_type;
  
  inline void  insert(const elem_type val) const;
  inline void  end_of_row()                const;
  inline      ~mat_injector();
  
  
  private:
  
  inline mat_injector(T1& in_X, const elem_type val);
  inline mat_injector(T1& in_X, const injector_end_of_row<>& x);
  
  T1&           X;
  mutable uword n_rows;
  
  mutable podarray< mat_injector_row<elem_type>* >* AA;
  mutable podarray< mat_injector_row<elem_type>* >* BB;
  
  friend class Mat<elem_type>;
  friend class Row<elem_type>;
  friend class Col<elem_type>;
  };



//



template<typename oT>
class field_injector_row
  {
  public:
  
  inline      field_injector_row();
  inline     ~field_injector_row();
  
  inline void insert(const oT& val) const;
  
  mutable uword      n_cols;
  mutable field<oT>* AA;
  mutable field<oT>* BB;
  };

  
  
template<typename T1>
class field_injector
  {
  public:
  
  typedef typename T1::object_type object_type;
  
  inline void  insert(const object_type& val) const;
  inline void  end_of_row()                   const;
  inline      ~field_injector();
  
  
  private:
  
  inline field_injector(T1& in_X, const object_type& val);
  inline field_injector(T1& in_X, const injector_end_of_row<>& x);
  
  T1&           X;
  mutable uword n_rows;
  
  mutable podarray< field_injector_row<object_type>* >* AA;
  mutable podarray< field_injector_row<object_type>* >* BB;
  
  friend class field<object_type>;
  };



//! @}
