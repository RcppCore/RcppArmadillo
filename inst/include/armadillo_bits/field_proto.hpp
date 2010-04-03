// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Ian Cullinan (ian dot cullinan at nicta dot com dot au)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup field
//! @{



//! A lightweight 2D container for abitrary objects
//! (the objects must have a copy constructor)

template<typename oT>
class field
  {
  public:
  
  typedef oT object_type;
  
  const u32 n_rows;     //!< number of rows in the field (read-only)
  const u32 n_cols;     //!< number of columns in the field (read-only)
  const u32 n_elem;     //!< number of elements in the field (read-only)

  private:
  
  //! pointer to memory used by the object
  arma_aligned oT** mem;
  arma_aligned oT*  mem_local[ 16 ];
  //!< Internal memory, to avoid calling the 'new' operator for small amounts of memory
  
  
  public:
  
  inline ~field();
  inline  field();
  
  inline                  field(const field& x);
  inline const field& operator=(const field& x);
  
  inline                  field(const subview_field<oT>& x);
  inline const field& operator=(const subview_field<oT>& x);
  
  inline explicit field(const u32 n_elem_in);
  inline          field(const u32 n_rows_in, const u32 n_cols_in);
  
  inline void  set_size(const u32 n_obj_in);
  inline void  set_size(const u32 n_rows_in, const u32 n_cols_in);
  
  template<typename oT2>
  inline void copy_size(const field<oT2>& x);
  
  arma_inline       oT& operator[](const u32 i);
  arma_inline const oT& operator[](const u32 i) const;
  
  arma_inline       oT& operator()(const u32 i);
  arma_inline const oT& operator()(const u32 i) const;
  
  arma_inline       oT&         at(const u32 row, const u32 col);
  arma_inline const oT&         at(const u32 row, const u32 col) const;
  
  arma_inline       oT& operator()(const u32 row, const u32 col);
  arma_inline const oT& operator()(const u32 row, const u32 col) const;
  
  inline       subview_field<oT> row(const u32 row_num);
  inline const subview_field<oT> row(const u32 row_num) const;
  
  inline       subview_field<oT> col(const u32 col_num);
  inline const subview_field<oT> col(const u32 col_num) const;
  
  inline       subview_field<oT> rows(const u32 in_row1, const u32 in_row2);
  inline const subview_field<oT> rows(const u32 in_row1, const u32 in_row2) const;
  
  inline       subview_field<oT> cols(const u32 in_col1, const u32 in_col2);
  inline const subview_field<oT> cols(const u32 in_col1, const u32 in_col2) const;
  
  inline       subview_field<oT> subfield(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2);
  inline const subview_field<oT> subfield(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2) const;
  
  inline void print(const std::string extra_text = "") const;
  inline void print(std::ostream& user_stream, const std::string extra_text = "") const;
  
  inline void fill(const oT& x);
  
  inline void reset();
  inline void reset_objects();
  
  inline void save(const std::string   name, const file_type type = arma_binary) const;
  inline void save(      std::ostream& os,   const file_type type = arma_binary) const;
  
  inline void load(const std::string   name, const file_type type = auto_detect);
  inline void load(      std::istream& is,   const file_type type = auto_detect);
  
  
  private:
  
  inline void init(const field<oT>& x);
  inline void init(const u32 n_rows_in, const u32 n_cols_in);
  
  inline void delete_objects();
  inline void create_objects();
  
  friend class field_aux;
  friend class subview_field<oT>;
  };



class field_aux
  {
  public:
  
  template<typename oT> inline static void reset_objects(field< oT >& x);
  template<typename eT> inline static void reset_objects(field< Mat<eT> >& x);
  template<typename eT> inline static void reset_objects(field< Col<eT> >& x);
  template<typename eT> inline static void reset_objects(field< Row<eT> >& x);
  template<typename eT> inline static void reset_objects(field< Cube<eT> >& x);
                        inline static void reset_objects(field< std::string >& x);
  
  
  template<typename oT> inline static void save(const field< oT >& x,          const std::string& name, const file_type type);
  template<typename oT> inline static void save(const field< oT >& x,          std::ostream& os,        const file_type type);
  template<typename oT> inline static void load(      field< oT >& x,          const std::string& name, const file_type type);
  template<typename oT> inline static void load(      field< oT >& x,          std::istream& is,        const file_type type);

  template<typename eT> inline static void save(const field< Mat<eT> >& x,     const std::string& name, const file_type type);
  template<typename eT> inline static void save(const field< Mat<eT> >& x,     std::ostream& os,        const file_type type);
  template<typename eT> inline static void load(      field< Mat<eT> >& x,     const std::string& name, const file_type type);
  template<typename eT> inline static void load(      field< Mat<eT> >& x,     std::istream& is,        const file_type type);
  
  template<typename eT> inline static void save(const field< Col<eT> >& x,     const std::string& name, const file_type type);
  template<typename eT> inline static void save(const field< Col<eT> >& x,     std::ostream& os,        const file_type type);
  template<typename eT> inline static void load(      field< Col<eT> >& x,     const std::string& name, const file_type type);
  template<typename eT> inline static void load(      field< Col<eT> >& x,     std::istream& is,        const file_type type);
  
  template<typename eT> inline static void save(const field< Row<eT> >& x,     const std::string& name, const file_type type);
  template<typename eT> inline static void save(const field< Row<eT> >& x,     std::ostream& os,        const file_type type);
  template<typename eT> inline static void load(      field< Row<eT> >& x,     const std::string& name, const file_type type);
  template<typename eT> inline static void load(      field< Row<eT> >& x,     std::istream& is,        const file_type type);

  template<typename eT> inline static void save(const field< Cube<eT> >& x,    const std::string& name, const file_type type);
  template<typename eT> inline static void save(const field< Cube<eT> >& x,    std::ostream& os,        const file_type type);
  template<typename eT> inline static void load(      field< Cube<eT> >& x,    const std::string& name, const file_type type);
  template<typename eT> inline static void load(      field< Cube<eT> >& x,    std::istream& is,        const file_type type);
  
                        inline static void save(const field< std::string >& x, const std::string& name, const file_type type);
                        inline static void save(const field< std::string >& x, std::ostream& os,        const file_type type);
                        inline static void load(      field< std::string >& x, const std::string& name, const file_type type);
                        inline static void load(      field< std::string >& x, std::istream& is,        const file_type type);
  
  };


//! @}
