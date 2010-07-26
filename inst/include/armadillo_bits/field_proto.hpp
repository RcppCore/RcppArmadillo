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
  
  arma_aligned oT** mem;             //!< pointer to memory used by the object
  arma_aligned oT*  mem_local[ 16 ]; //!< Internal memory, to avoid calling the 'new' operator for small amounts of memory
  
  
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
  
  
  inline bool save(const std::string   name, const file_type type = arma_binary, const bool print_status = true) const;
  inline bool save(      std::ostream& os,   const file_type type = arma_binary, const bool print_status = true) const;
  
  inline bool load(const std::string   name, const file_type type = auto_detect, const bool print_status = true);
  inline bool load(      std::istream& is,   const file_type type = auto_detect, const bool print_status = true);
  
  
  inline bool quiet_save(const std::string   name, const file_type type = arma_binary) const;
  inline bool quiet_save(      std::ostream& os,   const file_type type = arma_binary) const;
  
  inline bool quiet_load(const std::string   name, const file_type type = auto_detect);
  inline bool quiet_load(      std::istream& is,   const file_type type = auto_detect);
  
  
  // iterators
  
  class iterator
    {
    public:
    
    inline iterator(field<oT>& in_M, const bool at_end = false);
    
    inline oT& operator* ();
    
    inline iterator& operator++();
    inline void      operator++(int);
    
    inline iterator& operator--();
    inline void      operator--(int);
    
    inline bool operator!=(const iterator& X) const;
    inline bool operator==(const iterator& X) const;
    
    arma_aligned field<oT>& M;
    arma_aligned u32        i;
    };
  
  
  class const_iterator
    {
    public:
    
    const_iterator(const field<oT>& in_M, const bool at_end = false);
    const_iterator(const iterator& X);
    
    inline const oT& operator*() const;
    
    inline const_iterator& operator++();
    inline void            operator++(int);
    
    inline const_iterator& operator--();
    inline void            operator--(int);
    
    inline bool operator!=(const const_iterator& X) const;
    inline bool operator==(const const_iterator& X) const;
    
    arma_aligned const field<oT>& M;
    arma_aligned       u32        i;
    };
  
  inline       iterator begin();
  inline const_iterator begin() const;
  
  inline       iterator end();
  inline const_iterator end()   const;
  
  
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
  
  
  template<typename oT> inline static bool save(const field< oT >& x,       const std::string& name, const file_type type, std::string& err_msg);
  template<typename oT> inline static bool save(const field< oT >& x,       std::ostream& os,        const file_type type, std::string& err_msg);
  template<typename oT> inline static bool load(      field< oT >& x,       const std::string& name, const file_type type, std::string& err_msg);
  template<typename oT> inline static bool load(      field< oT >& x,       std::istream& is,        const file_type type, std::string& err_msg);

  template<typename eT> inline static bool save(const field< Mat<eT> >& x,  const std::string& name, const file_type type, std::string& err_msg);
  template<typename eT> inline static bool save(const field< Mat<eT> >& x,  std::ostream& os,        const file_type type, std::string& err_msg);
  template<typename eT> inline static bool load(      field< Mat<eT> >& x,  const std::string& name, const file_type type, std::string& err_msg);
  template<typename eT> inline static bool load(      field< Mat<eT> >& x,  std::istream& is,        const file_type type, std::string& err_msg);
  
  template<typename eT> inline static bool save(const field< Col<eT> >& x,  const std::string& name, const file_type type, std::string& err_msg);
  template<typename eT> inline static bool save(const field< Col<eT> >& x,  std::ostream& os,        const file_type type, std::string& err_msg);
  template<typename eT> inline static bool load(      field< Col<eT> >& x,  const std::string& name, const file_type type, std::string& err_msg);
  template<typename eT> inline static bool load(      field< Col<eT> >& x,  std::istream& is,        const file_type type, std::string& err_msg);
  
  template<typename eT> inline static bool save(const field< Row<eT> >& x,  const std::string& name, const file_type type, std::string& err_msg);
  template<typename eT> inline static bool save(const field< Row<eT> >& x,  std::ostream& os,        const file_type type, std::string& err_msg);
  template<typename eT> inline static bool load(      field< Row<eT> >& x,  const std::string& name, const file_type type, std::string& err_msg);
  template<typename eT> inline static bool load(      field< Row<eT> >& x,  std::istream& is,        const file_type type, std::string& err_msg);

  template<typename eT> inline static bool save(const field< Cube<eT> >& x, const std::string& name, const file_type type, std::string& err_msg);
  template<typename eT> inline static bool save(const field< Cube<eT> >& x, std::ostream& os,        const file_type type, std::string& err_msg);
  template<typename eT> inline static bool load(      field< Cube<eT> >& x, const std::string& name, const file_type type, std::string& err_msg);
  template<typename eT> inline static bool load(      field< Cube<eT> >& x, std::istream& is,        const file_type type, std::string& err_msg);
  
  inline static bool save(const field< std::string >& x, const std::string& name, const file_type type, std::string& err_msg);
  inline static bool save(const field< std::string >& x, std::ostream& os,        const file_type type, std::string& err_msg);
  inline static bool load(      field< std::string >& x, const std::string& name, const file_type type, std::string& err_msg);
  inline static bool load(      field< std::string >& x, std::istream& is,        const file_type type, std::string& err_msg);
  
  };


//! @}
