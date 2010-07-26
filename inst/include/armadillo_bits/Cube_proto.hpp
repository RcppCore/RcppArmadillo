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


//! \addtogroup Cube
//! @{

//! Dense cube class

template<typename eT>
class Cube : public BaseCube< eT, Cube<eT> >
  {
  public:
  
  typedef eT                                elem_type; //!< the type of elements stored in the cube
  typedef typename get_pod_type<eT>::result pod_type;  //!< if eT is non-complex, pod_type is same as eT. otherwise, pod_type is the underlying type used by std::complex
  
  const u32  n_rows;       //!< number of rows in each slice (read-only)
  const u32  n_cols;       //!< number of columns in each slice (read-only)
  const u32  n_elem_slice; //!< number of elements in each slice (read-only)
  const u32  n_slices;     //!< number of slices in the cube (read-only)
  const u32  n_elem;       //!< number of elements in the cube (read-only)
  const bool use_aux_mem;  //!< true if externally managed memory is being used (read-only)
  
  arma_aligned const Mat<eT>** const mat_ptrs; //!< pointer to an array containing pointers to Mat instances (one for each slice)
  arma_aligned const eT*       const mem;      //!< pointer to the memory used by the cube (memory is read-only)
  
  protected:
  arma_aligned Mat<eT>* mat_ptrs_local[ 16 ];
  arma_aligned eT            mem_local[ 64 ];
  
  
  public:
  
  inline ~Cube();
  inline  Cube();
  
  inline Cube(const u32 in_rows, const u32 in_cols, const u32 in_slices);
  
  inline Cube(      eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const u32 aux_n_slices, const bool copy_aux_mem = true);
  inline Cube(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols, const u32 aux_n_slices);
  
  arma_inline const Cube&  operator=(const eT val);
  arma_inline const Cube& operator+=(const eT val);
  arma_inline const Cube& operator-=(const eT val);
  arma_inline const Cube& operator*=(const eT val);
  arma_inline const Cube& operator/=(const eT val);
  
  inline                   Cube(const Cube& m);
  inline const Cube&  operator=(const Cube& m);
  inline const Cube& operator+=(const Cube& m);
  inline const Cube& operator-=(const Cube& m);
  inline const Cube& operator%=(const Cube& m);
  inline const Cube& operator/=(const Cube& m);

  template<typename T1, typename T2>
  inline explicit Cube(const BaseCube<pod_type,T1>& A, const BaseCube<pod_type,T2>& B);

  inline                   Cube(const subview_cube<eT>& X);
  inline const Cube&  operator=(const subview_cube<eT>& X);
  inline const Cube& operator+=(const subview_cube<eT>& X);
  inline const Cube& operator-=(const subview_cube<eT>& X);
  inline const Cube& operator%=(const subview_cube<eT>& X);
  inline const Cube& operator/=(const subview_cube<eT>& X);

  arma_inline       Mat<eT>& slice(const u32 in_slice);
  arma_inline const Mat<eT>& slice(const u32 in_slice) const;
  
  arma_inline       subview_cube<eT> slices(const u32 in_slice1, const u32 in_slice2);
  arma_inline const subview_cube<eT> slices(const u32 in_slice1, const u32 in_slice2) const;
  
  arma_inline       subview_cube<eT> subcube(const u32 in_row1, const u32 in_col1, const u32 in_slice1, const u32 in_row2, const u32 in_col2, const u32 in_slice2);
  arma_inline const subview_cube<eT> subcube(const u32 in_row1, const u32 in_col1, const u32 in_slice1, const u32 in_row2, const u32 in_col2, const u32 in_slice2) const;


  template<typename T1, typename op_type> inline                   Cube(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Cube&  operator=(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Cube& operator+=(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Cube& operator-=(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Cube& operator%=(const OpCube<T1, op_type>& X);
  template<typename T1, typename op_type> inline const Cube& operator/=(const OpCube<T1, op_type>& X);
  
  template<typename T1, typename eop_type> inline                   Cube(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Cube&  operator=(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Cube& operator+=(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Cube& operator-=(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Cube& operator%=(const eOpCube<T1, eop_type>& X);
  template<typename T1, typename eop_type> inline const Cube& operator/=(const eOpCube<T1, eop_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline                   Cube(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Cube&  operator=(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Cube& operator+=(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Cube& operator-=(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Cube& operator%=(const GlueCube<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const Cube& operator/=(const GlueCube<T1, T2, glue_type>& X);
  
  template<typename T1, typename T2, typename eglue_type> inline                   Cube(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Cube&  operator=(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Cube& operator+=(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Cube& operator-=(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Cube& operator%=(const eGlueCube<T1, T2, eglue_type>& X);
  template<typename T1, typename T2, typename eglue_type> inline const Cube& operator/=(const eGlueCube<T1, T2, eglue_type>& X);
  
  arma_inline eT& operator[] (const u32 i);
  arma_inline eT  operator[] (const u32 i) const;
  arma_inline eT& operator() (const u32 i);
  arma_inline eT  operator() (const u32 i) const;
  
  arma_inline eT& at         (const u32 in_row, const u32 in_col, const u32 in_slice);
  arma_inline eT  at         (const u32 in_row, const u32 in_col, const u32 in_slice) const;
  arma_inline eT& operator() (const u32 in_row, const u32 in_col, const u32 in_slice);
  arma_inline eT  operator() (const u32 in_row, const u32 in_col, const u32 in_slice) const;

  arma_inline const Cube& operator++();
  arma_inline void        operator++(int);

  arma_inline const Cube& operator--();
  arma_inline void        operator--(int);

  arma_inline bool is_finite() const;

  arma_inline       eT* memptr();
  arma_inline const eT* memptr() const;

  arma_inline       eT* slice_memptr(const u32 slice);
  arma_inline const eT* slice_memptr(const u32 slice) const;

  arma_inline       eT* slice_colptr(const u32 in_slice, const u32 in_col);
  arma_inline const eT* slice_colptr(const u32 in_slice, const u32 in_col) const;
  
  inline void print(const std::string extra_text = "") const;
  inline void print(std::ostream& user_stream, const std::string extra_text = "") const;

  inline void raw_print(const std::string extra_text = "") const;
  inline void raw_print(std::ostream& user_stream, const std::string extra_text = "") const;

  inline void  set_size(const u32 in_rows, const u32 in_cols, const u32 in_slices);
  
  template<typename eT2> inline void copy_size(const Cube<eT2>& m);

  inline void fill(const eT val);
  
  inline void zeros();
  inline void zeros(const u32 in_rows, const u32 in_cols, const u32 in_slices);
  
  inline void ones();
  inline void ones(const u32 in_rows, const u32 in_cols, const u32 in_slices);
  
  inline void reset();
  
  
  inline bool save(const std::string   name, const file_type type = arma_binary, const bool print_status = true) const;
  inline bool save(      std::ostream& os,   const file_type type = arma_binary, const bool print_status = true) const;
  
  inline bool load(const std::string   name, const file_type type = auto_detect, const bool print_status = true);
  inline bool load(      std::istream& is,   const file_type type = auto_detect, const bool print_status = true);
  
  inline bool quiet_save(const std::string   name, const file_type type = arma_binary) const;
  inline bool quiet_save(      std::ostream& os,   const file_type type = arma_binary) const;
  
  inline bool quiet_load(const std::string   name, const file_type type = auto_detect);
  inline bool quiet_load(      std::istream& is,   const file_type type = auto_detect);
  
  
  // iterators
  
  typedef       eT*       iterator;
  typedef const eT* const_iterator;
  
  typedef       eT*       slice_iterator;
  typedef const eT* const_slice_iterator;
  
  inline       iterator begin();
  inline const_iterator begin() const;
  
  inline       iterator end();
  inline const_iterator end()   const;
  
  inline       slice_iterator begin_slice(const u32 slice_num);
  inline const_slice_iterator begin_slice(const u32 slice_num) const;
  
  inline       slice_iterator end_slice(const u32 slice_num);
  inline const_slice_iterator end_slice(const u32 slice_num)   const;
  
  
  protected:
  
  inline void init(const u32 in_rows, const u32 in_cols, const u32 in_slices);
  inline void init(const Cube& x);
  
  inline void delete_mat();
  inline void create_mat();
  };



class Cube_aux
  {
  public:

  template<typename eT> arma_inline static void prefix_pp(Cube<eT>& x);
  template<typename T>  arma_inline static void prefix_pp(Cube< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void postfix_pp(Cube<eT>& x);
  template<typename T>  arma_inline static void postfix_pp(Cube< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void prefix_mm(Cube<eT>& x);
  template<typename T>  arma_inline static void prefix_mm(Cube< std::complex<T> >& x);
  
  template<typename eT> arma_inline static void postfix_mm(Cube<eT>& x);
  template<typename T>  arma_inline static void postfix_mm(Cube< std::complex<T> >& x);
  };



//! @}
