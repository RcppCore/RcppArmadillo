// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup Row
//! @{

//! Class for row vectors (matrices with only one row)

template<typename eT>
class Row : public Mat<eT>
  {
  public:
  
  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  
  
  inline          Row();
  inline explicit Row(const u32 N);
  inline          Row(const u32 in_rows, const u32 in_cols);
  
  inline                  Row(const char*        text);
  inline const Row& operator=(const char*        text);
  inline                  Row(const std::string& text);
  inline const Row& operator=(const std::string& text);
  
  inline const Row& operator=(const eT val);

  template<typename T1> inline                   Row(const Base<eT,T1>& X);
  template<typename T1> inline const Row&  operator=(const Base<eT,T1>& X);
  
  inline Row(      eT* aux_mem, const u32 aux_length, const bool copy_aux_mem = true, const bool strict = true);
  inline Row(const eT* aux_mem, const u32 aux_length);
  
  template<typename T1, typename T2>
  inline explicit Row(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);
  
  template<typename T1> inline                  Row(const BaseCube<eT,T1>& X);
  template<typename T1> inline const Row& operator=(const BaseCube<eT,T1>& X);
  
  inline                  Row(const subview_cube<eT>& X);
  inline const Row& operator=(const subview_cube<eT>& X);
  
  inline mat_injector<Row> operator<<(const eT val);
  
  arma_inline eT& col(const u32 col_num);
  arma_inline eT  col(const u32 col_num) const;
  
  arma_inline       subview_row<eT> cols(const u32 in_col1, const u32 in_col2);
  arma_inline const subview_row<eT> cols(const u32 in_col1, const u32 in_col2) const;
  
  arma_inline       subview_row<eT> subvec(const u32 in_col1, const u32 in_col2);
  arma_inline const subview_row<eT> subvec(const u32 in_col1, const u32 in_col2) const;
  
  
  inline void shed_col (const u32 col_num);
  inline void shed_cols(const u32 in_col1, const u32 in_col2);
  
                        inline void insert_cols(const u32 col_num, const u32 N, const bool set_to_zero = true);
  template<typename T1> inline void insert_cols(const u32 col_num, const Base<eT,T1>& X);
  
  
  typedef       eT*       row_iterator;
  typedef const eT* const_row_iterator;
  
  inline       row_iterator begin_row(const u32 row_num);
  inline const_row_iterator begin_row(const u32 row_num) const;
  
  inline       row_iterator end_row  (const u32 row_num);
  inline const_row_iterator end_row  (const u32 row_num) const;
  
  
  template<u32 fixed_n_elem>
  class fixed : public Row<eT>
    {
    private:
    
    arma_aligned eT mem_local_extra[ (fixed_n_elem > arma_config::mat_prealloc) ? fixed_n_elem : 1 ];
    
    arma_inline void mem_setup();
    
    
    public:
    
    inline fixed() { mem_setup(); }
    
    inline                fixed(const char*        text) { mem_setup(); Row<eT>::operator=(text);               }
    inline const Row& operator=(const char*        text) {              Row<eT>::operator=(text); return *this; }
    inline                fixed(const std::string& text) { mem_setup(); Row<eT>::operator=(text);               }
    inline const Row& operator=(const std::string& text) {              Row<eT>::operator=(text); return *this; }
    
    inline const Row& operator=(const eT val) { Row<eT>::operator=(val); return *this; }
    
    template<typename T1>
    inline fixed(const Base<eT,T1>& A) { mem_setup(); Row<eT>::operator=(A.get_ref()); }
    
    template<typename T1>
    inline const Row& operator=(const Base<eT,T1>& A) { Row<eT>::operator=(A.get_ref()); return *this; }
    
    template<typename T1, typename T2>
    inline explicit fixed(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B) { mem_setup(); Row<eT>::init(A,B); }
    
    inline fixed(      eT* aux_mem, const bool copy_aux_mem = true);
    inline fixed(const eT* aux_mem);
    
    inline                fixed(const subview_cube<eT>& X) { mem_setup(); Row<eT>::operator=(X);               }
    inline const Row& operator=(const subview_cube<eT>& X) {              Row<eT>::operator=(X); return *this; }
    };
  
  
  #ifdef ARMA_EXTRA_ROW_PROTO
    #include ARMA_INCFILE_WRAP(ARMA_EXTRA_ROW_PROTO)
  #endif
  
  };



//! @}
