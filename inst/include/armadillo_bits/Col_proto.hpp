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


//! \addtogroup Col
//! @{

//! Class for column vectors (matrices with only column)

template<typename eT>
class Col : public Mat<eT>, public BaseVec< eT, Col<eT> >
  {
  public:
  
  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  
  
  inline          Col();
  inline explicit Col(const u32 n_elem);
  inline          Col(const u32 in_rows, const u32 in_cols);
  
  inline                  Col(const char*        text);
  inline const Col& operator=(const char*        text);
  inline                  Col(const std::string& text);
  inline const Col& operator=(const std::string& text);
  
  inline const Col& operator=(const eT val);
    
  template<typename T1> inline                   Col(const Base<eT,T1>& X);
  template<typename T1> inline const Col&  operator=(const Base<eT,T1>& X);
  
  inline Col(      eT* aux_mem, const u32 aux_length, const bool copy_aux_mem = true, const bool strict = true);
  inline Col(const eT* aux_mem, const u32 aux_length);
  
  template<typename T1, typename T2>
  inline explicit Col(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);
  
  inline                  Col(const subview_cube<eT>& X);
  inline const Col& operator=(const subview_cube<eT>& X);
  
  inline mat_injector<Col> operator<<(const eT val);
  
  arma_inline eT& row(const u32 row_num);
  arma_inline eT  row(const u32 row_num) const;
  
  arma_inline       subview_col<eT> rows(const u32 in_row1, const u32 in_row2);
  arma_inline const subview_col<eT> rows(const u32 in_row1, const u32 in_row2) const;
  
  
  inline void shed_row (const u32 row_num);
  inline void shed_rows(const u32 in_row1, const u32 in_row2);
  
                        inline void insert_rows(const u32 row_num, const u32 N, const bool set_to_zero = true);
  template<typename T1> inline void insert_rows(const u32 row_num, const Base<eT,T1>& X);
  
  
  typedef       eT*       row_iterator;
  typedef const eT* const_row_iterator;
  
  inline       row_iterator begin_row(const u32 row_num);
  inline const_row_iterator begin_row(const u32 row_num) const;
  
  inline       row_iterator end_row  (const u32 row_num);
  inline const_row_iterator end_row  (const u32 row_num) const;
  
  
  template<u32 fixed_n_elem>
  class fixed : public Col<eT>
    {
    private:
    
    arma_aligned eT mem_local_extra[ ( fixed_n_elem > Mat_prealloc::mem_n_elem ) ? fixed_n_elem : 1 ];
    
    arma_inline void mem_setup();
    arma_inline void swap_rows_cols() { access::rw(Mat<eT>::n_cols) = fixed_n_elem; access::rw(Mat<eT>::n_rows) = 1; }
    
    public:
    
    inline fixed() { mem_setup(); }
    
    inline                fixed(const char*        text) { mem_setup(); swap_rows_cols(); Col<eT>::operator=(text);               }
    inline const Col& operator=(const char*        text) {              swap_rows_cols(); Col<eT>::operator=(text); return *this; }
    inline                fixed(const std::string& text) { mem_setup(); swap_rows_cols(); Col<eT>::operator=(text);               }
    inline const Col& operator=(const std::string& text) {              swap_rows_cols(); Col<eT>::operator=(text); return *this; }
    
    inline const Col& operator=(const eT val) { Col<eT>::operator=(val); return *this; }
    
    template<typename T1>
    inline fixed(const Base<eT,T1>& A) { mem_setup(); Col<eT>::operator=(A.get_ref()); }
    
    template<typename T1>
    inline const Col& operator=(const Base<eT,T1>& A) { Col<eT>::operator=(A.get_ref()); return *this; }
    
    template<typename T1, typename T2>
    inline explicit fixed(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B) { mem_setup(); Col<eT>::init(A,B); }
    
    inline                fixed(const subview_cube<eT>& X) { mem_setup(); Col<eT>::operator=(X);               }
    inline const Col& operator=(const subview_cube<eT>& X) {              Col<eT>::operator=(X); return *this; }
    };
  
  
  #ifdef ARMA_EXTRA_COL_PROTO
    #include ARMA_INCFILE_WRAP(ARMA_EXTRA_COL_PROTO)
  #endif
  
  };



//! @}
