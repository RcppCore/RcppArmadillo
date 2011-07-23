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


//! \addtogroup Col
//! @{

//! Class for column vectors (matrices with only one column)

template<typename eT>
class Col : public Mat<eT>
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
  
  template<typename T1> inline                  Col(const BaseCube<eT,T1>& X);
  template<typename T1> inline const Col& operator=(const BaseCube<eT,T1>& X);
  
  inline                  Col(const subview_cube<eT>& X);
  inline const Col& operator=(const subview_cube<eT>& X);
  
  inline mat_injector<Col> operator<<(const eT val);
  
  arma_inline eT& row(const u32 row_num);
  arma_inline eT  row(const u32 row_num) const;
  
  arma_inline       subview_col<eT> rows(const u32 in_row1, const u32 in_row2);
  arma_inline const subview_col<eT> rows(const u32 in_row1, const u32 in_row2) const;
  
  arma_inline       subview_col<eT> subvec(const u32 in_row1, const u32 in_row2);
  arma_inline const subview_col<eT> subvec(const u32 in_row1, const u32 in_row2) const;
  
  arma_inline       subview_col<eT> subvec(const span& row_span);
  arma_inline const subview_col<eT> subvec(const span& row_span) const;
  
  //arma_inline       subview_col<eT> operator()(const span& row_span);
  //arma_inline const subview_col<eT> operator()(const span& row_span) const;
  
  
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
    
    static const bool use_extra = (fixed_n_elem > arma_config::mat_prealloc);
    
    arma_aligned eT mem_local_extra[ (use_extra) ? fixed_n_elem : 1 ];
    
    arma_inline void mem_setup();
    arma_inline void change_to_row();
    
    
    public:
    
    static const u32 n_rows = fixed_n_elem;
    static const u32 n_cols = 1;
    static const u32 n_elem = fixed_n_elem;
    
    arma_inline fixed();
    arma_inline fixed(const fixed<fixed_n_elem>& X);
         inline fixed(const subview_cube<eT>& X);
    
    template<typename T1>              inline fixed(const Base<eT,T1>& A);
    template<typename T1, typename T2> inline fixed(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);
    
    inline fixed(      eT* aux_mem, const bool copy_aux_mem = true);
    inline fixed(const eT* aux_mem);
    
    inline fixed(const char*        text);
    inline fixed(const std::string& text);
    
    template<typename T1> inline const Col& operator=(const Base<eT,T1>& A);
    
    inline const Col& operator=(const eT val);
    inline const Col& operator=(const char*        text);
    inline const Col& operator=(const std::string& text);
    inline const Col& operator=(const subview_cube<eT>& X);
    
    inline       subview_row<eT> operator()(const u32   row_num,  const span& col_span);
    inline const subview_row<eT> operator()(const u32   row_num,  const span& col_span) const;
    
    inline       subview_col<eT> operator()(const span& row_span, const u32   col_num );
    inline const subview_col<eT> operator()(const span& row_span, const u32   col_num ) const;
    
    inline       subview<eT>     operator()(const span& row_span, const span& col_span);
    inline const subview<eT>     operator()(const span& row_span, const span& col_span) const;
    
    arma_inline arma_warn_unused eT& operator[] (const u32 i);
    arma_inline arma_warn_unused eT  operator[] (const u32 i) const;
    arma_inline arma_warn_unused eT& at         (const u32 i);
    arma_inline arma_warn_unused eT  at         (const u32 i) const;
    arma_inline arma_warn_unused eT& operator() (const u32 i);
    arma_inline arma_warn_unused eT  operator() (const u32 i) const;
    
    arma_inline arma_warn_unused eT& at         (const u32 in_row, const u32 in_col);
    arma_inline arma_warn_unused eT  at         (const u32 in_row, const u32 in_col) const;
    arma_inline arma_warn_unused eT& operator() (const u32 in_row, const u32 in_col);
    arma_inline arma_warn_unused eT  operator() (const u32 in_row, const u32 in_col) const;
    
    arma_hot inline const Col<eT>& fill(const eT val);
    arma_hot inline const Col<eT>& zeros();
    arma_hot inline const Col<eT>& ones();
    };
  
  
  #ifdef ARMA_EXTRA_COL_PROTO
    #include ARMA_INCFILE_WRAP(ARMA_EXTRA_COL_PROTO)
  #endif
  
  };



//! @}
