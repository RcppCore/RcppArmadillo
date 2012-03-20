// Copyright (C) 2012 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup subview_elem2
//! @{



template<typename eT, typename T1, typename T2>
class subview_elem2 : public Base<eT, subview_elem2<eT,T1,T2> >
  {
  public:    
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool is_row = false;
  static const bool is_col = false;
  
  arma_aligned const Mat<eT>& m;
  
  arma_aligned const Base<uword,T1>& base_ri;
  arma_aligned const Base<uword,T2>& base_ci;
  
  const bool all_rows;
  const bool all_cols;
  
  
  protected:
  
  arma_inline subview_elem2(const Mat<eT>& in_m, const Base<uword,T1>& in_ri, const Base<uword,T1>& in_ci, const bool in_all_rows, const bool in_all_cols);
  
  
  public:
  
  inline ~subview_elem2();
  
  template<typename op_type>
  inline void inplace_op(const eT val);
  
  template<typename op_type, typename expr>
  inline void inplace_op(const Base<eT,expr>& x);
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  
  // deliberately returning void
  template<typename T3, typename T4> inline void operator_equ(const subview_elem2<eT,T3,T4>& x);
  template<typename T3, typename T4> inline void operator=   (const subview_elem2<eT,T3,T4>& x);
                                     inline void operator=   (const subview_elem2<eT,T1,T2>& x);
  
  template<typename T3, typename T4> inline void operator+=  (const subview_elem2<eT,T3,T4>& x);
  template<typename T3, typename T4> inline void operator-=  (const subview_elem2<eT,T3,T4>& x);
  template<typename T3, typename T4> inline void operator%=  (const subview_elem2<eT,T3,T4>& x);
  template<typename T3, typename T4> inline void operator/=  (const subview_elem2<eT,T3,T4>& x);
  
  template<typename expr> inline void operator=  (const Base<eT,expr>& x);
  template<typename expr> inline void operator+= (const Base<eT,expr>& x);
  template<typename expr> inline void operator-= (const Base<eT,expr>& x);
  template<typename expr> inline void operator%= (const Base<eT,expr>& x);
  template<typename expr> inline void operator/= (const Base<eT,expr>& x);
  
  inline static void extract(Mat<eT>& out, const subview_elem2& in);
  
  inline static void  plus_inplace(Mat<eT>& out, const subview_elem2& in);
  inline static void minus_inplace(Mat<eT>& out, const subview_elem2& in);
  inline static void schur_inplace(Mat<eT>& out, const subview_elem2& in);
  inline static void   div_inplace(Mat<eT>& out, const subview_elem2& in);
  
  
  
  private:
  
  friend class Mat<eT>;
  subview_elem2();
  };



//! @}
