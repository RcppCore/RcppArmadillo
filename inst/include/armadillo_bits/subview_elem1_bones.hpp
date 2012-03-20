// Copyright (C) 2010-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup subview_elem1
//! @{



template<typename eT, typename T1>
class subview_elem1 : public Base<eT, subview_elem1<eT,T1> >
  {
  public:    
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  static const bool is_row = false;
  static const bool is_col = true;
  
  arma_aligned const Mat<eT>&        m;
  arma_aligned const Base<uword,T1>& a;
  
  
  protected:
  
  arma_inline subview_elem1(const Mat<eT>& in_m, const Base<uword,T1>& in_a);
  
  
  public:
  
  inline ~subview_elem1();
  
  template<typename op_type>              inline void inplace_op(const eT                    val);
  template<typename op_type, typename T2> inline void inplace_op(const subview_elem1<eT,T2>& x  );
  template<typename op_type, typename T2> inline void inplace_op(const Base<eT,T2>&          x  );
  
  arma_inline const Op<subview_elem1<eT,T1>,op_htrans>  t() const;
  arma_inline const Op<subview_elem1<eT,T1>,op_htrans> ht() const;
  arma_inline const Op<subview_elem1<eT,T1>,op_strans> st() const;
  
  inline void fill(const eT val);
  inline void zeros();
  inline void ones();
  
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  
  // deliberately returning void
  template<typename T2> inline void operator_equ(const subview_elem1<eT,T2>& x);
  template<typename T2> inline void operator=   (const subview_elem1<eT,T2>& x);
                        inline void operator=   (const subview_elem1<eT,T1>& x);
  template<typename T2> inline void operator+=  (const subview_elem1<eT,T2>& x);
  template<typename T2> inline void operator-=  (const subview_elem1<eT,T2>& x);
  template<typename T2> inline void operator%=  (const subview_elem1<eT,T2>& x);
  template<typename T2> inline void operator/=  (const subview_elem1<eT,T2>& x);
  
  template<typename T2> inline void operator=  (const Base<eT,T2>& x);
  template<typename T2> inline void operator+= (const Base<eT,T2>& x);
  template<typename T2> inline void operator-= (const Base<eT,T2>& x);
  template<typename T2> inline void operator%= (const Base<eT,T2>& x);
  template<typename T2> inline void operator/= (const Base<eT,T2>& x);
  
  inline static void extract(Mat<eT>& out, const subview_elem1& in);
  
  template<typename op_type> inline static void mat_inplace_op(Mat<eT>& out, const subview_elem1& in);
  
  inline static void  plus_inplace(Mat<eT>& out, const subview_elem1& in);
  inline static void minus_inplace(Mat<eT>& out, const subview_elem1& in);
  inline static void schur_inplace(Mat<eT>& out, const subview_elem1& in);
  inline static void   div_inplace(Mat<eT>& out, const subview_elem1& in);
  
  
  
  private:
  
  friend class Mat<eT>;
  subview_elem1();
  };



//! @}
