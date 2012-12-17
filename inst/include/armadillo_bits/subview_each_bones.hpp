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


//! \addtogroup subview_each
//! @{



template<typename parent, unsigned int mode>
class subview_each_common
  {
  public:
  
  typedef typename parent::elem_type eT;
  
  
  protected:
  
  parent& p;
  
  arma_inline subview_each_common(parent& in_p);
  
  arma_inline const Mat<typename parent::elem_type>& get_mat_ref_helper(const Mat    <typename parent::elem_type>& X) const;
  arma_inline const Mat<typename parent::elem_type>& get_mat_ref_helper(const subview<typename parent::elem_type>& X) const;
  
  arma_inline const Mat<typename parent::elem_type>& get_mat_ref() const;
  
  inline void check_size(const Mat<typename parent::elem_type>& A) const;
  
  arma_cold inline const std::string incompat_size_string(const Mat<typename parent::elem_type>& A) const;
  
  
  private:
  
  subview_each_common();
  };




template<typename parent, unsigned int mode>
class subview_each1 : public subview_each_common<parent, mode>
  {
  protected:
  
  arma_inline subview_each1(parent& in_p);
  
  
  public:
  
  typedef typename parent::elem_type eT;
  
  inline ~subview_each1();
  
  // deliberately returning void
  template<typename T1> inline void operator=  (const Base<eT,T1>& x);
  template<typename T1> inline void operator+= (const Base<eT,T1>& x);
  template<typename T1> inline void operator-= (const Base<eT,T1>& x);
  template<typename T1> inline void operator%= (const Base<eT,T1>& x);
  template<typename T1> inline void operator/= (const Base<eT,T1>& x);
  
  
  private:
  
  friend class Mat<eT>;
  friend class subview<eT>;
  
  subview_each1();
  };



template<typename parent, unsigned int mode, typename TB>
class subview_each2 : public subview_each_common<parent, mode>
  {
  protected:
  
  const Base<uword, TB>& base_indices;
  
  inline subview_each2(parent& in_p, const Base<uword, TB>& in_indices);
  
  inline void check_indices(const Mat<uword>& indices) const;
  
  
  public:
  
  typedef typename parent::elem_type eT;
  
  inline ~subview_each2();
  
  // deliberately returning void
  template<typename T1> inline void operator=  (const Base<eT,T1>& x);
  template<typename T1> inline void operator+= (const Base<eT,T1>& x);
  template<typename T1> inline void operator-= (const Base<eT,T1>& x);
  template<typename T1> inline void operator%= (const Base<eT,T1>& x);
  template<typename T1> inline void operator/= (const Base<eT,T1>& x);
  
  // TODO: add handling of scalars
  
  
  private:
  
  friend class Mat<eT>;
  friend class subview<eT>;
  
  subview_each2();
  };



//! @}
