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


//! \addtogroup Proxy
//! @{



template<typename T1>
class Proxy
  {
  public:
  inline Proxy(const T1& A)
    {
    arma_type_check< is_arma_type<T1>::value == false >::apply();
    }
  };



// ea_type is the "element accessor" type,
// which can provide access to elements via operator[]

template<typename eT>
class Proxy< Mat<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  
  arma_aligned const Mat<eT>& Q;
  
  inline explicit Proxy(const Mat<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename eT>
class Proxy< Col<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Col<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  
  arma_aligned const Col<eT>& Q;
  
  inline explicit Proxy(const Col<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename eT>
class Proxy< Row<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Row<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  
  arma_aligned const Row<eT>& Q;
  
  inline explicit Proxy(const Row<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename T1, typename op_type>
class Proxy< Op<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  
  arma_aligned const Mat<elem_type> Q;
  
  inline explicit Proxy(const Op<T1, op_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename T1, typename T2, typename glue_type>
class Proxy< Glue<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  
  arma_aligned const Mat<elem_type> Q;
  
  inline explicit Proxy(const Glue<T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }

  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename eT>
class Proxy< subview<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview<eT>                              stored_type;
  typedef const subview<eT>&                       ea_type;
  
  arma_aligned const subview<eT>& Q;
  
  inline explicit Proxy(const subview<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q; }
  };




template<typename eT>
class Proxy< diagview<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef diagview<eT>                             stored_type;
  typedef const diagview<eT>&                      ea_type;
  
  arma_aligned const diagview<eT>& Q;
  
  inline explicit Proxy(const diagview<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q; }
  };




template<typename T1, typename eop_type>
class Proxy< eOp<T1, eop_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eOp<T1, eop_type>                        stored_type;
  typedef const eOp<T1, eop_type>&                 ea_type;
  
  arma_aligned const eOp<T1, eop_type>& Q;
  
  inline explicit Proxy(const eOp<T1, eop_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.get_n_rows(); }
  arma_inline u32 get_n_cols() const { return Q.get_n_cols(); }
  arma_inline u32 get_n_elem() const { return Q.get_n_elem(); }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q; }
  };



template<typename T1, typename T2, typename eglue_type>
class Proxy< eGlue<T1, T2, eglue_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eGlue<T1, T2, eglue_type>                stored_type;
  typedef const eGlue<T1, T2, eglue_type>&         ea_type;
  
  arma_aligned const eGlue<T1, T2, eglue_type>& Q;
  
  inline explicit Proxy(const eGlue<T1, T2, eglue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.get_n_rows(); }
  arma_inline u32 get_n_cols() const { return Q.get_n_cols(); }
  arma_inline u32 get_n_elem() const { return Q.get_n_elem(); }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q; }
  };



template<typename out_eT, typename T1, typename op_type>
class Proxy< mtOp<out_eT, T1, op_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Mat<out_eT>                  stored_type;
  typedef          const elem_type*             ea_type;
  
  arma_aligned const Mat<out_eT> Q;
  
  inline explicit Proxy(const mtOp<out_eT, T1, op_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];          }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row,col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename out_eT, typename T1, typename T2, typename glue_type>
class Proxy< mtGlue<out_eT, T1, T2, glue_type > >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Mat<out_eT>                  stored_type;
  typedef          const elem_type*             ea_type;
  
  arma_aligned const Mat<out_eT> Q;
  
  inline explicit Proxy(const mtGlue<out_eT, T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows() const { return Q.n_rows; }
  arma_inline u32 get_n_cols() const { return Q.n_cols; }
  arma_inline u32 get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];          }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row,col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



//! @}
