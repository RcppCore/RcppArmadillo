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



template<typename eT>
class Proxy< Mat<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  
  const Mat<eT>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const Mat<eT>& A)
    : Q(A)
    , n_rows(A.n_rows)
    , n_cols(A.n_cols)
    , n_elem(A.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  inline explicit Proxy(const u32 in_n_rows, const u32 in_n_cols)
    : Q(Q)
    , n_rows(in_n_rows)
    , n_cols(in_n_cols)
    , n_elem(in_n_rows*in_n_cols)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  };



template<typename eT>
class Proxy< Col<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Col<eT>                                  stored_type;
  
  const Col<eT>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const Col<eT>& A)
    : Q(A)
    , n_rows(A.n_rows)
    , n_cols(A.n_cols)
    , n_elem(A.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  inline explicit Proxy(const u32 in_n_rows, const u32 in_n_cols)
    : Q(Q)
    , n_rows(in_n_rows)
    , n_cols(in_n_cols)
    , n_elem(in_n_rows*in_n_cols)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  };



template<typename eT>
class Proxy< Row<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Row<eT>                                  stored_type;
  
  const Row<eT>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const Row<eT>& A)
    : Q(A)
    , n_rows(A.n_rows)
    , n_cols(A.n_cols)
    , n_elem(A.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  inline explicit Proxy(const u32 in_n_rows, const u32 in_n_cols)
    : Q(Q)
    , n_rows(in_n_rows)
    , n_cols(in_n_cols)
    , n_elem(in_n_rows*in_n_cols)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  };



template<typename T1, typename op_type>
class Proxy< Op<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  
  const Mat<elem_type> Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const Op<T1, op_type>& A)
    : Q(A)
    , n_rows(Q.n_rows)
    , n_cols(Q.n_cols)
    , n_elem(Q.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  };



template<typename T1, typename T2, typename glue_type>
class Proxy< Glue<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  
  const Mat<elem_type> Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const Glue<T1, T2, glue_type>& A)
    : Q(A)
    , n_rows(Q.n_rows)
    , n_cols(Q.n_cols)
    , n_elem(Q.n_elem)
    {
    arma_extra_debug_sigprint();
    }

  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  };



template<typename eT>
class Proxy< subview<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview<eT>                              stored_type;
  
  const subview<eT>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const subview<eT>& A)
    : Q(A)
    , n_rows(A.n_rows)
    , n_cols(A.n_cols)
    , n_elem(A.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  };




template<typename eT>
class Proxy< diagview<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef diagview<eT>                             stored_type;
  
  const diagview<eT>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const diagview<eT>& A)
    : Q(A)
    , n_rows(A.n_rows)
    , n_cols(A.n_cols)
    , n_elem(A.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return Q[i];           }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return Q.at(row, col); }
  };




template<typename T1, typename eop_type>
class Proxy< eOp<T1, eop_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eOp<T1, eop_type>                        stored_type;
  
  const eOp<T1, eop_type>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const eOp<T1, eop_type>& A)
    : Q(A)
    , n_rows(A.P.n_rows)
    , n_cols(A.P.n_cols)
    , n_elem(A.P.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return eop_type::get_elem(Q, i);       }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return eop_type::get_elem(Q, row,col); }
  };



template<typename T1, typename T2, typename eglue_type>
class Proxy< eGlue<T1, T2, eglue_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eGlue<T1, T2, eglue_type>                stored_type;
  
  const eGlue<T1, T2, eglue_type>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  inline explicit Proxy(const eGlue<T1, T2, eglue_type>& A)
    : Q(A)
    , n_rows(A.P1.n_rows)
    , n_cols(A.P1.n_cols)
    , n_elem(A.P1.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                  const { return eglue_type::get_elem(Q, i);        }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return eglue_type::get_elem(Q, row, col); }
  };



//! @}
