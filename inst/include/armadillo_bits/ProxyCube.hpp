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


//! \addtogroup ProxyCube
//! @{



template<typename T1>
class ProxyCube
  {
  public:
  inline ProxyCube(const T1& A)
    {
    arma_type_check< is_arma_type<T1>::value == false >::apply();
    }
  };



template<typename eT>
class ProxyCube< Cube<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<eT>                                 stored_type;
  
  
  const Cube<eT>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem_slice;
  const u32 n_slices;
  const u32 n_elem;
  
  inline explicit ProxyCube(const Cube<eT>& A)
    : Q           (A)
    , n_rows      (A.n_rows)
    , n_cols      (A.n_cols)
    , n_elem_slice(A.n_elem_slice)
    , n_slices    (A.n_slices)
    , n_elem      (A.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  inline explicit ProxyCube(const u32 in_n_rows, const u32 in_n_cols, const u32 in_n_slices)
    : Q           (Q)
    , n_rows      (in_n_rows)
    , n_cols      (in_n_cols)
    , n_elem_slice(in_n_rows*in_n_cols)
    , n_slices    (in_n_slices)
    , n_elem      (in_n_rows*in_n_cols*in_n_slices)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  };



template<typename T1, typename op_type>
class ProxyCube< OpCube<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<elem_type>                          stored_type;
  
  const Cube<elem_type> Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem_slice;
  const u32 n_slices;
  const u32 n_elem;
  
  inline explicit ProxyCube(const OpCube<T1, op_type>& A)
    : Q           (A)
    , n_rows      (Q.n_rows)
    , n_cols      (Q.n_cols)
    , n_elem_slice(Q.n_elem_slice)
    , n_slices    (Q.n_slices)
    , n_elem      (Q.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  };



template<typename T1, typename T2, typename glue_type>
class ProxyCube< GlueCube<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<elem_type>                          stored_type;
  
  const Cube<elem_type> Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem_slice;
  const u32 n_slices;
  const u32 n_elem;
  
  inline explicit ProxyCube(const GlueCube<T1, T2, glue_type>& A)
    : Q           (A)
    , n_rows      (Q.n_rows)
    , n_cols      (Q.n_cols)
    , n_elem_slice(Q.n_elem_slice)
    , n_slices    (Q.n_slices)
    , n_elem      (Q.n_elem)
    {
    arma_extra_debug_sigprint();
    }

  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  };



template<typename eT>
class ProxyCube< subview_cube<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_cube<eT>                         stored_type;
  
  const subview_cube<eT>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem_slice;
  const u32 n_slices;
  const u32 n_elem;
  
  inline explicit ProxyCube(const subview_cube<eT>& A)
    : Q           (A)
    , n_rows      (A.n_rows)
    , n_cols      (A.n_cols)
    , n_elem_slice(A.n_elem_slice)
    , n_slices    (A.n_slices)
    , n_elem      (A.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  };




template<typename T1, typename eop_type>
class ProxyCube< eOpCube<T1, eop_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eOpCube<T1, eop_type>                    stored_type;
  
  const eOpCube<T1, eop_type>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem_slice;
  const u32 n_slices;
  const u32 n_elem;
  
  inline explicit ProxyCube(const eOpCube<T1, eop_type>& A)
    : Q           (A)
    , n_rows      (A.P.n_rows)
    , n_cols      (A.P.n_cols)
    , n_elem_slice(A.P.n_elem_slice)
    , n_slices    (A.P.n_slices)
    , n_elem      (A.P.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return eop_type::get_elem(Q, i);               }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return eop_type::get_elem(Q, row, col, slice); }
  };



template<typename T1, typename T2, typename eglue_type>
class ProxyCube< eGlueCube<T1, T2, eglue_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eGlueCube<T1, T2, eglue_type>            stored_type;
  
  const eGlueCube<T1, T2, eglue_type>& Q;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem_slice;
  const u32 n_slices;
  const u32 n_elem;
  
  inline explicit ProxyCube(const eGlueCube<T1, T2, eglue_type>& A)
    : Q           (A)
    , n_rows      (A.P1.n_rows)
    , n_cols      (A.P1.n_cols)
    , n_elem_slice(A.P1.n_elem_slice)
    , n_slices    (A.P1.n_slices)
    , n_elem      (A.P1.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return eglue_type::get_elem(Q, i);               }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return eglue_type::get_elem(Q, row, col, slice); }
  };



//! @}
