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


//! \addtogroup diagmat_proxy
//! @{



template<typename T1>
class diagmat_proxy
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline diagmat_proxy(const Base<typename T1::elem_type,T1>& X)
    : P       (X.get_ref())
    , P_is_vec( (P.n_rows == 1) || (P.n_cols == 1) )
    , n_elem  ( P_is_vec ? P.n_elem : P.n_rows )
    {
    arma_extra_debug_sigprint();
    
    arma_debug_check( (P_is_vec == false) && (P.n_rows != P.n_cols), "diagmat(): only vectors and square matrices are accepted" );
    }
  
  
  arma_inline elem_type operator[] (const u32 i)                  const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }
  
  
  const Proxy<T1> P;
  const bool      P_is_vec;
  const u32       n_elem;
  };



template<typename eT>
class diagmat_proxy< Mat<eT> >
  {
  public:
  
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  
  inline diagmat_proxy(const Mat<eT>& X)
    : P(X)
    , P_is_vec( (P.n_rows == 1) || (P.n_cols == 1) )
    , n_elem( P_is_vec ? P.n_elem : P.n_rows )
    {
    arma_extra_debug_sigprint();
    arma_debug_check( (P_is_vec == false) && (P.n_rows != P.n_cols), "diagmat(): only vectors and square matrices are accepted" );
    }
  
  
  arma_inline elem_type operator[] (const u32 i)                  const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }

  const Mat<eT>& P;
  const bool     P_is_vec;
  const u32      n_elem;
  };



template<typename eT>
class diagmat_proxy< Row<eT> >
  {
  public:
  
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  
  inline diagmat_proxy(const Row<eT>& X)
    : P(X)
    , P_is_vec(true)
    , n_elem(P.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  
  arma_inline elem_type operator[] (const u32 i)                  const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }


  const Row<eT>& P;
  const bool     P_is_vec;
  const u32      n_elem;
  };



template<typename eT>
class diagmat_proxy< Col<eT> >
  {
  public:
  
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  
  inline diagmat_proxy(const Col<eT>& X)
    : P(X)
    , P_is_vec(true)
    , n_elem(P.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  
  arma_inline elem_type operator[] (const u32 i)                  const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }
  

  const Col<eT>& P;
  const bool     P_is_vec;
  const u32      n_elem;
  };



template<typename T1>
class diagmat_proxy_check
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline diagmat_proxy_check(const Base<typename T1::elem_type,T1>& X, const Mat<typename T1::elem_type>& out)
    : P(X.get_ref())
    , P_is_vec( (P.n_rows == 1) || (P.n_cols == 1) )
    , n_elem( P_is_vec ? P.n_elem : P.n_rows )
    {
    arma_extra_debug_sigprint();
    arma_debug_check( (P_is_vec == false) && (P.n_rows != P.n_cols), "diagmat(): only vectors and square matrices are accepted" );
    }
  
  
  arma_inline elem_type operator[] (const u32 i)                  const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }
  

  const Mat<elem_type> P;
  const bool           P_is_vec;
  const u32            n_elem;
  };



template<typename eT>
class diagmat_proxy_check< Mat<eT> >
  {
  public:
  
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  
  inline diagmat_proxy_check(const Mat<eT>& X, const Mat<eT>& out)
    : P_local ( (&X == &out) ? new Mat<eT>(X) : 0  )
    , P       ( (&X == &out) ? (*P_local)     : X  )
    , P_is_vec( (P.n_rows == 1) || (P.n_cols == 1) )
    , n_elem  ( P_is_vec ? P.n_elem : P.n_rows     )
    {
    arma_extra_debug_sigprint();
    
    arma_debug_check( (P_is_vec == false) && (P.n_rows != P.n_cols), "diagmat(): only vectors and square matrices are accepted" );
    }
  
  inline ~diagmat_proxy_check()
    {
    if(P_local)
      {
      delete P_local;
      }
    }
  
  
  arma_inline elem_type operator[] (const u32 i)                  const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }
  

  const Mat<eT>* P_local;
  const Mat<eT>& P;
  const bool     P_is_vec;
  const u32      n_elem;
  };



template<typename eT>
class diagmat_proxy_check< Row<eT> >
  {
  public:
  
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline diagmat_proxy_check(const Row<eT>& X, const Mat<eT>& out)
    : P_local ( (&X == reinterpret_cast<const Row<eT>*>(&out)) ? new Row<eT>(X) : 0 )
    , P       ( (&X == reinterpret_cast<const Row<eT>*>(&out)) ? (*P_local)     : X )
    , P_is_vec(true)
    , n_elem  (P.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline ~diagmat_proxy_check()
    {
    if(P_local)
      {
      delete P_local;
      }
    }
  
  
  arma_inline elem_type operator[] (const u32 i)                  const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }
  
  
  const Row<eT>* P_local;
  const Row<eT>& P;
  const bool     P_is_vec;
  const u32      n_elem;
  };






template<typename eT>
class diagmat_proxy_check< Col<eT> >
  {
  public:
  
  typedef          eT                              elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline diagmat_proxy_check(const Col<eT>& X, const Mat<eT>& out)
    : P_local ( (&X == reinterpret_cast<const Col<eT>*>(&out)) ? new Col<eT>(X) : 0 )
    , P       ( (&X == reinterpret_cast<const Col<eT>*>(&out)) ? (*P_local)     : X )
    , P_is_vec(true)
    , n_elem  (P.n_elem)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline ~diagmat_proxy_check()
    {
    if(P_local)
      {
      delete P_local;
      }
    }
  
  
  arma_inline elem_type operator[] (const u32 i)                  const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const u32 row, const u32 col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }
  
  
  const Col<eT>* P_local;
  const Col<eT>& P;
  const bool     P_is_vec;
  const u32      n_elem;
  };



//! @}
