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


//! \addtogroup diagmat_proxy
//! @{



template<typename T1>
class diagmat_proxy
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  
  inline diagmat_proxy(const Base<typename T1::elem_type,T1>& X)
    : P       ( X.get_ref() )
    , P_is_vec( (P.get_n_rows() == 1) || (P.get_n_cols() == 1) )
    , n_elem  ( P_is_vec ? P.get_n_elem() : (std::min)(P.get_n_elem(), P.get_n_rows()) )
    {
    arma_extra_debug_sigprint();
    
    arma_debug_check
      (
      (P_is_vec == false) && (P.get_n_rows() != P.get_n_cols()),
      "diagmat(): only vectors and square matrices are accepted"
      );
    }
  
  
  arma_inline
  elem_type
  operator[](const uword i) const
    {
    if( (Proxy<T1>::prefer_at_accessor == true) || (P_is_vec == false) )
      {
      return P.at(i,i);
      }
    else
      {
      return P[i];
      }
    }
  
  
  arma_inline
  elem_type
  at(const uword row, const uword col) const
    {
    if(row == col)
      {
      if( (Proxy<T1>::prefer_at_accessor == true) || (P_is_vec == false) )
        {
        return P.at(row,row);
        }
      else
        {
        return P[row];
        }
      }
    else
      {
      return elem_type(0);
      }
    }
  
  
  const Proxy<T1> P;
  const bool      P_is_vec;
  const uword     n_elem;
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
    , n_elem( P_is_vec ? P.n_elem : (std::min)(P.n_elem, P.n_rows) )
    {
    arma_extra_debug_sigprint();
    
    arma_debug_check
      (
      (P_is_vec == false) && (P.n_rows != P.n_cols),
      "diagmat(): only vectors and square matrices are accepted"
      );
    }
  
  
  arma_inline elem_type operator[] (const uword i)                    const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const uword row, const uword col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }

  const Mat<eT>& P;
  const bool     P_is_vec;
  const uword    n_elem;
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
  
  
  arma_inline elem_type operator[] (const uword i)                    const { return P[i];                                                                }
  arma_inline elem_type at         (const uword row, const uword col) const { return (row == col) ? P[row] : elem_type(0); }


  const Row<eT>& P;
  const bool     P_is_vec;
  const uword    n_elem;
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
  
  
  arma_inline elem_type operator[] (const uword i)                    const { return P[i];                                 }
  arma_inline elem_type at         (const uword row, const uword col) const { return (row == col) ? P[row] : elem_type(0); }
  

  const Col<eT>& P;
  const bool     P_is_vec;
  const uword    n_elem;
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
    , n_elem( P_is_vec ? P.n_elem : (std::min)(P.n_elem, P.n_rows) )
    {
    arma_extra_debug_sigprint();
    arma_ignore(out);
    
    arma_debug_check
      (
      (P_is_vec == false) && (P.n_rows != P.n_cols),
      "diagmat(): only vectors and square matrices are accepted"
      );
    }
  
  
  arma_inline elem_type operator[] (const uword i)                    const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const uword row, const uword col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }
  

  const Mat<elem_type> P;
  const bool           P_is_vec;
  const uword          n_elem;
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
    , n_elem  ( P_is_vec ? P.n_elem : (std::min)(P.n_elem, P.n_rows) )
    {
    arma_extra_debug_sigprint();
    
    arma_debug_check
      (
      (P_is_vec == false) && (P.n_rows != P.n_cols),
      "diagmat(): only vectors and square matrices are accepted"
      );
    }
  
  inline ~diagmat_proxy_check()
    {
    if(P_local)
      {
      delete P_local;
      }
    }
  
  
  arma_inline elem_type operator[] (const uword i)                    const { return P_is_vec ? P[i] : P.at(i,i);                                         }
  arma_inline elem_type at         (const uword row, const uword col) const { return (row == col) ? ( P_is_vec ? P[row] : P.at(row,row) ) : elem_type(0); }
  

  const Mat<eT>* P_local;
  const Mat<eT>& P;
  const bool     P_is_vec;
  const uword    n_elem;
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
  
  
  arma_inline elem_type operator[] (const uword i)                    const { return P[i];                                 }
  arma_inline elem_type at         (const uword row, const uword col) const { return (row == col) ? P[row] : elem_type(0); }
  
  
  const Row<eT>* P_local;
  const Row<eT>& P;
  const bool     P_is_vec;
  const uword    n_elem;
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
  
  
  arma_inline elem_type operator[] (const uword i)                    const { return P[i];                                 }
  arma_inline elem_type at         (const uword row, const uword col) const { return (row == col) ? P[row] : elem_type(0); }
  
  
  const Col<eT>* P_local;
  const Col<eT>& P;
  const bool     P_is_vec;
  const uword    n_elem;
  };



//! @}
