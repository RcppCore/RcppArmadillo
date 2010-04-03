// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_pinv
//! @{



template<typename eT>
inline
void
op_pinv::direct_pinv(Mat<eT>& out, const Mat<eT>& A, eT tol)
  {
  arma_extra_debug_sigprint();
  
  const u32 n_rows = A.n_rows;
  const u32 n_cols = A.n_cols;
  
  // SVD decomposition 
  Mat<eT> U;
  Col<eT> s;
  Mat<eT> V;
  
  const bool status = (n_cols > n_rows) ? svd(U,s,V,trans(A)) : svd(U,s,V,A);
  
  if(status == false)
    {
    out.set_size(0,0);
    return;
    }
  
  // set tolerance to default if it hasn't been specified as an argument
  if(tol == eT(0))
    {
    tol = (std::max)(n_rows,n_cols) * eop_aux::direct_eps(max(s));
    }
   
  // count non zero valued elements in s
  u32 count = 0; 
  for(u32 i = 0; i < s.n_elem; ++i)
    {
    if(s[i] > tol)
      {
      ++count;
      }
    }
  
  
  if(count != 0)
    {
    // reduce the length of s in order to contain only the values above tolerance
    s = s.rows(0,count-1);
    
    // set the elements of s equal to their reciprocals
    s = eT(1) / s;
    
    if(A.n_cols <= A.n_rows)
      {
      out = V.cols(0,count-1) * diagmat(s) * trans(U.cols(0,count-1));
      }
    else
      {
      out = U.cols(0,count-1) * diagmat(s) * trans(V.cols(0,count-1));
      }
    }
  else
    {
    out.zeros(n_cols, n_rows);
    }
  }



template<typename T>
inline
void
op_pinv::direct_pinv(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& A, T tol)
  {
  arma_extra_debug_sigprint();
  
  const u32 n_rows = A.n_rows;
  const u32 n_cols = A.n_cols;
  
  typedef typename std::complex<T> eT;
 
  // SVD decomposition 
  Mat<eT> U;
  Col< T> s;
  Mat<eT> V;
  
  const bool status = (n_cols > n_rows) ? svd(U,s,V,htrans(A)) : svd(U,s,V,A);
  
  if(status == false)
    {
    out.set_size(0,0);
    return;
    }
 
  // set tolerance to default if it hasn't been specified as an argument 
  if(tol == T(0))
    {
    tol = (std::max)(n_rows,n_cols) * eop_aux::direct_eps(max(s));
    }
  
  
  // count non zero valued elements in s
  u32 count = 0;
  for(u32 i = 0; i < s.n_elem; ++i)
    {
    if(s[i] > tol)
      {
      ++count;
      }
    }
  
  if(count != 0)
    {
    // reduce the length of s in order to contain only the values above tolerance
    s = s.rows(0,count-1);

    // set the elements of s equal to their reciprocals
    s = T(1) / s;
    
    if(n_rows >= n_cols)
      {
      out = V.cols(0,count-1) * diagmat(s) * htrans(U.cols(0,count-1));
      }
    else
      {
      out = U.cols(0,count-1) * diagmat(s) * htrans(V.cols(0,count-1));
      }
    }
  else
    {
    out.zeros(n_cols, n_rows);
    }
  }



template<typename T1>
inline
void
op_pinv::apply(Mat<typename T1::pod_type>& out, const Op<T1,op_pinv>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type eT;
  
  const eT tol = in.aux; 
  
  arma_debug_check((tol < eT(0)), "pinv(): tol must be >= 0");
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& A     = tmp.M;
  
  op_pinv::direct_pinv(out, A, tol);
  }



template<typename T1>
inline
void
op_pinv::apply(Mat< std::complex<typename T1::pod_type> >& out, const Op<T1,op_pinv>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  typedef typename std::complex<typename T1::pod_type> eT;
  
  const T tol = in.aux.real();
  
  arma_debug_check((tol < T(0)), "pinv(): tol must be >= 0");
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& A     = tmp.M;
  
  op_pinv::direct_pinv(out, A, tol);
  }



//! @}
