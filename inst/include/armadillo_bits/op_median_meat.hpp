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


//! \addtogroup op_median
//! @{



//! find the median value of a std::vector (contents is modified)
template<typename eT>
inline 
eT
op_median::direct_median(std::vector<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  std::sort(X.begin(), X.end());
  
  const u32 n_elem = X.size();
  const u32 half   = n_elem/2;
  
  if((n_elem % 2) == 0)
    {
    return (X[half-1] + X[half]) / eT(2);
    }
  else
    {
    return X[half];
    }
  }



template<typename eT>
inline 
eT
op_median::direct_median(const eT* X, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  std::vector<eT> tmp(X, X+n_elem);
  return op_median::direct_median(tmp);
  }



template<typename eT>
inline 
eT
op_median::direct_median(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  std::vector<eT> tmp(X.n_elem);
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    tmp[i] = X[i];
    }
  
  return op_median::direct_median(tmp);
  }



template<typename eT>
inline 
eT
op_median::direct_median(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  std::vector<eT> tmp(X.n_elem);
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    tmp[i] = X[i];
    }
  
  return op_median::direct_median(tmp);
  }



//! \brief
//! For each row or for each column, find the median value.
//! The result is stored in a dense matrix that has either one column or one row.
//! The dimension, for which the medians are found, is set via the median() function.
template<typename T1>
inline
void
op_median::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_median>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem == 0), "median(): given matrix has no elements" );
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "median(): incorrect usage. dim must be 0 or 1");
  
  
  if(dim == 0)  // column-wise
    {
    arma_extra_debug_print("op_median::apply(), dim = 0");
    
    out.set_size(1, X.n_cols);
    
    std::vector<eT> tmp_vec(X.n_rows);
    
    for(u32 col=0; col<X.n_cols; ++col)
      {
      const eT* colmem = X.colptr(col);
      
      for(u32 row=0; row<X.n_rows; ++row)
        {
        tmp_vec[row] = colmem[row];
        }
      
      out[col] = op_median::direct_median(tmp_vec);
      }
    }
  else
  if(dim == 1)  // row-wise
    {
    arma_extra_debug_print("op_median::apply(), dim = 1");
  
    out.set_size(X.n_rows, 1);
    
    std::vector<eT> tmp_vec(X.n_cols);
    
    for(u32 row=0; row<X.n_rows; ++row)
      {
      for(u32 col=0; col<X.n_cols; ++col)
        {
        tmp_vec[col] = X.at(row,col);
        }
  
      out[row] = op_median::direct_median(tmp_vec);
      }
    }
  
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index(u32& out_index1, u32& out_index2, std::vector< arma_cx_median_packet<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  std::sort(X.begin(), X.end());
  
  const u32 n_elem = X.size();
  const u32 half   = n_elem/2;
  
  if((n_elem % 2) == 0)
    {
    out_index1 = X[half-1].index;
    out_index2 = X[half].index;
    }
  else
    {
    out_index1 = X[half].index;
    out_index2 = X[half].index;
    }
  
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index(u32& out_index1, u32& out_index2, const std::complex<T>* X, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  std::vector< arma_cx_median_packet<T> > tmp(n_elem);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    tmp[i].val   = std::abs(X[i]);
    tmp[i].index = i;
    }
  
  op_median::direct_cx_median_index(out_index1, out_index2, tmp);
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index(u32& out_index1, u32& out_index2, const subview< std::complex<T> >&X)
  {
  arma_extra_debug_sigprint();
  
  const u32 n_elem = X.n_elem;
  
  std::vector< arma_cx_median_packet<T> > tmp(n_elem);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    tmp[i].val   = std::abs(X[i]);
    tmp[i].index = i;
    }
  
  op_median::direct_cx_median_index(out_index1, out_index2, tmp);
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index(u32& out_index1, u32& out_index2, const diagview< std::complex<T> >&X)
  {
  arma_extra_debug_sigprint();
  
  const u32 n_elem = X.n_elem;
  
  std::vector< arma_cx_median_packet<T> > tmp(n_elem);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    tmp[i].val   = std::abs(X[i]);
    tmp[i].index = i;
    }
  
  op_median::direct_cx_median_index(out_index1, out_index2, tmp);
  }



//! Implementation for complex numbers
template<typename T, typename T1>
inline
void
op_median::apply(Mat< std::complex<T> >& out, const Op<T1,op_median>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& X = tmp.M;
  
  arma_debug_check( (X.n_elem == 0), "median(): given matrix has no elements" );
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "median(): incorrect usage. dim must be 0 or 1");
  
  
  if(dim == 0)  // column-wise
    {
    arma_extra_debug_print("op_median::apply(), dim = 0");
    
    out.set_size(1, X.n_cols);
    
    std::vector< arma_cx_median_packet<T> > tmp_vec(X.n_rows);
    
    for(u32 col=0; col<X.n_cols; ++col)
      {
      const eT* colmem = X.colptr(col);
      
      for(u32 row=0; row<X.n_rows; ++row)
        {
        tmp_vec[row].val   = std::abs(colmem[row]);
        tmp_vec[row].index = row;
        }
      
      u32 index1;
      u32 index2;
      op_median::direct_cx_median_index(index1, index2, tmp_vec);
      
      out[col] = (colmem[index1] + colmem[index2]) / T(2);
      }
    }
  else
  if(dim == 1)  // row-wise
    {
    arma_extra_debug_print("op_median::apply(), dim = 1");
  
    out.set_size(X.n_rows, 1);
    
    std::vector< arma_cx_median_packet<T> > tmp_vec(X.n_cols);
    
    for(u32 row=0; row<X.n_rows; ++row)
      {
      for(u32 col=0; col<X.n_cols; ++col)
        {
        tmp_vec[col].val   = std::abs(X.at(row,col));
        tmp_vec[row].index = col;
        }
  
      u32 index1;
      u32 index2;
      op_median::direct_cx_median_index(index1, index2, tmp_vec);
      
      out[row] = ( X.at(row,index1) + X.at(row,index2) ) / T(2);
      }
    }
  
  }



//! @}
