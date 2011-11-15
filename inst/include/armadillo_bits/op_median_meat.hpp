// Copyright (C) 2009-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2011 Conrad Sanderson
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



template<typename eT>
arma_inline
eT
op_median::robust_mean(const eT A, const eT B)
  {
  return A + (B - A)/eT(2);
  }



//! find the median value of a std::vector (contents is modified)
template<typename eT>
inline 
eT
op_median::direct_median(std::vector<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  const uword n_elem = X.size();
  const uword half   = n_elem/2;
  
  std::sort(X.begin(), X.end());
  
  if((n_elem % 2) == 0)
    {
    return op_median::robust_mean(X[half-1], X[half]);
    }
  else
    {
    return X[half];
    }
  }



template<typename eT>
inline 
eT
op_median::direct_median(const eT* X, const uword n_elem)
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
  
  const uword X_n_elem = X.n_elem;
  
  std::vector<eT> tmp(X_n_elem);
  
  for(uword i=0; i<X_n_elem; ++i)
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
  
  const uword X_n_elem = X.n_elem;
  
  std::vector<eT> tmp(X_n_elem);
  
  for(uword i=0; i<X_n_elem; ++i)
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
  const Mat<eT>&     X = tmp.M;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword dim = in.aux_uword_a;
  arma_debug_check( (dim > 1), "median(): incorrect usage. dim must be 0 or 1");
  
  if(dim == 0)  // in each column
    {
    arma_extra_debug_print("op_median::apply(), dim = 0");
    
    arma_debug_check( (X_n_rows == 0), "median(): given object has zero rows" );

    out.set_size(1, X_n_cols);
    
    std::vector<eT> tmp_vec(X_n_rows);
      
    for(uword col=0; col<X_n_cols; ++col)
      {
      const eT* colmem = X.colptr(col);
      
      for(uword row=0; row<X_n_rows; ++row)
        {
        tmp_vec[row] = colmem[row];
        }
      
      out[col] = op_median::direct_median(tmp_vec);
      }
    }
  else
  if(dim == 1)  // in each row
    {
    arma_extra_debug_print("op_median::apply(), dim = 1");
    
    arma_debug_check( (X_n_cols == 0), "median(): given object has zero columns" );

    out.set_size(X_n_rows, 1);
    
    std::vector<eT> tmp_vec(X_n_cols);
      
    for(uword row=0; row<X_n_rows; ++row)
      {
      for(uword col=0; col<X_n_cols; ++col)
        {
        tmp_vec[col] = X.at(row,col);
        }
      
      out[row] =  op_median::direct_median(tmp_vec);
      }
    }
  }



template<typename T>
arma_inline
std::complex<T>
op_median::robust_mean(const std::complex<T>& A, const std::complex<T>& B)
  {
  return A + (B - A)/T(2);
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index
  (
  uword& out_index1, 
  uword& out_index2, 
  std::vector< arma_cx_median_packet<T> >& X
  )
  {
  arma_extra_debug_sigprint();
  
  const uword n_elem = X.size();
  const uword half   = n_elem/2;
  
  std::sort(X.begin(), X.end());
  
  if((n_elem % 2) == 0)
    {
    out_index1 = X[half-1].index;
    out_index2 = X[half  ].index;
    }
  else
    {
    out_index1 = X[half].index;
    out_index2 = out_index1;
    }
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index
  (
  uword& out_index1, 
  uword& out_index2, 
  const std::complex<T>* X, 
  const uword n_elem
  )
  {
  arma_extra_debug_sigprint();
  
  std::vector< arma_cx_median_packet<T> > tmp(n_elem);
  
  for(uword i=0; i<n_elem; ++i)
    {
    tmp[i].val   = std::abs(X[i]);
    tmp[i].index = i;
    }
  
  op_median::direct_cx_median_index(out_index1, out_index2, tmp);
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index
  (
  uword& out_index1, 
  uword& out_index2, 
  const subview< std::complex<T> >&X
  )
  {
  arma_extra_debug_sigprint();
  
  const uword n_elem = X.n_elem;
  
  std::vector< arma_cx_median_packet<T> > tmp(n_elem);
  
  for(uword i=0; i<n_elem; ++i)
    {
    tmp[i].val   = std::abs(X[i]);
    tmp[i].index = i;
    }
  
  op_median::direct_cx_median_index(out_index1, out_index2, tmp);
  }



template<typename T>
inline 
void
op_median::direct_cx_median_index
  (
  uword& out_index1, 
  uword& out_index2, 
  const diagview< std::complex<T> >&X
  )
  {
  arma_extra_debug_sigprint();
  
  const uword n_elem = X.n_elem;
  
  std::vector< arma_cx_median_packet<T> > tmp(n_elem);
  
  for(uword i=0; i<n_elem; ++i)
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
  
  arma_type_check(( is_same_type<eT, typename T1::elem_type>::value == false ));
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>&     X = tmp.M;
  
  const uword X_n_rows = X.n_rows;
  const uword X_n_cols = X.n_cols;
  
  const uword dim = in.aux_uword_a;
  arma_debug_check( (dim > 1), "median(): incorrect usage. dim must be 0 or 1");
  
  if(dim == 0)  // in each column
    {
    arma_extra_debug_print("op_median::apply(), dim = 0");
    
    arma_debug_check( (X_n_rows == 0), "median(): given object has zero rows" );

    out.set_size(1, X_n_cols);
    
    std::vector< arma_cx_median_packet<T> > tmp_vec(X_n_rows);
    
    for(uword col=0; col<X_n_cols; ++col)
      {
      const eT* colmem = X.colptr(col);
      
      for(uword row=0; row<X_n_rows; ++row)
        {
        tmp_vec[row].val   = std::abs(colmem[row]);
        tmp_vec[row].index = row;
        }
      
      uword index1;
      uword index2;
      op_median::direct_cx_median_index(index1, index2, tmp_vec);
        
      out[col] = op_median::robust_mean(colmem[index1], colmem[index2]);
      }
    }
  else
  if(dim == 1)  // in each row
    {
    arma_extra_debug_print("op_median::apply(), dim = 1");
    
    arma_debug_check( (X_n_cols == 0), "median(): given object has zero columns" );

    out.set_size(X_n_rows, 1);
    
    std::vector< arma_cx_median_packet<T> > tmp_vec(X_n_cols);
    
    for(uword row=0; row<X_n_rows; ++row)
      {
      for(uword col=0; col<X_n_cols; ++col)
        {
        tmp_vec[col].val   = std::abs(X.at(row,col));
        tmp_vec[row].index = col;
        }
      
      uword index1;
      uword index2;
      op_median::direct_cx_median_index(index1, index2, tmp_vec);
      
      out[row] = op_median::robust_mean( X.at(row,index1), X.at(row,index2) );
      }
    }
  }



//! @}

