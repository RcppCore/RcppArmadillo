// Copyright (C) 2009-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2012 Conrad Sanderson
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
      
      out[row] = op_median::direct_median(tmp_vec);
      }
    }
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
        
      out[col] = op_mean::robust_mean(colmem[index1], colmem[index2]);
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
        tmp_vec[col].index = col;
        }
      
      uword index1;
      uword index2;
      op_median::direct_cx_median_index(index1, index2, tmp_vec);
      
      out[row] = op_mean::robust_mean( X.at(row,index1), X.at(row,index2) );
      }
    }
  }



template<typename T1>
inline
typename T1::elem_type
op_median::median_vec
  (
  const T1& X,
  const typename arma_not_cx<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;

  const Proxy<T1> P(X);
  
  const uword n_elem = P.get_n_elem();
  
  arma_debug_check( (n_elem == 0), "median(): given object has no elements" );
  
  std::vector<eT> tmp_vec(n_elem);
    
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    ea_type A = P.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      tmp_vec[i] = A[i];
      }
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_cols == 1)
      {
      for(uword row=0; row < n_rows; ++row)
        {
        tmp_vec[row] = P.at(row,0);
        }
      }
    else
    if(n_rows == 1)
      {
      for(uword col=0; col < n_cols; ++col)
        {
        tmp_vec[col] = P.at(0,col);
        }
      }
    else
      {
      arma_stop("op_median::median_vec(): expected a vector" );
      }
    }
  
  return op_median::direct_median(tmp_vec);
  }



template<typename T1>
inline
typename T1::elem_type
op_median::median_vec
  (
  const T1& X,
  const typename arma_cx_only<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  const Proxy<T1> P(X);
  
  const uword n_elem = P.get_n_elem();
  
  arma_debug_check( (n_elem == 0), "median(): given object has no elements" );
  
  std::vector< arma_cx_median_packet<T> > tmp_vec(n_elem);
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typedef typename Proxy<T1>::ea_type ea_type;
    
    ea_type A = P.get_ea();
    
    for(uword i=0; i<n_elem; ++i)
      {
      tmp_vec[i].val   = std::abs( A[i] );
      tmp_vec[i].index = i;
      }
    
    uword index1;
    uword index2;
    op_median::direct_cx_median_index(index1, index2, tmp_vec);
    
    return op_mean::robust_mean( A[index1], A[index2] );
    }
  else
    {
    const uword n_rows = P.get_n_rows();
    const uword n_cols = P.get_n_cols();
    
    if(n_cols == 1)
      {
      for(uword row=0; row < n_rows; ++row)
        {
        tmp_vec[row].val   = std::abs( P.at(row,0) );
        tmp_vec[row].index = row;
        }
      
      uword index1;
      uword index2;
      op_median::direct_cx_median_index(index1, index2, tmp_vec);
      
      return op_mean::robust_mean( P.at(index1,0), P.at(index2,0) );
      }
    else
    if(n_rows == 1)
      {
      for(uword col=0; col < n_cols; ++col)
        {
        tmp_vec[col].val   = std::abs( P.at(0,col) );
        tmp_vec[col].index = col;
        }
      
      uword index1;
      uword index2;
      op_median::direct_cx_median_index(index1, index2, tmp_vec);
      
      return op_mean::robust_mean( P.at(0,index1), P.at(0,index2) );
      }
    else
      {
      arma_stop("op_median::median_vec(): expected a vector" );
      
      return eT(0);
      }
    }
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
    return op_mean::robust_mean(X[half-1], X[half]);
    }
  else
    {
    return X[half];
    }
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



//! @}

