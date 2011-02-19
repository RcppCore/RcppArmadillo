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


//! \addtogroup op_sort
//! @{



template<typename eT>
class arma_ascend_sort_helper
  {
  public:
  
  arma_inline
  bool
  operator() (eT a, eT b) const
    {
    return (a < b);
    }
  };
  


template<typename eT>
class arma_descend_sort_helper
  {
  public:
  
  arma_inline
  bool
  operator() (eT a, eT b) const
    {
    return (a > b);
    }
  };
  


template<typename T>
class arma_ascend_sort_helper< std::complex<T> >
  {
  public:
  
  typedef typename std::complex<T> eT;
  
  inline
  bool
  operator() (const eT& a, const eT& b) const
    {
    return (std::abs(a) < std::abs(b));
    }
  };



template<typename T>
class arma_descend_sort_helper< std::complex<T> >
  {
  public:
  
  typedef typename std::complex<T> eT;
  
  inline
  bool
  operator() (const eT& a, const eT& b) const
    {
    return (std::abs(a) > std::abs(b));
    }
  };



template<typename eT>
inline 
void
op_sort::direct_sort(eT* X, const u32 n_elem, const u32 sort_type)
  {
  arma_extra_debug_sigprint();
  
  if(sort_type == 0)
    {
    arma_ascend_sort_helper<eT> comparator;
    
    std::sort(&X[0], &X[n_elem], comparator);
    }
  else
    {
    arma_descend_sort_helper<eT> comparator;
    
    std::sort(&X[0], &X[n_elem], comparator);
    }
  }



template<typename eT>
inline 
void
op_sort::copy_row(eT* X, const Mat<eT>& A, const u32 row)
  {
  const u32 N = A.n_cols;
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    X[i] = A.at(row,i);
    X[j] = A.at(row,j);
    }
  
  if(i < N)
    {
    X[i] = A.at(row,i);
    }
  }



template<typename eT>
inline 
void
op_sort::copy_row(Mat<eT>& A, const eT* X, const u32 row)
  {
  const u32 N = A.n_cols;
  
  u32 i,j;
  
  for(i=0, j=1; j<N; i+=2, j+=2)
    {
    A.at(row,i) = X[i]; 
    A.at(row,j) = X[j];
    }
  
  if(i < N)
    {
    A.at(row,i) = X[i];
    }
  }



template<typename T1>
inline
void
op_sort::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sort>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  const u32 sort_type = in.aux_u32_a;
  const u32 dim       = in.aux_u32_b;
  
  arma_debug_check( (sort_type > 1),          "sort(): incorrect usage. sort_type must be 0 or 1");
  arma_debug_check( (dim > 1),                "sort(): incorrect usage. dim must be 0 or 1"      );
  arma_debug_check( (X.is_finite() == false), "sort(): given object has non-finite elements"     );
  
  if( (X.n_rows * X.n_cols) <= 1 )
    {
    out = X;
    return;
    }
  
  
  if(dim == 0)  // sort the contents of each column
    {
    arma_extra_debug_print("op_sort::apply(), dim = 0");
    
    out = X;
    
    const u32 n_rows = out.n_rows;
    const u32 n_cols = out.n_cols;
      
    for(u32 col=0; col < n_cols; ++col)
      {
      op_sort::direct_sort( out.colptr(col), n_rows, sort_type );
      }
    }
  else
  if(dim == 1)  // sort the contents of each row
    {
    if(X.n_rows == 1)  // a row vector
      {
      arma_extra_debug_print("op_sort::apply(), dim = 1, vector specific");
      
      out = X;
      op_sort::direct_sort(out.memptr(), out.n_elem, sort_type);
      }
    else  // not a row vector
      {
      arma_extra_debug_print("op_sort::apply(), dim = 1, generic");
      
      out.copy_size(X);
      
      const u32 n_rows = out.n_rows;
      const u32 n_cols = out.n_cols;
      
      podarray<eT> tmp_array(n_cols);
      
      for(u32 row=0; row < n_rows; ++row)
        {
        op_sort::copy_row(tmp_array.memptr(), X, row);
        
        op_sort::direct_sort( tmp_array.memptr(), n_cols, sort_type );
        
        op_sort::copy_row(out, tmp_array.memptr(), row);
        }
      }
    }
  
  }


//! @}
