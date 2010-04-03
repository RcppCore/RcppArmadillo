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


//! \addtogroup op_sort
//! @{


// using qsort() rather than std::sort() for now.
// std::sort() will be used when a Random Access Iterator wrapper for plain arrays is ready,
// otherwise using std::sort() would currently entail copying elements to/from std::vector

template<typename eT>
class arma_qsort_helper
  {
  public:
  
  static
  int
  ascend_compare(const void* A_orig, const void* B_orig)
    {
    const eT& A = *(static_cast<const eT*>(A_orig));
    const eT& B = *(static_cast<const eT*>(B_orig));
    
    if(A < B)
      {
      return -1;
      }
    else
    if(A > B)
      {
      return +1;
      }
    else
      {
      return 0;
      }
    }
  
  
  
  static
  int
  descend_compare(const void* A_orig, const void* B_orig)
    {
    const eT& A = *(static_cast<const eT*>(A_orig));
    const eT& B = *(static_cast<const eT*>(B_orig));
    
    if(A < B)
      {
      return +1;
      }
    else
    if(A > B)
      {
      return -1;
      }
    else
      {
      return 0;
      }
    }
  
  
  };



//template<>
template<typename T>
class arma_qsort_helper< std::complex<T> >
  {
  public:
  
  typedef typename std::complex<T> eT;
  
  
  static
  int
  ascend_compare(const void* A_orig, const void* B_orig)
    {
    const eT& A = *(static_cast<const eT*>(A_orig));
    const eT& B = *(static_cast<const eT*>(B_orig));
    
    const T abs_A = std::abs(A);
    const T abs_B = std::abs(B);
    
    if(abs_A < abs_B)
      {
      return -1;
      }
    else
    if(abs_A > abs_B)
      {
      return +1;
      }
    else
      {
      return 0;
      }
    }
  
  
  
  static
  int
  descend_compare(const void* A_orig, const void* B_orig)
    {
    const eT& A = *(static_cast<const eT*>(A_orig));
    const eT& B = *(static_cast<const eT*>(B_orig));
    
    const T abs_A = std::abs(A);
    const T abs_B = std::abs(B);
    
    if(abs_A < abs_B)
      {
      return +1;
      }
    else
    if(abs_A > abs_B)
      {
      return -1;
      }
    else
      {
      return 0;
      }
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
    std::qsort(X, n_elem, sizeof(eT), arma_qsort_helper<eT>::ascend_compare);
    }
  else
    {
    std::qsort(X, n_elem, sizeof(eT), arma_qsort_helper<eT>::descend_compare);
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
  
  arma_debug_check( (X.is_finite() == false), "sort(): given object has non-finite elements"     );
  arma_debug_check( (sort_type > 1),          "sort(): incorrect usage. sort_type must be 0 or 1");
  arma_debug_check( (dim > 1),                "sort(): incorrect usage. dim must be 0 or 1"      );
  
  
  if(dim == 0)  // column-wise
    {
    arma_extra_debug_print("op_sort::apply(), dim = 0");
    
    out = X;
    
    for(u32 col=0; col<out.n_cols; ++col)
      {
      op_sort::direct_sort( out.colptr(col), out.n_rows, sort_type );
      }
    }
  else
  if(dim == 1)  // row-wise
    {
    if(X.n_rows != 1)  // not a row vector
      {
      arma_extra_debug_print("op_sort::apply(), dim = 1, generic");
      
      //out.set_size(X.n_rows, X.n_cols);
      out.copy_size(X);
      
      podarray<eT> tmp_array(X.n_cols);
      
      for(u32 row=0; row<out.n_rows; ++row)
        {
        
        for(u32 col=0; col<out.n_cols; ++col)
          {
          tmp_array[col] = X.at(row,col);
          }
        
        op_sort::direct_sort( tmp_array.memptr(), out.n_cols, sort_type );
        
        for(u32 col=0; col<out.n_cols; ++col)
          {
          out.at(row,col) = tmp_array[col];
          }
        
        }
      }
    else  // a row vector
      {
      arma_extra_debug_print("op_sort::apply(), dim = 1, vector specific");
      
      out = X;
      op_sort::direct_sort(out.memptr(), out.n_elem, sort_type);
      }
    }
  
  }


//! @}
