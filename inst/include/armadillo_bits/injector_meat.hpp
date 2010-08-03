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


//! \addtogroup injector
//! @{



template<typename eT>
inline
mat_injector_row<eT>::mat_injector_row()
  : n_cols(0)
  {
  arma_extra_debug_sigprint();
  
  A.set_size( podarray_prealloc_n_elem::val );
  }



template<typename eT>
inline
void
mat_injector_row<eT>::insert(const eT val) const
  {
  arma_extra_debug_sigprint();
  
  if(n_cols < A.n_elem)
    {
    A[n_cols] = val;
    ++n_cols;
    }
  else
    {
    B.set_size(2 * A.n_elem);
    
    syslib::copy_elem(B.memptr(), A.memptr(), n_cols);
    
    B[n_cols] = val;
    ++n_cols;
    
    std::swap( access::rw(A.mem),    access::rw(B.mem)    );
    std::swap( access::rw(A.n_elem), access::rw(B.n_elem) );
    }
  }



//
//
//



template<typename T1>
inline
mat_injector<T1>::mat_injector(T1& in_X, const typename mat_injector<T1>::elem_type val)
  : X(in_X)
  , n_rows(1)
  {
  arma_extra_debug_sigprint();
  
  typedef typename mat_injector<T1>::elem_type eT;
  
  AA = new podarray< mat_injector_row<eT>* >;
  BB = new podarray< mat_injector_row<eT>* >;
  
  podarray< mat_injector_row<eT>* >& A = *AA;
  
  A.set_size(n_rows);
  
  for(u32 row=0; row<n_rows; ++row)
    {
    A[row] = new mat_injector_row<eT>;
    }
  
  (*(A[0])).insert(val);
  }



template<typename T1>
inline
mat_injector<T1>::mat_injector(T1& in_X, const injector_helper x)
  : X(in_X)
  , n_rows(1)
  {
  arma_extra_debug_sigprint();
  
  typedef typename mat_injector<T1>::elem_type eT;
  
  AA = new podarray< mat_injector_row<eT>* >;
  BB = new podarray< mat_injector_row<eT>* >;
  
  podarray< mat_injector_row<eT>* >& A = *AA;
  
  A.set_size(n_rows);
  
  for(u32 row=0; row<n_rows; ++row)
    {
    A[row] = new mat_injector_row<eT>;
    }
  
  if(x == endr)
    {
    (*this).end_of_row();
    }
  }



template<typename T1>
inline
mat_injector<T1>::~mat_injector()
  {
  arma_extra_debug_sigprint();
  
  typedef typename mat_injector<T1>::elem_type eT;
  
  podarray< mat_injector_row<eT>* >& A = *AA;
  
  if(n_rows > 0)
    {
    u32 max_n_cols = (*(A[0])).n_cols;
    
    for(u32 row=1; row<n_rows; ++row)
      {
      const u32 n_cols = (*(A[row])).n_cols;
      
      if(max_n_cols < n_cols)
        {
        max_n_cols = n_cols;
        }
      }
    
    const u32 max_n_rows = ((*(A[n_rows-1])).n_cols == 0) ? n_rows-1 : n_rows;
    
    if(is_Mat_only<T1>::value == true)
      {
      X.set_size(max_n_rows, max_n_cols);
      
      for(u32 row=0; row<max_n_rows; ++row)
        {
        const u32 n_cols = (*(A[row])).n_cols;
        
        for(u32 col=0; col<n_cols; ++col)
          {
          X.at(row,col) = (*(A[row])).A[col];
          }
        
        for(u32 col=n_cols; col<max_n_cols; ++col)
          {
          X.at(row,col) = eT(0);
          }
        }
      }
    else
    if(is_Row<T1>::value == true)
      {
      arma_debug_check( (max_n_rows > 1), "operator<<: incompatible dimensions" );
      
      const u32 n_cols = (*(A[0])).n_cols;
      
      X.set_size(1, n_cols);
      
      syslib::copy_elem( X.memptr(), (*(A[0])).A.memptr(), n_cols );
      }
    else
    if(is_Col<T1>::value == true)
      {
      const bool is_vec = ( (max_n_rows == 1) || (max_n_cols == 1) );
      
      arma_debug_check( (is_vec == false), "operator<<: incompatible dimensions" );
      
      const u32 n_elem = (std::max)(max_n_rows, max_n_cols);
      
      X.set_size(n_elem, 1);
      
      u32 i = 0;
      for(u32 row=0; row<max_n_rows; ++row)
        {
        const u32 n_cols = (*(A[0])).n_cols;
        
        for(u32 col=0; col<n_cols; ++col)
          {
          X[i] = (*(A[row])).A[col];
          ++i;
          }
        
        for(u32 col=n_cols; col<max_n_cols; ++col)
          {
          X[i] = eT(0);
          ++i;
          }
        }
      }
    }
  
  for(u32 row=0; row<n_rows; ++row)
    {
    delete A[row];
    }
    
  delete AA;
  delete BB;
  }



template<typename T1>
inline
void
mat_injector<T1>::insert(const typename mat_injector<T1>::elem_type val) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename mat_injector<T1>::elem_type eT;
  
  podarray< mat_injector_row<eT>* >& A = *AA;
  
  (*(A[n_rows-1])).insert(val);
  }




template<typename T1>
inline
void
mat_injector<T1>::end_of_row() const
  {
  arma_extra_debug_sigprint();
  
  typedef typename mat_injector<T1>::elem_type eT;
  
  podarray< mat_injector_row<eT>* >& A = *AA;
  podarray< mat_injector_row<eT>* >& B = *BB;
  
  B.set_size( n_rows+1 );
  
  syslib::copy_elem(B.memptr(), A.memptr(), n_rows);
  
  for(u32 row=n_rows; row<(n_rows+1); ++row)
    {
    B[row] = new mat_injector_row<eT>;
    }
  
  std::swap(AA, BB);
  
  n_rows += 1;
  }



template<typename T1>
arma_inline
const mat_injector<T1>&
operator<<(const mat_injector<T1>& ref, const typename mat_injector<T1>::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  ref.insert(val);
  
  return ref;
  }



template<typename T1>
arma_inline
const mat_injector<T1>&
operator<<(const mat_injector<T1>& ref, const injector_helper x)
  {
  arma_extra_debug_sigprint();
  
  if(x == endr)
    {
    ref.end_of_row();
    }
  
  return ref;
  }




//
//
//



template<typename oT>
inline
field_injector_row<oT>::field_injector_row()
  : n_cols(0)
  {
  arma_extra_debug_sigprint();
  
  AA = new field<oT>;
  BB = new field<oT>;
  
  field<oT>& A = *AA;
  
  A.set_size( field_prealloc_n_elem::val );
  }



template<typename oT>
inline
field_injector_row<oT>::~field_injector_row()
  {
  arma_extra_debug_sigprint();
  
  delete AA;
  delete BB;
  }



template<typename oT>
inline
void
field_injector_row<oT>::insert(const oT& val) const
  {
  arma_extra_debug_sigprint();
  
  field<oT>& A = *AA;
  field<oT>& B = *BB;
  
  if(n_cols < A.n_elem)
    {
    A[n_cols] = val;
    ++n_cols;
    }
  else
    {
    B.set_size(2 * A.n_elem);
    
    for(u32 i=0; i<n_cols; ++i)
      {
      B[i] = A[i];
      }
    
    B[n_cols] = val;
    ++n_cols;
    
    std::swap(AA, BB);
    }
  }



//
//
//


template<typename T1>
inline
field_injector<T1>::field_injector(T1& in_X, const typename field_injector<T1>::object_type& val)
  : X(in_X)
  , n_rows(1)
  {
  arma_extra_debug_sigprint();
  
  typedef typename field_injector<T1>::object_type oT;
  
  AA = new podarray< field_injector_row<oT>* >;
  BB = new podarray< field_injector_row<oT>* >;
  
  podarray< field_injector_row<oT>* >& A = *AA;
  
  A.set_size(n_rows);
  
  for(u32 row=0; row<n_rows; ++row)
    {
    A[row] = new field_injector_row<oT>;
    }
  
  (*(A[0])).insert(val);
  }



template<typename T1>
inline
field_injector<T1>::field_injector(T1& in_X, const injector_helper x)
  : X(in_X)
  , n_rows(1)
  {
  arma_extra_debug_sigprint();
  
  typedef typename field_injector<T1>::object_type oT;
  
  AA = new podarray< field_injector_row<oT>* >;
  BB = new podarray< field_injector_row<oT>* >;
  
  podarray< field_injector_row<oT>* >& A = *AA;
  
  A.set_size(n_rows);
  
  for(u32 row=0; row<n_rows; ++row)
    {
    A[row] = new field_injector_row<oT>;
    }
  
  if(x == endr)
    {
    (*this).end_of_row();
    }
  }



template<typename T1>
inline
field_injector<T1>::~field_injector()
  {
  arma_extra_debug_sigprint();
  
  typedef typename field_injector<T1>::object_type oT;
  
  podarray< field_injector_row<oT>* >& A = *AA;
  
  if(n_rows > 0)
    {
    u32 max_n_cols = (*(A[0])).n_cols;
    
    for(u32 row=1; row<n_rows; ++row)
      {
      const u32 n_cols = (*(A[row])).n_cols;
      
      if(max_n_cols < n_cols)
        {
        max_n_cols = n_cols;
        }
      }
      
    const u32 max_n_rows = ((*(A[n_rows-1])).n_cols == 0) ? n_rows-1 : n_rows;
    
    X.set_size(max_n_rows, max_n_cols);
    
    for(u32 row=0; row<max_n_rows; ++row)
      {
      const u32 n_cols = (*(A[row])).n_cols;
      
      for(u32 col=0; col<n_cols; ++col)
        {
        const field<oT>& tmp = *((*(A[row])).AA);
        X.at(row,col) = tmp[col];
        }
      
      for(u32 col=n_cols; col<max_n_cols; ++col)
        {
        X.at(row,col) = oT();
        }
      }
    }
  
  
  for(u32 row=0; row<n_rows; ++row)
    {
    delete A[row];
    }
  
  delete AA;
  delete BB;
  }



template<typename T1>
inline
void
field_injector<T1>::insert(const typename field_injector<T1>::object_type& val) const
  {
  arma_extra_debug_sigprint();
  
  typedef typename field_injector<T1>::object_type oT;
  
  podarray< field_injector_row<oT>* >& A = *AA;
  
  (*(A[n_rows-1])).insert(val);
  }




template<typename T1>
inline
void
field_injector<T1>::end_of_row() const
  {
  arma_extra_debug_sigprint();
  
  typedef typename field_injector<T1>::object_type oT;
  
  podarray< field_injector_row<oT>* >& A = *AA;
  podarray< field_injector_row<oT>* >& B = *BB;
  
  B.set_size( n_rows+1 );
  
  for(u32 row=0; row<n_rows; ++row)
    {
    B[row] = A[row];
    }
  
  for(u32 row=n_rows; row<(n_rows+1); ++row)
    {
    B[row] = new field_injector_row<oT>;
    }
  
  std::swap(AA, BB);
  
  n_rows += 1;
  }



template<typename T1>
arma_inline
const field_injector<T1>&
operator<<(const field_injector<T1>& ref, const typename field_injector<T1>::object_type& val)
  {
  arma_extra_debug_sigprint();
  
  ref.insert(val);
  
  return ref;
  }



template<typename T1>
arma_inline
const field_injector<T1>&
operator<<(const field_injector<T1>& ref, const injector_helper x)
  {
  arma_extra_debug_sigprint();
  
  if(x == endr)
    {
    ref.end_of_row();
    }
  
  return ref;
  }



//! @}
