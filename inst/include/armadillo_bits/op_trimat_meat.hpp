// Copyright (C) 2010-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2012 Conrad Sanderson
// Copyright (C) 2011      Ryan Curtin
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_trimat
//! @{



template<typename eT>
inline
void
op_trimat::fill_zeros(Mat<eT>& out, const bool upper)
  {
  arma_extra_debug_sigprint();
  
  const uword N = out.n_rows;
  
  if(upper)
    {
    // upper triangular: set all elements below the diagonal to zero
    
    for(uword i=0; i<N; ++i)
      {
      eT* data = out.colptr(i);
      
      arrayops::inplace_set( &data[i+1], eT(0), (N-(i+1)) );
      }
    }
  else
    {
    // lower triangular: set all elements above the diagonal to zero
    
    for(uword i=1; i<N; ++i)
      {
      eT* data = out.colptr(i);
      
      arrayops::inplace_set( data, eT(0), i );
      }
    }
  }



template<typename T1>
inline
void
op_trimat::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trimat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.is_square() == false), "trimatu()/trimatl(): given matrix must be square" );
  
  const uword N     = A.n_rows;
  const bool  upper = (in.aux_uword_a == 0);
  
  if(&out != &A)
    {
    out.copy_size(A);
    
    if(upper)
      {
      // upper triangular: copy the diagonal and the elements above the diagonal
      for(uword i=0; i<N; ++i)
        {
        const eT* A_data   = A.colptr(i);
              eT* out_data = out.colptr(i);
        
        arrayops::copy( out_data, A_data, i+1 );
        }
      }
    else
      {
      // lower triangular: copy the diagonal and the elements below the diagonal
      for(uword i=0; i<N; ++i)
        {
        const eT* A_data   = A.colptr(i);
              eT* out_data = out.colptr(i);
        
        arrayops::copy( &out_data[i], &A_data[i], N-i );
        }
      }
    }
  
  op_trimat::fill_zeros(out, upper);
  }



template<typename T1>
inline
void
op_trimat::apply(Mat<typename T1::elem_type>& out, const Op<Op<T1, op_htrans>, op_trimat>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m.m);
  const Mat<eT>& A = tmp.M;
  
  const bool upper = (in.aux_uword_a == 0);
  
  op_trimat::apply_htrans(out, A, upper);
  }



template<typename eT>
inline
void
op_trimat::apply_htrans
  (
        Mat<eT>& out,
  const Mat<eT>& A,
  const bool     upper,
  const typename arma_not_cx<eT>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  // This specialisation is for trimatl(trans(X)) = trans(trimatu(X)) and also
  // trimatu(trans(X)) = trans(trimatl(X)).  We want to avoid the creation of an
  // extra temporary.
  
  // It doesn't matter if the input and output matrices are the same; we will
  // pull data from the upper or lower triangular to the lower or upper
  // triangular (respectively) and then set the rest to 0, so overwriting issues
  // aren't present.
  
  arma_debug_check( (A.is_square() == false), "trimatu()/trimatl(): given matrix must be square" );
  
  const uword N = A.n_rows;
  
  if(&out != &A)
    {
    out.copy_size(A);
    }
  
  // We can't really get away with any array copy operations here,
  // unfortunately...
  
  if(upper)
    {
    // Upper triangular: but since we're transposing, we're taking the lower
    // triangular and putting it in the upper half.
    for(uword row = 0; row < N; ++row)
      {
      eT* out_colptr = out.colptr(row);
      
      for(uword col = 0; col <= row; ++col)
        {
        //out.at(col, row) = A.at(row, col);
        out_colptr[col] = A.at(row, col);
        }
      }
    }
  else
    {
    // Lower triangular: but since we're transposing, we're taking the upper
    // triangular and putting it in the lower half.
    for(uword row = 0; row < N; ++row)
      {
      for(uword col = row; col < N; ++col)
        {
        out.at(col, row) = A.at(row, col);
        }
      }
    }
  
  op_trimat::fill_zeros(out, upper);
  }



template<typename eT>
inline
void
op_trimat::apply_htrans
  (
        Mat<eT>& out,
  const Mat<eT>& A,
  const bool     upper,
  const typename arma_cx_only<eT>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check( (A.is_square() == false), "trimatu()/trimatl(): given matrix must be square" );
  
  const uword N = A.n_rows;
  
  if(&out != &A)
    {
    out.copy_size(A);
    }
  
  if(upper)
    {
    // Upper triangular: but since we're transposing, we're taking the lower
    // triangular and putting it in the upper half.
    for(uword row = 0; row < N; ++row)
      {
      eT* out_colptr = out.colptr(row);
      
      for(uword col = 0; col <= row; ++col)
        {
        //out.at(col, row) = std::conj( A.at(row, col) );
        out_colptr[col] = std::conj( A.at(row, col) );
        }
      }
    }
  else
    {
    // Lower triangular: but since we're transposing, we're taking the upper
    // triangular and putting it in the lower half.
    for(uword row = 0; row < N; ++row)
      {
      for(uword col = row; col < N; ++col)
        {
        out.at(col, row) = std::conj( A.at(row, col) );
        }
      }
    }
  
  op_trimat::fill_zeros(out, upper);
  }



//! @}
