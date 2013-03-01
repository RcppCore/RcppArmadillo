// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup glue_join
//! @{



template<typename T1, typename T2>
inline
void
glue_join::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_join>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;

  const unwrap<T1> A_tmp(X.A);
  const unwrap<T2> B_tmp(X.B);
  
  const Mat<eT>& A = A_tmp.M;
  const Mat<eT>& B = B_tmp.M;
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword B_n_rows = B.n_rows;
  const uword B_n_cols = B.n_cols;
  
  const uword join_type = X.aux_uword;
  
  if(join_type == 0)
    {
    arma_debug_check
      (
      ( (A_n_cols != B_n_cols) && ( (A_n_rows > 0) || (A_n_cols > 0) ) && ( (B_n_rows > 0) || (B_n_cols > 0) ) ),
      "join_cols(): number of columns must be the same"
      );
    }
  else
    {
    arma_debug_check
      (
      ( (A_n_rows != B.n_rows) && ( (A_n_rows > 0) || (A_n_cols > 0) ) && ( (B_n_rows > 0) || (B_n_cols > 0) ) ),
      "join_rows(): number of rows must be the same"
      );
    }
  
  
  if( (&out != &A) && (&out != &B) )
    {
    if(join_type == 0)   // join columns (i.e. result matrix has more rows)
      {
      out.set_size( A_n_rows + B_n_rows, (std::max)(A_n_cols, B_n_cols) );
      
      if( out.n_elem > 0 )
        {
        if(A.is_empty() == false)
          { 
          out.submat(0,        0,   A_n_rows-1, out.n_cols-1) = A;
          }
          
        if(B.is_empty() == false)
          {
          out.submat(A_n_rows, 0, out.n_rows-1, out.n_cols-1) = B;
          }
        }
      }
    else   // join rows  (i.e. result matrix has more columns)
      {
      out.set_size( (std::max)(A_n_rows, B_n_rows), A_n_cols + B_n_cols );
      
      if( out.n_elem > 0 )
        {
        if(A.is_empty() == false)
          {
          out.submat(0, 0,        out.n_rows-1,   A.n_cols-1) = A;
          }
        
        if(B.is_empty() == false)
          {
          out.submat(0, A_n_cols, out.n_rows-1, out.n_cols-1) = B;
          }
        }
      }
    }
  else  // we have aliasing
    {
    Mat<eT> C;
    
    if(join_type == 0)
      {
      C.set_size( A_n_rows + B_n_rows, (std::max)(A_n_cols, B_n_cols) );
      
      if( C.n_elem > 0 )
        {
        if(A.is_empty() == false)
          {
          C.submat(0,        0, A_n_rows-1, C.n_cols-1) = A;
          }
        
        if(B.is_empty() == false)
          {
          C.submat(A_n_rows, 0, C.n_rows-1, C.n_cols-1) = B;
          }
        }
      }
    else
      {
      C.set_size( (std::max)(A_n_rows, B_n_rows), A_n_cols + B_n_cols );
      
      if( C.n_elem > 0 )
        {
        if(A.is_empty() == false)
          {
          C.submat(0, 0,        C.n_rows-1, A_n_cols-1) = A;
          }
        
        if(B.is_empty() == false)
          {
          C.submat(0, A_n_cols, C.n_rows-1, C.n_cols-1) = B;
          }
        }
      }
    
    out.steal_mem(C);
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_join::apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_join>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;

  const unwrap_cube<T1> A_tmp(X.A);
  const unwrap_cube<T2> B_tmp(X.B);
  
  const Cube<eT>& A = A_tmp.M;
  const Cube<eT>& B = B_tmp.M;
  
  if(A.n_elem == 0)
    {
    out = B;
    return;
    }
  
  if(B.n_elem == 0)
    {
    out = A;
    return;
    }
  
  
  arma_debug_check( ( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) ), "join_slices(): size of slices must be the same" );
  
  
  if( (&out != &A) && (&out != &B) )
    {
    out.set_size(A.n_rows, A.n_cols, A.n_slices + B.n_slices);
    
    out.slices(0,          A.n_slices-1  ) = A;
    out.slices(A.n_slices, out.n_slices-1) = B;
    }
  else  // we have aliasing
    {
    Cube<eT> C(A.n_rows, A.n_cols, A.n_slices + B.n_slices);
    
    C.slices(0,          A.n_slices-1) = A;
    C.slices(A.n_slices, C.n_slices-1) = B;
    
    out.steal_mem(C);
    }
  
  }



//! @}
