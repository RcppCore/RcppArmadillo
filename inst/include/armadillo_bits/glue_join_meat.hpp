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
  
  const u32 join_type = X.aux_u32;
  
  
  if(join_type == 0)
    {
    arma_debug_check( (A.n_cols != B.n_cols), "join_cols(): number of columns must be the same" );
    }
  else
    {
    arma_debug_check( (A.n_rows != B.n_rows), "join_rows(): number of rows must be the same" );
    }
  
  
  if( (&out != &A) && (&out != &B) )
    {
    if(join_type == 0)   // join columns (i.e. result matrix has more rows)
      {
      out.set_size(A.n_rows + B.n_rows, A.n_cols);
      
      if( out.n_elem > 0 )
        {
        out.submat(0,        0,   A.n_rows-1, out.n_cols-1) = A;
        out.submat(A.n_rows, 0, out.n_rows-1, out.n_cols-1) = B;
        }
      }
    else   // join rows  (i.e. result matrix has more columns)
      {
      out.set_size(A.n_rows, A.n_cols + B.n_cols);
      
      if( out.n_elem > 0 )
        {
        out.submat(0, 0,        out.n_rows-1,   A.n_cols-1) = A;
        out.submat(0, A.n_cols, out.n_rows-1, out.n_cols-1) = B;
        }
      }
    }
  else  // we have aliasing
    {
    Mat<eT> C;
    
    if(join_type == 0)
      {
      C.set_size(A.n_rows + B.n_rows, A.n_cols);
      
      if( C.n_elem > 0 )
        {
        C.submat(0,        0, A.n_rows-1, C.n_cols-1) = A;
        C.submat(A.n_rows, 0, C.n_rows-1, C.n_cols-1) = B;
        }
      }
    else
      {
      C.set_size(A.n_rows, A.n_cols + B.n_cols);
      
      if( C.n_elem > 0 )
        {
        C.submat(0, 0,        C.n_rows-1, A.n_cols-1) = A;
        C.submat(0, A.n_cols, C.n_rows-1, C.n_cols-1) = B;
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
