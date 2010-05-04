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
      
      out.submat(0,        0, A.n_rows-1,   out.n_cols-1) = A;
      out.submat(A.n_rows, 0, out.n_rows-1, out.n_cols-1) = B;
      }
    else   // join rows  (i.e. result matrix has more columns)
      {
      out.set_size(A.n_rows, A.n_cols + B.n_cols);
      
      out.submat(0, 0,        out.n_rows-1, A.n_cols-1  ) = A;
      out.submat(0, A.n_cols, out.n_rows-1, out.n_cols-1) = B;
      }
    }
  else  // we have aliasing
    {
    Mat<eT> C;
    
    if(join_type == 0)
      {
      C.set_size(A.n_rows + B.n_rows, A.n_cols);
      
      C.submat(0,        0, A.n_rows-1, C.n_cols-1) = A;
      C.submat(A.n_rows, 0, C.n_rows-1, C.n_cols-1) = B;
      }
    else
      {
      C.set_size(A.n_rows, A.n_cols + B.n_cols);
      
      C.submat(0, 0,        C.n_rows-1, A.n_cols-1) = A;
      C.submat(0, A.n_cols, C.n_rows-1, C.n_cols-1) = B;
      }
    
    if(C.n_elem > sizeof(C.mem_local)/sizeof(eT))
      {
      out.reset();
      
      access::rw(out.n_elem) = C.n_elem;
      access::rw(out.n_rows) = C.n_rows;
      access::rw(out.n_cols) = C.n_cols;
      access::rw(out.mem   ) = C.mem;
      
      access::rw(C.n_elem) = 0;
      access::rw(C.n_rows) = 0;
      access::rw(C.n_cols) = 0;
      }
    else
      {
      out = C;
      }
    }
  
  }



//! @}
