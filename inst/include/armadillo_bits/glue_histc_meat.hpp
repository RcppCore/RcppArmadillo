// Copyright (C) 2015 Conrad Sanderson
// Copyright (C) 2015 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



template<typename T1, typename T2>
inline
void
glue_histc::apply(Mat<uword>& C, const mtGlue<uword,T1,T2,glue_histc>& expr)
  {
  arma_extra_debug_sigprint();
  
  const uword dim = expr.aux_uword;
  
  arma_debug_check( (dim > 1), "histc(): parameter 'dim' must be 0 or 1" );
  
  const unwrap_check_mixed<T1> tmpA(expr.A, C);
  const unwrap_check_mixed<T2> tmpB(expr.B, C);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmpA.M;
  const Mat<eT>& B = tmpB.M;
  
  arma_debug_check( ((B.is_vec() == false) && (B.is_empty() == false)), "histc(): parameter 'edges' is not a vector" );
  
  const uword A_n_rows = A.n_rows;
  const uword A_n_cols = A.n_cols;
  
  const uword A_n_elem = A.n_elem;
  const uword B_n_elem = B.n_elem;
  
  if( B_n_elem == uword(0) )  { C.reset(); return; }
  
  // NOTE: the "dim" parameter is currently ignored for vectors
  
  uword C_n_rows = uword(0);
  uword C_n_cols = uword(0);
  
  if( (A.n_elem == uword(1)) && (A.vec_state == uword(0)) )
    {
    if(C.vec_state == uword(1))
      {
      C_n_rows = B_n_elem;
      C_n_cols = uword(1);
      }
    else
      {
      C_n_rows = uword(1);
      C_n_cols = B_n_elem;
      }
    }
  else
  if( A.is_vec() || (A.vec_state > 0) )
    {
    if(A.vec_state == uword(2))
      {
      C_n_rows = uword(1);
      C_n_cols = B_n_elem;
      }
    else
    if(A.vec_state == uword(1))
      {
      C_n_rows = B_n_elem;
      C_n_cols = uword(1);
      }
    else
    if(A.is_rowvec())
      {
      C_n_rows = uword(1);
      C_n_cols = B_n_elem;
      }
    else
    if(A.is_colvec())
      {
      C_n_rows = B_n_elem;
      C_n_cols = uword(1);
      }
    }
  else
    {
    if(dim == uword(0))
      {
      C_n_rows = B_n_elem;
      C_n_cols = A_n_cols;
      }
    else
    if(dim == uword(1))
      {
      C_n_rows = A_n_rows;
      C_n_cols = B_n_elem;
      }
    }
  
  
  C.zeros(C_n_rows, C_n_cols);
  
  
  const eT*   B_mem       = B.memptr();
  const uword B_n_elem_m1 = B_n_elem - 1;
  
  if( A.is_vec() || (A.vec_state > 0) )
    {
    const eT*    A_mem = A.memptr();
          uword* C_mem = C.memptr();
    
    for(uword j=0; j < A_n_elem; ++j)
      {
      const eT x = A_mem[j];
      
      for(uword i=0; i < B_n_elem_m1; ++i)
        {
             if( (B_mem[i]           <= x) && (x < B_mem[i+1]) )  { C_mem[i]++;           break; }
        else if(  B_mem[B_n_elem_m1] == x                      )  { C_mem[B_n_elem_m1]++; break; }    // for compatibility with Matlab
        }
      }
    }
  else
  if(dim == uword(1))
    {
    for(uword row=0; row < A_n_rows; ++row)
      {
      for(uword col=0; col < A_n_cols; ++col)
        {
        const eT x = A.at(row,col);
        
        for(uword i=0; i < B_n_elem_m1; ++i)
          {
               if( (B_mem[i]            <= x) && (x < B_mem[i+1]) )  { C.at(row,i)++;           break; }
          else if(  B_mem[B_n_elem_m1]  == x                      )  { C.at(row,B_n_elem_m1)++; break; }   // for compatibility with Matlab
          }
        }
      }
    }
  else
  if(dim == uword(0))
    {
    for(uword col=0; col < A_n_cols; ++col)
      {
      const eT*    A_coldata = A.colptr(col);
            uword* C_coldata = C.colptr(col);
      
      for(uword row=0; row < A_n_rows; ++row)
        {
        const eT x = A_coldata[row];
        
        for(uword i=0; i < B_n_elem_m1; ++i)
          {
               if( (B_mem[i]           <= x) && (x < B_mem[i+1]) )  { C_coldata[i]++;           break; }
          else if(  B_mem[B_n_elem_m1] == x                      )  { C_coldata[B_n_elem_m1]++; break; }    // for compatibility with Matlab
          }
        }
      }
    }
  }
