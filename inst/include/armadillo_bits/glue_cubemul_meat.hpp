// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup glue_cubemul
//! @{



template<typename T1, typename T2>
inline
void
glue_cubemul::apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_cubemul>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const strip_trans<T1> SA(X.A);
  const strip_trans<T2> SB(X.B);
  
  constexpr bool do_strans_A = strip_trans<T1>::do_strans;
  constexpr bool do_htrans_A = strip_trans<T1>::do_htrans;
  
  constexpr bool do_strans_B = strip_trans<T2>::do_strans;
  constexpr bool do_htrans_B = strip_trans<T2>::do_htrans;
  
  const unwrap_cube<typename strip_trans<T1>::stored_type> UA(SA.M);
  const unwrap_cube<typename strip_trans<T2>::stored_type> UB(SB.M);
  
  const Cube<eT>& A = UA.M;
  const Cube<eT>& B = UB.M;
  
  if(UA.is_alias(out) || UB.is_alias(out))
    {
    Cube<eT> tmp;
    
    glue_cubemul::apply_noalias<eT, do_strans_A, do_htrans_A, do_strans_B, do_htrans_B>(tmp, A, B);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_cubemul::apply_noalias<eT, do_strans_A, do_htrans_A, do_strans_B, do_htrans_B>(out, A, B);
    }
  }



template<typename eT, bool do_strans_A, bool do_htrans_A, bool do_strans_B, bool do_htrans_B>
inline
void
glue_cubemul::apply_noalias(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B)
  {
  arma_debug_sigprint();
  
  const uword A_n_slices = A.n_slices;
  
  arma_conform_check( (A_n_slices != B.n_slices), "cubemul(): given cubes must have the same number of slices" );
  
  constexpr bool do_xtrans_A = (do_strans_A || do_htrans_A);
  constexpr bool do_xtrans_B = (do_strans_B || do_htrans_B);
  
  const uword final_A_n_rows = (do_xtrans_A == false) ? A.n_rows : A.n_cols;
  const uword final_A_n_cols = (do_xtrans_A == false) ? A.n_cols : A.n_rows;
  
  const uword final_B_n_rows = (do_xtrans_B == false) ? B.n_rows : B.n_cols;
  const uword final_B_n_cols = (do_xtrans_B == false) ? B.n_cols : B.n_rows;
  
  if( (arma_config::check_conform) && (final_A_n_cols != final_B_n_rows) )
    {
    std::ostringstream tmp;
    
    tmp << "cubemul(): incompatible dimensions: " << final_A_n_rows << 'x' << final_A_n_cols << 'x' << A.n_slices << " and " << final_B_n_rows << 'x' << final_B_n_cols << 'x' << B.n_slices;
    
    arma_stop_logic_error(tmp.str());
    }
  
  out.set_size(final_A_n_rows, final_B_n_cols, A_n_slices);
  
  if(out.is_empty())  { return; }
  
  for(uword s=0; s < A_n_slices; ++s)
    {
    const Mat<eT> A_slice_s(const_cast<eT*>(A.slice_memptr(s)), A.n_rows, A.n_cols, false, true);
    const Mat<eT> B_slice_s(const_cast<eT*>(B.slice_memptr(s)), B.n_rows, B.n_cols, false, true);
    
    Mat<eT> out_slice_s(out.slice_memptr(s), final_A_n_rows, final_B_n_cols, false, true);
    
    if( (do_xtrans_A == false) && (do_xtrans_B == false) )
      {
      out_slice_s = A_slice_s * B_slice_s;
      }
    else
    if( (do_xtrans_A == false) && (do_xtrans_B == true ) )
      {
           if(do_strans_B)  { out_slice_s = A_slice_s * strans(B_slice_s); }
      else if(do_htrans_B)  { out_slice_s = A_slice_s * htrans(B_slice_s); }
      }
    else
    if( (do_xtrans_A == true ) && (do_xtrans_B == false) )
      {
           if(do_strans_A)  { out_slice_s = strans(A_slice_s) * B_slice_s; }
      else if(do_htrans_A)  { out_slice_s = htrans(A_slice_s) * B_slice_s; }
      }
    else
    if( (do_xtrans_A == true ) && (do_xtrans_B == true ) )
      {
           if( (do_strans_A) && (do_strans_B) )  { out_slice_s = strans(A_slice_s) * strans(B_slice_s); }
      else if( (do_strans_A) && (do_htrans_B) )  { out_slice_s = strans(A_slice_s) * htrans(B_slice_s); }
      else if( (do_htrans_A) && (do_strans_B) )  { out_slice_s = htrans(A_slice_s) * strans(B_slice_s); }
      else if( (do_htrans_A) && (do_htrans_B) )  { out_slice_s = htrans(A_slice_s) * htrans(B_slice_s); }
      }
    }
  }



//



template<typename T1, typename T2>
inline
Cube<typename T1::elem_type>
glue_cubemul::apply(const BaseCube<typename T1::elem_type,T1>& expr_A, const Base<typename T1::elem_type,T2>& expr_B)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const strip_trans<T1> SA(expr_A.get_ref());
  const strip_trans<T2> SB(expr_B.get_ref());
  
  constexpr bool do_strans_A = strip_trans<T1>::do_strans;
  constexpr bool do_htrans_A = strip_trans<T1>::do_htrans;
  
  constexpr bool do_strans_B = strip_trans<T2>::do_strans;
  constexpr bool do_htrans_B = strip_trans<T2>::do_htrans;
  
  const  unwrap_cube<typename strip_trans<T1>::stored_type> UA(SA.M);
  const quasi_unwrap<typename strip_trans<T2>::stored_type> UB(SB.M);
  
  const Cube<eT>& A = UA.M;
  const  Mat<eT>& B = UB.M;
  
  constexpr bool do_xtrans_A = (do_strans_A || do_htrans_A);
  constexpr bool do_xtrans_B = (do_strans_B || do_htrans_B);
  
  const uword final_A_n_rows = (do_xtrans_A == false) ? A.n_rows : A.n_cols;
  const uword final_A_n_cols = (do_xtrans_A == false) ? A.n_cols : A.n_rows;
  
  const uword final_B_n_rows = (do_xtrans_B == false) ? B.n_rows : B.n_cols;
  const uword final_B_n_cols = (do_xtrans_B == false) ? B.n_cols : B.n_rows;
  
  if( (arma_config::check_conform) && (final_A_n_cols != final_B_n_rows) )
    {
    std::ostringstream tmp;
    
    tmp << "cubemul(): incompatible dimensions: " << final_A_n_rows << 'x' << final_A_n_cols << 'x' << A.n_slices << " and " << final_B_n_rows << 'x' << final_B_n_cols;
    
    arma_stop_logic_error(tmp.str());
    }
  
  Cube<eT> out(final_A_n_rows, final_B_n_cols, A.n_slices, arma_nozeros_indicator());
  
  if(out.is_empty() == false)
    {
    for(uword s=0; s < A.n_slices; ++s)
      {
      const Mat<eT> A_slice_s(const_cast<eT*>(A.slice_memptr(s)), A.n_rows, A.n_cols, false, true);
      
      Mat<eT> out_slice_s(out.slice_memptr(s), final_A_n_rows, final_B_n_cols, false, true);
      
      if( (do_xtrans_A == false) && (do_xtrans_B == false) )
        {
        out_slice_s = A_slice_s * B;
        }
      else
      if( (do_xtrans_A == false) && (do_xtrans_B == true ) )
        {
             if(do_strans_B)  { out_slice_s = A_slice_s * strans(B); }
        else if(do_htrans_B)  { out_slice_s = A_slice_s * htrans(B); }
        }
      else
      if( (do_xtrans_A == true ) && (do_xtrans_B == false) )
        {
             if(do_strans_A)  { out_slice_s = strans(A_slice_s) * B; }
        else if(do_htrans_A)  { out_slice_s = htrans(A_slice_s) * B; }
        }
      else
      if( (do_xtrans_A == true ) && (do_xtrans_B == true ) )
        {
             if( (do_strans_A) && (do_strans_B) )  { out_slice_s = strans(A_slice_s) * strans(B); }
        else if( (do_strans_A) && (do_htrans_B) )  { out_slice_s = strans(A_slice_s) * htrans(B); }
        else if( (do_htrans_A) && (do_strans_B) )  { out_slice_s = htrans(A_slice_s) * strans(B); }
        else if( (do_htrans_A) && (do_htrans_B) )  { out_slice_s = htrans(A_slice_s) * htrans(B); }
        }
      }
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
Cube<typename T1::elem_type>
glue_cubemul::apply(const Base<typename T1::elem_type,T1>& expr_A, const BaseCube<typename T1::elem_type,T2>& expr_B)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const strip_trans<T1> SA(expr_A.get_ref());
  const strip_trans<T2> SB(expr_B.get_ref());
  
  constexpr bool do_strans_A = strip_trans<T1>::do_strans;
  constexpr bool do_htrans_A = strip_trans<T1>::do_htrans;
  
  constexpr bool do_strans_B = strip_trans<T2>::do_strans;
  constexpr bool do_htrans_B = strip_trans<T2>::do_htrans;
  
  const quasi_unwrap<typename strip_trans<T1>::stored_type> UA(SA.M);
  const  unwrap_cube<typename strip_trans<T2>::stored_type> UB(SB.M);
  
  const  Mat<eT>& A = UA.M;
  const Cube<eT>& B = UB.M;
  
  constexpr bool do_xtrans_A = (do_strans_A || do_htrans_A);
  constexpr bool do_xtrans_B = (do_strans_B || do_htrans_B);
  
  const uword final_A_n_rows = (do_xtrans_A == false) ? A.n_rows : A.n_cols;
  const uword final_A_n_cols = (do_xtrans_A == false) ? A.n_cols : A.n_rows;
  
  const uword final_B_n_rows = (do_xtrans_B == false) ? B.n_rows : B.n_cols;
  const uword final_B_n_cols = (do_xtrans_B == false) ? B.n_cols : B.n_rows;
  
  if( (arma_config::check_conform) && (final_A_n_cols != final_B_n_rows) )
    {
    std::ostringstream tmp;
    
    tmp << "cubemul(): incompatible dimensions: " << final_A_n_rows << 'x' << final_A_n_cols << " and " << final_B_n_rows << 'x' << final_B_n_cols << 'x' << B.n_slices;
    
    arma_stop_logic_error(tmp.str());
    }
  
  Cube<eT> out(final_A_n_rows, final_B_n_cols, B.n_slices, arma_nozeros_indicator());
  
  if(out.is_empty() == false)
    {
    for(uword s=0; s < B.n_slices; ++s)
      {
      const Mat<eT> B_slice_s(const_cast<eT*>(B.slice_memptr(s)), B.n_rows, B.n_cols, false, true);
      
      Mat<eT> out_slice_s(out.slice_memptr(s), final_A_n_rows, final_B_n_cols, false, true);
      
      if( (do_xtrans_A == false) && (do_xtrans_B == false) )
        {
        out_slice_s = A * B_slice_s;
        }
      else
      if( (do_xtrans_A == false) && (do_xtrans_B == true ) )
        {
             if(do_strans_B)  { out_slice_s = A * strans(B_slice_s); }
        else if(do_htrans_B)  { out_slice_s = A * htrans(B_slice_s); }
        }
      else
      if( (do_xtrans_A == true ) && (do_xtrans_B == false) )
        {
             if(do_strans_A)  { out_slice_s = strans(A) * B_slice_s; }
        else if(do_htrans_A)  { out_slice_s = htrans(A) * B_slice_s; }
        }
      else
      if( (do_xtrans_A == true ) && (do_xtrans_B == true ) )
        {
             if( (do_strans_A) && (do_strans_B) )  { out_slice_s = strans(A) * strans(B_slice_s); }
        else if( (do_strans_A) && (do_htrans_B) )  { out_slice_s = strans(A) * htrans(B_slice_s); }
        else if( (do_htrans_A) && (do_strans_B) )  { out_slice_s = htrans(A) * strans(B_slice_s); }
        else if( (do_htrans_A) && (do_htrans_B) )  { out_slice_s = htrans(A) * htrans(B_slice_s); }
        }
      }
    }
  
  return out;
  }



//! @}
