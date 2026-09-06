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



//! \addtogroup op_permute
//! @{



template<typename T1>
inline
void
op_permute::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_permute>& in)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> U(in.m);
  const Cube<eT>& X   = U.M;
  
  const uword d0 = in.aux_uword_a;
  const uword d1 = in.aux_uword_b;
  const uword d2 = in.aux_uword_c;
  
  const bool mapping_ok = (d0 <= uword(2)) && (d1 <= uword(2)) && (d2 <= uword(2)) && (d0 != d1) && (d0 != d2) && (d1 != d2);
  
  arma_conform_check( (mapping_ok == false), "permute(): each mapping dimension must be unique and one of: 0, 1, 2" );
  
  if(U.is_alias(out))
    {
    Cube<eT> tmp;
    
    op_permute::apply_noalias(tmp, X, d0, d1, d2);
    
    out.steal_mem(tmp);
    }
  else
    {
    op_permute::apply_noalias(out, X, d0, d1, d2);
    }
  }



template<typename eT>
inline
void
op_permute::apply_noalias(Cube<eT>& Y, const Cube<eT>& X, const uword d0, const uword d1, const uword d2)
  {
  arma_debug_sigprint();
  
  const uword X_nr = X.n_rows;
  const uword X_nc = X.n_cols;
  const uword X_ns = X.n_slices;

  const uword Y_nr = (d0 == uword(0)) ? X_nr : ((d0 == uword(1)) ? X_nc : X_ns);
  const uword Y_nc = (d1 == uword(0)) ? X_nr : ((d1 == uword(1)) ? X_nc : X_ns);
  const uword Y_ns = (d2 == uword(0)) ? X_nr : ((d2 == uword(1)) ? X_nc : X_ns);
  
  Y.set_size(Y_nr, Y_nc, Y_ns);
  
  if(Y.is_empty())  { return; }
  
  // {0,1,2} -> identity
  // {0,2,1} -> swap columns/slices
  // {1,0,2} -> transpose each slice
  // {1,2,0} -> cyclic permutation
  // {2,0,1} -> cyclic permutation
  // {2,1,0} -> reverse/transpose-style permutation
  
  if( (d0 == uword(0)) && (d1 == uword(1)) && (d2 == uword(2)) )
    {
    arma_debug_print("mapping {0,1,2}");
    
    Y = X;
    }
  else
  if( (d0 == uword(0)) && (d1 == uword(2)) && (d2 == uword(1)) )
    {
    arma_debug_print("mapping {0,2,1}");
    
    // for(uword s=0; s < X_ns; ++s)
    // for(uword c=0; c < X_nc; ++c)
    //   {
    //   Y.slice(c).col(s) = X.slice(s).col(c);
    //   }
    
    for(uword s=0; s < X_ns; ++s)
      {
      const Mat<eT> X_slice_s(const_cast<eT*>(X.slice_memptr(s)), X_nr, X_nc, false, true);
      
      for(uword c=0; c < X_nc; ++c)
        {
        Mat<eT> Y_slice_c(Y.slice_memptr(c), Y_nr, Y_nc, false, true);
        
        arrayops::copy(Y_slice_c.colptr(s), X_slice_s.colptr(c), X_nr);
        }
      }
    }
  else
  if( (d0 == uword(1)) && (d1 == uword(0)) && (d2 == uword(2)) )
    {
    arma_debug_print("mapping {1,0,2}");
    
    // for(uword s=0; s < X_ns; ++s)
    //   {
    //   Y.slice(s) = X.slice(s).st();
    //   }
    
    for(uword s=0; s < X_ns; ++s)
      {
      const Mat<eT> X_slice_s(const_cast<eT*>(X.slice_memptr(s)), X_nr, X_nc, false, true);
            Mat<eT> Y_slice_s(                Y.slice_memptr(s) , Y_nr, Y_nc, false, true);
      
      op_strans::apply_mat_noalias(Y_slice_s, X_slice_s);
      }
    }
  else
  if( (d0 == uword(1)) && (d1 == uword(2)) && (d2 == uword(0)) )
    {
    arma_debug_print("mapping {1,2,0}");
    
    // better-than-nothing implementation
    
    for(uword s=0; s < X_ns; ++s)
    for(uword c=0; c < X_nc; ++c)
    for(uword r=0; r < X_nr; ++r)
      {
      Y.at(c, s, r) = X.at(r, c, s);
      }
    }
  else
  if( (d0 == uword(2)) && (d1 == uword(0)) && (d2 == uword(1)) )
    {
    arma_debug_print("mapping {2,0,1}");
    
    // better-than-nothing implementation
    
    for(uword s=0; s < X_ns; ++s)
    for(uword c=0; c < X_nc; ++c)
    for(uword r=0; r < X_nr; ++r)
      {
      Y.at(s, r, c) = X.at(r, c, s);
      }
    }
  else
  if( (d0 == uword(2)) && (d1 == uword(1)) && (d2 == uword(0)) )
    {
    arma_debug_print("mapping {2,1,0}");
    
    // better-than-nothing implementation
    
    for(uword s=0; s < X_ns; ++s)
    for(uword c=0; c < X_nc; ++c)
    for(uword r=0; r < X_nr; ++r)
      {
      Y.at(s, c, r) = X.at(r, c, s);
      }
    }
  else
    {
    arma_stop_logic_error("permute(): unknown permutation");
    }
  }



//! @}
