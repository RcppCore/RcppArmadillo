// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup unwrap_cube
//! @{



template<typename T1>
struct unwrap_cube
  {
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_cube(const T1& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Cube<eT> M;
  
  template<typename eT2>
  constexpr bool is_alias(const Cube<eT2>&) const { return false; }
  };



template<typename eT>
struct unwrap_cube< Cube<eT> >
  {
  inline
  unwrap_cube(const Cube<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Cube<eT>& M;
  
  template<typename eT2>
  arma_inline bool is_alias(const Cube<eT2>& X) const { return (void_ptr(&M) == void_ptr(&X)); }
  };



//
//
//



template<typename T1>
struct unwrap_cube_check
  {
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_cube_check(const T1& A, const Cube<eT>&)
    : M(A)
    {
    arma_extra_debug_sigprint();
    
    arma_type_check(( is_arma_cube_type<T1>::value == false ));
    }
  
  inline
  unwrap_cube_check(const T1& A, const bool)
    : M(A)
    {
    arma_extra_debug_sigprint();
    
    arma_type_check(( is_arma_cube_type<T1>::value == false ));
    }
  
  const Cube<eT> M;
  };



template<typename eT>
struct unwrap_cube_check< Cube<eT> >
  {
  inline
  unwrap_cube_check(const Cube<eT>& A, const Cube<eT>& B)
    : M_local( (&A == &B) ? new Cube<eT>(A) : nullptr )
    , M      ( (&A == &B) ? (*M_local)      : A       )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_cube_check(const Cube<eT>& A, const bool is_alias)
    : M_local( is_alias ? new Cube<eT>(A) : nullptr )
    , M      ( is_alias ? (*M_local)      : A       )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_cube_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)  { delete M_local; }
    }
  
  
  // the order below is important
  const Cube<eT>* M_local;
  const Cube<eT>& M;
  };



//! @}
