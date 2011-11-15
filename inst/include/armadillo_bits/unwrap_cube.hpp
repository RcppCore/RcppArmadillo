// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup unwrap_cube
//! @{



template<typename T1>
class unwrap_cube
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_cube(const T1& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Cube<eT> M;
  };



template<typename eT>
class unwrap_cube< Cube<eT> >
  {
  public:

  inline
  unwrap_cube(const Cube<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Cube<eT>& M;
  };



//
//
//



template<typename T1>
class unwrap_cube_check
  {
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_cube_check(const T1& A, const Cube<eT>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    
    arma_type_check(( is_arma_cube_type<T1>::value == false ));
    }
  
  const Cube<eT> M;
  };



template<typename eT>
class unwrap_cube_check< Cube<eT> >
  {
  public:

  inline
  unwrap_cube_check(const Cube<eT>& A, const Cube<eT>& B)
    : M_local( (&A == &B) ? new Cube<eT>(A) : 0 )
    , M      ( (&A == &B) ? (*M_local)      : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_cube_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  // the order below is important
  const Cube<eT>* M_local;
  const Cube<eT>& M;
  
  };



//! @}
