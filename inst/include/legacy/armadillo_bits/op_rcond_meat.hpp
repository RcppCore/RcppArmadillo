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


//! \addtogroup op_rcond
//! @{



template<typename T1>
inline
typename T1::pod_type
op_rcond::apply(const Base<typename T1::elem_type, T1>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  typedef typename T1::pod_type   T;
  
  if(strip_trimat<T1>::do_trimat)
    {
    arma_debug_print("op_rcond::apply(): tri optimisation");
    
    const strip_trimat<T1> S(X.get_ref());
    
    const quasi_unwrap<typename strip_trimat<T1>::stored_type> U(S.M);
    
    arma_conform_check( (U.M.is_square() == false), "rcond(): matrix must be square sized" );
    
    const uword layout = (S.do_triu) ? uword(0) : uword(1);
    
    return auxlib::rcond_trimat(U.M, layout);
    }
  
  Mat<eT> A = X.get_ref();
  
  arma_conform_check( (A.is_square() == false), "rcond(): matrix must be square sized" );
  
  if(A.is_empty()) { return Datum<T>::inf; }
  
  if(is_op_diagmat<T1>::value || A.is_diagmat())
    {
    arma_debug_print("op_rcond::apply(): diag optimisation");
    
    const eT*   colmem = A.memptr();
    const uword N      = A.n_rows;
    
    T abs_min = Datum<T>::inf;
    T abs_max = T(0);
    
    for(uword i=0; i<N; ++i)
      {
      const T abs_val = std::abs(colmem[i]);
      
      abs_min = (abs_val < abs_min) ? abs_val : abs_min;
      abs_max = (abs_val > abs_max) ? abs_val : abs_max;
      
      colmem += N;
      }
    
    if((abs_min == T(0)) || (abs_max == T(0)))  { return T(0); }
    
    return T(abs_min / abs_max);
    }
  
  const bool is_triu =                     trimat_helper::is_triu(A);
  const bool is_tril = (is_triu) ? false : trimat_helper::is_tril(A);
  
  if(is_triu || is_tril)
    {
    arma_debug_print("op_rcond::apply(): tri optimisation");
    
    const uword layout = (is_triu) ? uword(0) : uword(1);
    
    return auxlib::rcond_trimat(A, layout);
    }
  
  if( (arma_config::optimise_sym) && (auxlib::crippled_lapack(A) == false) && ( is_sym_expr<T1>::eval(X.get_ref()) || sym_helper::is_approx_sym(A, uword(100)) ) )
    {
    arma_debug_print("op_rcond::apply(): symmetric/hermitian optimisation");
    
    return auxlib::rcond_sym(A);
    }
  
  return auxlib::rcond(A);
  }



//! @}
