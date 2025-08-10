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


//! \addtogroup fn_balance
//! @{



template<typename T1>
inline
typename enable_if2< is_blas_type<typename T1::elem_type>::value, bool >::result
balance(Col<typename T1::pod_type>& S, Col<uword>& P, Mat<typename T1::elem_type>& B, const Base<typename T1::elem_type,T1>& A, const char* method = "both")
  {
  arma_debug_sigprint();
  
  arma_conform_check( (void_ptr(&S)  == void_ptr(&B)), "eig_gen(): parameter 'S' is an alias of parameter 'B'" );
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  if( (sig != 'b') && (sig != 's') && (sig != 'p') )  { arma_stop_logic_error("balance(): unsupported method"); }
  
  const bool do_scale = (sig == 'b') || (sig == 's');
  const bool do_perm  = (sig == 'b') || (sig == 'p');
  
  const bool calc_SP = true;
  
  B = A.get_ref();
  
  if(B.is_square() == false)
    {
    B.soft_reset();
    
    arma_stop_logic_error("balance(): given matrix must be square sized");
    
    return false;
    }
  
  const bool status = auxlib::balance(S, P, B, calc_SP, do_scale, do_perm);
  
  if(status == false)
    {
    S.soft_reset();
    P.soft_reset();
    B.soft_reset();
    
    arma_warn(3, "balance(): transformation failed");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_blas_type<typename T1::elem_type>::value, bool >::result
balance(Mat<typename T1::elem_type>& B, const Base<typename T1::elem_type,T1>& A, const char* method = "both")
  {
  arma_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  const char sig = (method != nullptr) ? method[0] : char(0);
  
  if( (sig != 'b') && (sig != 's') && (sig != 'p') )  { arma_stop_logic_error("balance(): unsupported method"); }
  
  const bool do_scale = (sig == 'b') || (sig == 's');
  const bool do_perm  = (sig == 'b') || (sig == 'p');
  
  const bool calc_SP = false;
  
  B = A.get_ref();
  
  if(B.is_square() == false)
    {
    B.soft_reset();
    
    arma_stop_logic_error("balance(): given matrix must be square sized");
    
    return false;
    }
  
  Col<T>     S;
  Col<uword> P;
  
  const bool status = auxlib::balance(S, P, B, calc_SP, do_scale, do_perm);
  
  if(status == false)
    {
    B.soft_reset();
    
    arma_warn(3, "balance(): transformation failed");
    }
  
  return status;
  }



template<typename T1>
inline
typename enable_if2< is_blas_type<typename T1::elem_type>::value, Mat<typename T1::elem_type> >::result
balance(const Base<typename T1::elem_type,T1>& A, const char* method = "both")
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> B;
  
  const bool status = balance(B, A, method);
  
  if(status == false)  { arma_stop_runtime_error("balance(): transformation failed"); }
  
  return B;
  }



//! @}
