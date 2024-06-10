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


//! \addtogroup mtSpToDOp
//! @{



template<typename out_eT, typename T1, typename op_type>
class mtSpToDOp : public Base< out_eT, mtSpToDOp<out_eT, T1, op_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  
  typedef typename T1::elem_type                in_eT;
  
  static constexpr bool is_row  = op_type::template traits<T1>::is_row;
  static constexpr bool is_col  = op_type::template traits<T1>::is_col;
  static constexpr bool is_xvec = op_type::template traits<T1>::is_xvec;
  
  inline explicit mtSpToDOp(const T1& in_m);
  inline          mtSpToDOp(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b);
  inline         ~mtSpToDOp();
  
  arma_aligned const T1&   m;            //!< the operand; must be derived from SpBase
  arma_aligned       uword aux_uword_a;  //!< auxiliary data, uword format
  arma_aligned       uword aux_uword_b;  //!< auxiliary data, uword format
  };



//! @}
