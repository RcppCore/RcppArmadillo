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


//! \addtogroup mtSpReduceOp
//! @{


// NOTE: mtSpReduceOp is dedicated for reduction operations on sparse matrices,
// NOTE: including sum(), min(), max(), mean(), var(), stddev(),
// NOTE  where the entire sparse matrix is reduced to a vector.
// NOTE: 
// NOTE: Even though it would make more sense for mtSpReduceOp to be derived from Base
// NOTE: (as it's more efficient to store the resulting vectors in dense format),
// NOTE: mtSpReduceOp is derived from SpBase so the default user-accessible output is in sparse storage format.
// NOTE. This is to mimic the behaviour of Octave and to keep compatibility with existing user code.
// NOTE: 
// NOTE: However, for simplicity and efficiency, each mtSpReduceOp op_type::apply() function
// NOTE: must output to a dense matrix, ie. apply(Mat, ...).
// NOTE: The SpMat class handles mtSpReduceOp by converting the dense output to a sparse representation.
// NOTE: The Mat class has an explicit constructor to efficiently handle mtSpReduceOp.

template<typename out_eT, typename T1, typename op_type>
class mtSpReduceOp : public SpBase< out_eT, mtSpReduceOp<out_eT, T1, op_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  
  typedef typename T1::elem_type                in_eT;
  
  static constexpr bool is_row  = op_type::template traits<T1>::is_row;
  static constexpr bool is_col  = op_type::template traits<T1>::is_col;
  static constexpr bool is_xvec = op_type::template traits<T1>::is_xvec;
  
  inline explicit mtSpReduceOp(const T1& in_m);
  inline          mtSpReduceOp(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b);
  inline         ~mtSpReduceOp();
  
  arma_aligned const T1&   m;            //!< the operand; must be derived from SpBase
  arma_aligned       uword aux_uword_a;  //!< auxiliary data, uword format
  arma_aligned       uword aux_uword_b;  //!< auxiliary data, uword format
  };



//! @}
