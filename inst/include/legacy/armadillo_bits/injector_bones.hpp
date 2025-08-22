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


//! \addtogroup injector
//! @{



template<typename T1>
class mat_injector
  {
  public:
  
  typedef typename T1::elem_type elem_type;
  
  arma_cold inline void  insert(const elem_type val) const;
  arma_cold inline void  end_of_row()                const;
  arma_cold inline      ~mat_injector();
  
  
  private:
  
  inline mat_injector(T1& in_X, const elem_type val);
  inline mat_injector(T1& in_X, const injector_end_of_row<>&);
  
  T1& parent;
  
  mutable std::vector<elem_type> values;
  mutable std::vector<char>      rowend;
  
  friend class Mat<elem_type>;
  friend class Row<elem_type>;
  friend class Col<elem_type>;
  };



//



template<typename T1>
class field_injector
  {
  public:
  
  typedef typename T1::object_type object_type;
  
  arma_cold inline void  insert(const object_type& val) const;
  arma_cold inline void  end_of_row()                   const;
  arma_cold inline      ~field_injector();
  
  
  private:
  
  inline field_injector(T1& in_X, const object_type& val);
  inline field_injector(T1& in_X, const injector_end_of_row<>&);
  
  T1& parent;
  
  mutable std::vector<object_type> values;
  mutable std::vector<char>        rowend;
  
  friend class field<object_type>;
  };



//! @}
