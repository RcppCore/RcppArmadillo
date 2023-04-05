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
inline
mat_injector<T1>::mat_injector(T1& in_parent, const typename mat_injector<T1>::elem_type val)
  : parent(in_parent)
  {
  arma_extra_debug_sigprint();
  
  values.reserve(16);
  rowend.reserve(16);
  
  insert(val);
  }



template<typename T1>
inline
mat_injector<T1>::mat_injector(T1& in_parent, const injector_end_of_row<>&)
  : parent(in_parent)
  {
  arma_extra_debug_sigprint();
  
  values.reserve(16);
  rowend.reserve(16);
  
  end_of_row();
  }



template<typename T1>
inline
mat_injector<T1>::~mat_injector()
  {
  arma_extra_debug_sigprint();
  
  const uword N = values.size();
  
  if(N == 0)  { return; }
  
  uword n_rows = 1;
  uword n_cols = 0;
  
  for(uword i=0; i<N; ++i)  { n_rows += (rowend[i]) ? uword(1) : uword(0); }
  
  uword n_cols_in_row = 0;
  
  for(uword i=0; i<N; ++i)
    {
    if(rowend[i])
      {
      n_cols        = (std::max)(n_cols, n_cols_in_row);
      n_cols_in_row = 0;
      }
    else
      {
      ++n_cols_in_row;
      }
    }
  
  n_rows = (rowend[N-1]) ? (n_rows-1) : n_rows;
  n_cols = (std::max)(n_cols, n_cols_in_row);
  
  if(is_Row<T1>::value)
    {
    arma_debug_check( (n_rows > 1), "matrix initialisation: incompatible dimensions" );
    
    parent.zeros(1,n_cols);
    
    uword col = 0;
    
    for(uword i=0; i<N; ++i)
      {
      if(rowend[i])
        {
        break;
        }
      else
        {
        parent.at(col) = values[i];
        ++col;
        }
      }
    }
  else
  if(is_Col<T1>::value)
    {
    const bool is_vec = ((n_cols == 1) || (n_rows == 1));
    
    arma_debug_check( (is_vec == false), "matrix initialisation: incompatible dimensions" );
    
    if(n_cols == 1)
      {
      parent.zeros(n_rows,1);
      
      uword row = 0;
      
      for(uword i=0; i<N; ++i)
        {
        if(rowend[i])
          {
          if((i > 0) && rowend[i-1])  { ++row; }
          }
        else
          {
          parent.at(row) = values[i];
          ++row;
          }
        }
      }
    else
    if(n_rows == 1)
      {
      parent.zeros(n_cols,1);
      
      uword row = 0;
      
      for(uword i=0; i<N; ++i)
        {
        if(rowend[i])
          {
          break;
          }
        else
          {
          parent.at(row) = values[i];
          ++row;
          }
        }
      }
    }
  else
    {
    parent.zeros(n_rows,n_cols);
    
    uword row = 0;
    uword col = 0;
    
    for(uword i=0; i<N; ++i)
      {
      if(rowend[i])
        {
        ++row;
        col = 0;
        }
      else
        {
        parent.at(row,col) = values[i];
        ++col;
        }
      }
    }
  }



template<typename T1>
inline
void
mat_injector<T1>::insert(const typename mat_injector<T1>::elem_type val) const
  {
  arma_extra_debug_sigprint();
  
  values.push_back(val    );
  rowend.push_back(char(0));
  }




template<typename T1>
inline
void
mat_injector<T1>::end_of_row() const
  {
  arma_extra_debug_sigprint();
  
  typedef typename mat_injector<T1>::elem_type eT;
  
  values.push_back(  eT(0));
  rowend.push_back(char(1));
  }



template<typename T1>
inline
const mat_injector<T1>&
operator<<(const mat_injector<T1>& ref, const typename mat_injector<T1>::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  ref.insert(val);
  
  return ref;
  }



template<typename T1>
inline
const mat_injector<T1>&
operator<<(const mat_injector<T1>& ref, const injector_end_of_row<>&)
  {
  arma_extra_debug_sigprint();
  
  ref.end_of_row();
  
  return ref;
  }



//
//
//



template<typename T1>
inline
field_injector<T1>::field_injector(T1& in_parent, const typename field_injector<T1>::object_type& val)
  : parent(in_parent)
  {
  arma_extra_debug_sigprint();
  
  insert(val);
  }



template<typename T1>
inline
field_injector<T1>::field_injector(T1& in_parent, const injector_end_of_row<>&)
  : parent(in_parent)
  {
  arma_extra_debug_sigprint();
  
  end_of_row();
  }



template<typename T1>
inline
field_injector<T1>::~field_injector()
  {
  arma_extra_debug_sigprint();
  
  const uword N = values.size();
  
  if(N == 0)  { return; }
  
  uword n_rows = 1;
  uword n_cols = 0;
  
  for(uword i=0; i<N; ++i)  { n_rows += (rowend[i]) ? uword(1) : uword(0); }
  
  uword n_cols_in_row = 0;
  
  for(uword i=0; i<N; ++i)
    {
    if(rowend[i])
      {
      n_cols        = (std::max)(n_cols, n_cols_in_row);
      n_cols_in_row = 0;
      }
    else
      {
      ++n_cols_in_row;
      }
    }
  
  n_rows = (rowend[N-1]) ? (n_rows-1) : n_rows;
  n_cols = (std::max)(n_cols, n_cols_in_row);
  
  parent.set_size(n_rows,n_cols);
  
  uword row = 0;
  uword col = 0;
  
  for(uword i=0; i<N; ++i)
    {
    if(rowend[i])
      {
      ++row;
      col = 0;
      }
    else
      {
      parent.at(row,col) = std::move(values[i]);
      ++col;
      }
    }
  }



template<typename T1>
inline
void
field_injector<T1>::insert(const typename field_injector<T1>::object_type& val) const
  {
  arma_extra_debug_sigprint();
  
  values.push_back(val    );
  rowend.push_back(char(0));
  }




template<typename T1>
inline
void
field_injector<T1>::end_of_row() const
  {
  arma_extra_debug_sigprint();
  
  typedef typename field_injector<T1>::object_type oT;
  
  values.push_back(oT()   );
  rowend.push_back(char(1));
  }



template<typename T1>
inline
const field_injector<T1>&
operator<<(const field_injector<T1>& ref, const typename field_injector<T1>::object_type& val)
  {
  arma_extra_debug_sigprint();
  
  ref.insert(val);
  
  return ref;
  }



template<typename T1>
inline
const field_injector<T1>&
operator<<(const field_injector<T1>& ref, const injector_end_of_row<>&)
  {
  arma_extra_debug_sigprint();
  
  ref.end_of_row();
  
  return ref;
  }



//! @}
