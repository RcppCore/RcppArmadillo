// Copyright (C) 2012 Ryan Curtin
// Copyright (C) 2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup spop_max
//! @{



template<typename T1>
inline
void
spop_max::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_max>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  arma_debug_check((dim > 1), "max(): incorrect usage. dim must be 0 or 1");
  
  const SpProxy<T1> p(in.m);
  
  if(p.is_alias(out) == false)
    {
    spop_max::apply_noalias(out, p, dim);
    }
  else
    {
    SpMat<eT> tmp;
    
    spop_max::apply_noalias(tmp, p, dim);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1>
inline
void
spop_max::apply_noalias
  (
        SpMat<typename T1::elem_type>& result,
  const SpProxy<T1>&                   p,
  const uword                          dim,
  const typename arma_not_cx<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;

  if(dim == 0)
    {
    // maximum in each column
    result.set_size(1, p.get_n_cols());
    
    if(p.get_n_nonzero() == 0)
      {
      return;
      }

    typename SpProxy<T1>::const_iterator_type it = p.begin();

    uword cur_col = it.col();
    uword elem_in_col = 1;
    eT cur_max = (*it);
    ++it;

    while(it != p.end())
      {
      if(it.col() != cur_col)
        {
        // was the column full?
        if(elem_in_col == p.get_n_rows())
          {
          result.at(0, cur_col) = cur_max;
          }
        else
          {
          result.at(0, cur_col) = std::max(eT(0), cur_max);
          }

        cur_col = it.col();
        elem_in_col = 0;
        cur_max = (*it);
        }
      else
        {
        cur_max = std::max(cur_max, *it);
        }

      ++elem_in_col;
      ++it;
      }

    if(elem_in_col == p.get_n_rows())
      {
      result.at(0, cur_col) = cur_max;
      }
    else
      {
      result.at(0, cur_col) = std::max(eT(0), cur_max);
      }
    }
  else
    {
    // maximum in each row
    result.set_size(p.get_n_rows(), 1);
  
    if(p.get_n_nonzero() == 0)
      {
      return;
      }

    typename SpProxy<T1>::const_row_iterator_type it = p.begin_row();

    uword cur_row = it.row();
    uword elem_in_row = 1;
    eT cur_max = (*it);
    ++it;

    while(it.pos() < p.get_n_nonzero())
      {
      if(it.row() != cur_row)
        {
        // was the row full?
        if(elem_in_row == p.get_n_cols())
          {
          result.at(cur_row, 0) = cur_max;
          }
        else
          {
          result.at(cur_row, 0) = std::max(eT(0), cur_max);
          }

        cur_row = it.row();
        elem_in_row = 0;
        cur_max = (*it);
        }
      else
        {
        cur_max = std::max(cur_max, *it);
        }

      ++elem_in_row;
      ++it;
      }

    if(elem_in_row == p.get_n_cols())
      {
      result.at(cur_row, 0) = cur_max;
      }
    else
      {
      result.at(cur_row, 0) = std::max(eT(0), cur_max);
      }
    }
  }



template<typename T1>
inline
void
spop_max::apply_noalias
  (
        SpMat<typename T1::elem_type>& result,
  const SpProxy<T1>&                   p,
  const uword                          dim,
  const typename arma_cx_only<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type            eT;
  typedef typename get_pod_type<eT>::result  T;
  
  if(dim == 0)
    {
    // maximum in each column
    result.set_size(1, p.get_n_cols());
  
    if(p.get_n_nonzero() == 0)
      {
      return;
      }
    
    typename SpProxy<T1>::const_iterator_type it = p.begin();
    
    uword cur_col = it.col();
    uword elem_in_col = 1;
    
    eT cur_max_orig = *it;
     T cur_max_abs  = std::abs(cur_max_orig);
    
    ++it;
    
    while(it != p.end())
      {
      if(it.col() != cur_col)
        {
        // was the column full?
        if(elem_in_col == p.get_n_rows())
          {
          result.at(0, cur_col) = cur_max_orig;
          }
        else
          {
          eT val1 = eT(0);
          
          result.at(0, cur_col) = ( std::abs(val1) >= cur_max_abs ) ? val1 : cur_max_orig;
          }
        
        cur_col = it.col();
        elem_in_col = 0;
        
        cur_max_orig = *it;
        cur_max_abs  = std::abs(cur_max_orig);
        }
      else
        {
        eT val1_orig = *it;
         T val1_abs  = std::abs(val1_orig);
        
        if( val1_abs >= cur_max_abs )
          {
          cur_max_abs  = val1_abs;
          cur_max_orig = val1_orig;
          }
        }

      ++elem_in_col;
      ++it;
      }

    if(elem_in_col == p.get_n_rows())
      {
      result.at(0, cur_col) = cur_max_orig;
      }
    else
      {
      eT val1 = eT(0);
      
      result.at(0, cur_col) = ( std::abs(val1) >= cur_max_abs ) ? val1 : cur_max_orig;
      }
    }
  else
    {
    // maximum in each row
    result.set_size(p.get_n_rows(), 1);
  
    if(p.get_n_nonzero() == 0)
      {
      return;
      }
    
    typename SpProxy<T1>::const_row_iterator_type it = p.begin_row();
    
    uword cur_row = it.row();
    uword elem_in_row = 1;
    
    eT cur_max_orig = *it;
     T cur_max_abs  = std::abs(cur_max_orig);
    
    ++it;
    
    while(it.pos() < p.get_n_nonzero())
      {
      if(it.row() != cur_row)
        {
        // was the row full?
        if(elem_in_row == p.get_n_cols())
          {
          result.at(cur_row, 0) = cur_max_orig;
          }
        else
          {
          eT val1 = eT(0);
          
          result.at(cur_row, 0) = ( std::abs(val1) >= cur_max_abs ) ? val1 : cur_max_orig;
          }

        cur_row = it.row();
        elem_in_row = 0;
        
        cur_max_orig = *it;
        cur_max_abs  = std::abs(cur_max_orig);
        }
      else
        {
        eT val1_orig = *it;
         T val1_abs  = std::abs(val1_orig);
        
        if( val1_abs >= cur_max_abs )
          {
          cur_max_abs  = val1_abs;
          cur_max_orig = val1_orig;
          }
        }

      ++elem_in_row;
      ++it;
      }

    if(elem_in_row == p.get_n_cols())
      {
      result.at(cur_row, 0) = cur_max_orig;
      }
    else
      {
      eT val1 = eT(0);
      
      result.at(cur_row, 0) = ( std::abs(val1) >= cur_max_abs ) ? val1 : cur_max_orig;
      }
    }
  }



template<typename T1>
inline
typename T1::elem_type
spop_max::vector_max
  (
  const T1& x,
  const typename arma_not_cx<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> p(x);
  
  if(p.get_n_nonzero() == 0)
    {
    return eT(0);
    }
  
  if(SpProxy<T1>::must_use_iterator == false)
    {
    // direct access of values
    if(p.get_n_nonzero() == p.get_n_elem())
      {
      return op_max::direct_max(p.get_values(), p.get_n_nonzero());
      }
    else
      {
      return std::max(eT(0), op_max::direct_max(p.get_values(), p.get_n_nonzero()));
      }
    }
  else
    {
    // use iterator
    typename SpProxy<T1>::const_iterator_type it = p.begin();

    eT result = (*it);
    ++it;

    while(it != p.end())
      {
      if((*it) > result)
        {
        result = (*it);
        }

      ++it;
      }

    if(p.get_n_nonzero() == p.get_n_elem())
      {
      return result;
      }
    else
      {
      return std::max(eT(0), result);
      }
    }
  }



template<typename T1>
inline
typename T1::elem_type
spop_max::vector_max
  (
  const T1& x,
  const typename arma_cx_only<typename T1::elem_type>::result* junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type            eT;
  typedef typename get_pod_type<eT>::result  T;

  const SpProxy<T1> p(x);

  if(p.get_n_nonzero() == 0)
    {
    return eT(0);
    }

  if(SpProxy<T1>::must_use_iterator == false)
    {
    // direct access of values
    if(p.get_n_nonzero() == p.get_n_elem())
      {
      return op_max::direct_max(p.get_values(), p.get_n_nonzero());
      }
    else
      {
      const eT val1 = eT(0);
      const eT val2 = op_max::direct_max(p.get_values(), p.get_n_nonzero());
      
      return ( std::abs(val1) >= std::abs(val2) ) ? val1 : val2;
      }
    }
  else
    {
    // use iterator
    typename SpProxy<T1>::const_iterator_type it = p.begin();

    eT best_val_orig = *it;
     T best_val_abs  = std::abs(best_val_orig);
    
    ++it;
    
    while(it != p.end())
      {
      eT val_orig = *it;
       T val_abs  = std::abs(val_orig);
      
      if(val_abs > best_val_abs)
        {
        best_val_abs  = val_abs;
        best_val_orig = val_orig;
        }

      ++it;
      }

    if(p.get_n_nonzero() == p.get_n_elem())
      {
      return best_val_orig;
      }
    else
      {
      const eT val1 = eT(0);
      
      return ( std::abs(val1) >= best_val_abs ) ? val1 : best_val_orig;
      }
    }
  }



//! @}
