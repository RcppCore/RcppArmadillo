// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// Copyright (C) 2012 Ryan Curtin
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup operator_schur
//! @{


// operator %, which we define it to do a schur product (element-wise multiplication)


//! element-wise multiplication of user-accessible Armadillo objects with same element type
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value,
  const eGlue<T1, T2, eglue_schur>
  >::result
operator%
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return eGlue<T1, T2, eglue_schur>(X, Y);
  }



//! element-wise multiplication of user-accessible Armadillo objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value == false)),
  const mtGlue<typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_schur>
  >::result
operator%
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  return mtGlue<out_eT, T1, T2, glue_mixed_schur>( X, Y );
  }



//! element-wise multiplication of two sparse matrices
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  SpMat<typename T1::elem_type>
  >::result
operator%
  (
  const SpBase<typename T1::elem_type, T1>& x,
  const SpBase<typename T2::elem_type, T2>& y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;

  const SpProxy<T1> pa(x.get_ref());
  const SpProxy<T2> pb(y.get_ref());

  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise multiplication");

  SpMat<typename T1::elem_type> result(pa.get_n_rows(), pa.get_n_cols());
  
  if( (pa.get_n_nonzero() != 0) && (pb.get_n_nonzero() != 0) )
    {
    // Resize memory to correct size.
    result.mem_resize(n_unique(x, y, op_n_unique_mul()));
    
    // Now iterate across both matrices.
    typename SpProxy<T1>::const_iterator_type x_it = pa.begin();
    typename SpProxy<T2>::const_iterator_type y_it = pb.begin();
    
    typename SpProxy<T1>::const_iterator_type x_end = pa.end();
    typename SpProxy<T2>::const_iterator_type y_end = pb.end();
    
    uword cur_val = 0;
    while((x_it != x_end) || (y_it != y_end))
      {
      if(x_it == y_it)
        {
        const eT val = (*x_it) * (*y_it);
        
        if (val != eT(0))
          {
          access::rw(result.values[cur_val]) = val;
          access::rw(result.row_indices[cur_val]) = x_it.row();
          ++access::rw(result.col_ptrs[x_it.col() + 1]);
          ++cur_val;
          }
        
        ++x_it;
        ++y_it;
        }
      else
        {
        const uword x_it_row = x_it.row();
        const uword x_it_col = x_it.col();
        
        const uword y_it_row = y_it.row();
        const uword y_it_col = y_it.col();
        
        if((x_it_col < y_it_col) || ((x_it_col == y_it_col) && (x_it_row < y_it_row))) // if y is closer to the end
          {
          ++x_it;
          }
        else
          {
          ++y_it;
          }
        }
      }
    
    // Fix column pointers to be cumulative.
    for(uword c = 1; c <= result.n_cols; ++c)
      {
      access::rw(result.col_ptrs[c]) += result.col_ptrs[c - 1];
      }
    }
  
  return result;
  }



//! element-wise multiplication of one sparse and one dense object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  SpMat<typename T1::elem_type>
  >::result
operator%
  (
  const SpBase<typename T1::elem_type, T1>& x,
  const   Base<typename T2::elem_type, T2>& y
  )
  {
  // This operation is commutative.
  return (y % x);
  }



//! element-wise multiplication of one dense and one sparse object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  SpMat<typename T1::elem_type>
  >::result
operator%
  (
  const   Base<typename T1::elem_type, T1>& x,
  const SpBase<typename T2::elem_type, T2>& y
  )
  {
  arma_extra_debug_sigprint();

  const   Proxy<T1> pa(x.get_ref());
  const SpProxy<T2> pb(y.get_ref());

  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise multiplication");

  SpMat<typename T1::elem_type> result(pa.get_n_rows(), pa.get_n_cols());

  if(Proxy<T1>::prefer_at_accessor == false)
    {
    // use direct operator[] access
    // count new size
    uword new_n_nonzero = 0;
    typename SpProxy<T2>::const_iterator_type it = pb.begin();
    while(it != pb.end())
      {
      if(((*it) * pa[(it.col() * pa.get_n_rows()) + it.row()]) != 0)
        {
        ++new_n_nonzero;
        }

      ++it;
      }

    // Resize memory accordingly.
    result.mem_resize(new_n_nonzero);

    uword cur_val = 0;
    typename SpProxy<T2>::const_iterator_type it2 = pb.begin();
    while(it2 != pb.end())
      {
      const typename T1::elem_type val = (*it2) * pa[(it2.col() * pa.get_n_rows()) + it2.row()];
      if(val != 0)
        {
        access::rw(result.values[cur_val]) = val;
        access::rw(result.row_indices[cur_val]) = it2.row();
        ++access::rw(result.col_ptrs[it2.col() + 1]);
        ++cur_val;
        }

      ++it2;
      }
    }
  else
    {
    // use at() access
    // count new size
    uword new_n_nonzero = 0;
    typename SpProxy<T2>::const_iterator_type it = pb.begin();
    while(it != pb.end())
      {
      if(((*it) * pa.at(it.row(), it.col())) != 0)
        {
        ++new_n_nonzero;
        }

      ++it;
      }

    // Resize memory accordingly.
    result.mem_resize(new_n_nonzero);

    uword cur_val = 0;
    typename SpProxy<T2>::const_iterator_type it2 = pb.begin();
    while(it2 != pb.end())
      {
      const typename T1::elem_type val = (*it2) * pa.at(it2.row(), it2.col());
      if(val != 0)
        {
        access::rw(result.values[cur_val]) = val;
        access::rw(result.row_indices[cur_val]) = it2.row();
        ++access::rw(result.col_ptrs[it2.col() + 1]);
        ++cur_val;
        }

      ++it2;
      }
    }

  // Fix column pointers.
  for(uword c = 1; c <= result.n_cols; ++c)
    {
    access::rw(result.col_ptrs[c]) += result.col_ptrs[c - 1];
    }

  return result;
  }


//! @}
