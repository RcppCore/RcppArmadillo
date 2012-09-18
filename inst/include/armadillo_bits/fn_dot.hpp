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


//! \addtogroup fn_dot
//! @{


template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value,
  typename T1::elem_type
  >::result
dot
  (
  const T1& A,
  const T2& B
  )
  {
  arma_extra_debug_sigprint();
  
  return op_dot::apply(A,B);
  }



template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value,
  typename T1::elem_type
  >::result
norm_dot
  (
  const T1& A, 
  const T2& B
  )
  {
  arma_extra_debug_sigprint();
  
  return op_norm_dot::apply(A,B);
  }



//
// cdot



template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value && is_not_complex<typename T1::elem_type>::value,
  typename T1::elem_type
  >::result
cdot
  (
  const T1& A,
  const T2& B
  )
  {
  arma_extra_debug_sigprint();
  
  return op_dot::apply(A,B);
  }




template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename
enable_if2
  <
  is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value && is_complex<typename T1::elem_type>::value,
  typename T1::elem_type
  >::result
cdot
  (
  const T1& A,
  const T2& B
  )
  {
  arma_extra_debug_sigprint();
  
  return op_cdot::apply(A,B);
  }



// convert dot(htrans(x), y) to cdot(x,y)

template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename
enable_if2
  <
  is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value && is_complex<typename T1::elem_type>::value,
  typename T1::elem_type
  >::result
dot
  (
  const Op<T1, op_htrans>& A,
  const T2&                B
  )
  {
  arma_extra_debug_sigprint();
  
  return cdot(A.m, B);
  }



//! dot product of two sparse objects
template<typename T1, typename T2>
inline
arma_warn_unused
typename
enable_if2
  <(is_arma_sparse_type<T1>::value) && (is_arma_sparse_type<T2>::value) && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
   typename T1::elem_type
  >::result
dot
  (
  const SpBase<typename T1::elem_type, T1>& x,
  const SpBase<typename T2::elem_type, T2>& y
  )
  {
  arma_extra_debug_sigprint();

  const SpProxy<T1> pa(x.get_ref());
  const SpProxy<T2> pb(y.get_ref());

  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "dot()");

  typedef typename T1::elem_type eT;

  if((&(x.get_ref()) == &(y.get_ref())) && (SpProxy<T1>::must_use_iterator == false))
    {
    // We can do it directly!
    return op_dot::direct_dot_arma(pa.get_n_nonzero(), pa.get_values(), pa.get_values());
    }
  else
    {
    // Iterate over both objects and see when they are the same
    eT result = eT(0);

    typename SpProxy<T1>::const_iterator_type a_it = pa.begin();
    typename SpProxy<T2>::const_iterator_type b_it = pb.begin();

    while((a_it.pos() < pa.get_n_nonzero()) && (b_it.pos() < pb.get_n_nonzero()))
      {
      if(a_it == b_it)
        {
        result += (*a_it) * (*b_it);

        ++a_it;
        ++b_it;
        }
      else if((a_it.col() < b_it.col()) || ((a_it.col() == b_it.col()) && (a_it.row() < b_it.row())))
        {
        // a_it is "behind"
        ++a_it;
        }
      else
        {
        // b_it is "behind"
        ++b_it;
        }
      }

    return result;
    }
  }



//! dot product of one sparse and one dense object
template<typename T1, typename T2>
arma_inline
arma_warn_unused
typename
enable_if2
  <(is_arma_sparse_type<T1>::value) && (is_arma_type<T2>::value) && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
   typename T1::elem_type
  >::result
dot
  (
  const SpBase<typename T1::elem_type, T1>& x,
  const   Base<typename T2::elem_type, T2>& y
  )
  {
  // this is commutative
  return dot(y, x);
  }



//! dot product of one dense and one sparse object
template<typename T1, typename T2>
inline
arma_warn_unused
typename
enable_if2
  <(is_arma_type<T1>::value) && (is_arma_sparse_type<T2>::value) && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
   typename T1::elem_type
  >::result
dot
  (
  const   Base<typename T1::elem_type, T1>& x,
  const SpBase<typename T2::elem_type, T2>& y
  )
  {
  arma_extra_debug_sigprint();

  const Proxy<T1> pa(x.get_ref());
  const SpProxy<T2> pb(y.get_ref());

  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "dot()");

  typedef typename T1::elem_type eT;

  eT result = eT(0);

  typename SpProxy<T2>::const_iterator_type it = pb.begin();

  // prefer_at_accessor won't save us operations
  while(it.pos() < pb.get_n_nonzero())
    {
    result += (*it) * pa.at(it.row(), it.col());
    ++it;
    }

  return result;
  }


//! @}
