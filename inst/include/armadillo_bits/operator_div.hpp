// Copyright (C) 2009-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2012 Conrad Sanderson
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


//! \addtogroup operator_div
//! @{



//! Base / scalar
template<typename T1>
arma_inline
typename
enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_scalar_div_post> >::result
operator/
  (
  const T1&                    X,
  const typename T1::elem_type k
  )
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_div_post>(X, k);
  }



//! scalar / Base
template<typename T1>
arma_inline
typename
enable_if2< is_arma_type<T1>::value, const eOp<T1, eop_scalar_div_pre> >::result
operator/
  (
  const typename T1::elem_type k,
  const T1&                    X
  )
  {
  arma_extra_debug_sigprint();
  
  return eOp<T1, eop_scalar_div_pre>(X, k);
  }



//! complex scalar / non-complex Base
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_complex<typename T1::elem_type>::value == false),
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_pre>
  >::result
operator/
  (
  const std::complex<typename T1::pod_type>& k,
  const T1&                                  X
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_pre>('j', X, k);
  }



//! non-complex Base / complex scalar
template<typename T1>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_complex<typename T1::elem_type>::value == false),
  const mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_post>
  >::result
operator/
  (
  const T1&                                  X,
  const std::complex<typename T1::pod_type>& k
  )
  {
  arma_extra_debug_sigprint();
  
  return mtOp<typename std::complex<typename T1::pod_type>, T1, op_cx_scalar_div_post>('j', X, k);
  }



//! element-wise division of Base objects with same element type
template<typename T1, typename T2>
arma_inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  const eGlue<T1, T2, eglue_div>
  >::result
operator/
  (
  const T1& X,
  const T2& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return eGlue<T1, T2, eglue_div>(X, Y);
  }



//! element-wise division of Base objects with different element types
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_type<T2>::value && (is_same_type<typename T1::elem_type, typename T2::elem_type>::value == false)),
  const mtGlue<typename promote_type<typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_div>
  >::result
operator/
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
  
  return mtGlue<out_eT, T1, T2, glue_mixed_div>( X, Y );
  }



//! element-wise division of sparse matrix by scalar
template<typename T1>
inline
typename
enable_if2<is_arma_sparse_type<T1>::value, SpMat<typename T1::elem_type> >::result
operator/
  (
  const SpBase<typename T1::elem_type, T1>& X,
  const typename T1::elem_type y
  )
  {
  arma_extra_debug_sigprint();

  arma_debug_check(y == typename T1::elem_type(0), "element-wise division: division by zero");

  SpMat<typename T1::elem_type> result(X.get_ref());

  for(uword i = 0; i < result.n_nonzero; ++i)
    {
    access::rw(result.values[i]) /= y;
    }

  return result;
  }



// //! element-wise division of two sparse objects.  what a bad idea
// template<typename T1, typename T2>
// inline
// typename
// enable_if2
//   <
//   (is_arma_sparse_type<T1>::value && is_arma_sparse_type<T2>::value &&
// is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
//   SpMat<typename T1::elem_type>
//   >::result
// operator/
//   (
//   const SpBase<typename T1::elem_type, T1>& x,
//   const SpBase<typename T2::elem_type, T2>& y
//   )
//   {
//   arma_extra_debug_sigprint();
// 
//   const SpProxy<T1> pa(x.get_ref());
//   const SpProxy<T2> pb(y.get_ref());
// 
//   arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise division");
// 
//   SpMat<typename T1::elem_type> result(pa.get_n_rows(), pa.get_n_cols());
// 
//   // terrible
//   for(uword i = 0; i < result.n_elem; ++i)
//     {
//     result[i] = (pa[i] / pb[i]);
//     }
// 
//   return result;
//   }



//! element-wise division of one sparse and one dense object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_sparse_type<T1>::value && is_arma_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  SpMat<typename T1::elem_type>
  >::result
operator/
  (
  const SpBase<typename T1::elem_type, T1>& x,
  const   Base<typename T2::elem_type, T2>& y
  )
  {
  arma_extra_debug_sigprint();

  const SpProxy<T1> pa(x.get_ref());
  const   Proxy<T2> pb(y.get_ref());

  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise division");

  SpMat<typename T1::elem_type> result(pa.get_n_rows(), pa.get_n_cols());

  // The compiler should be smart enough to optimize out the inner if/else statement entirely
  typename SpProxy<T1>::const_iterator_type it = pa.begin();
  uword new_n_nonzero = 0;
  while(it != pa.end())
    {
    if(Proxy<T2>::prefer_at_accessor == false)
      {
      const typename T1::elem_type val = (*it) / pb[(it.col() * pb.get_n_rows()) + it.row()];
      if(val != 0)
        {
        ++new_n_nonzero;
        }
      }
    else
      {
      const typename T1::elem_type val = (*it) / pb.at(it.row(), it.col());
      if(val != 0)
        {
        ++new_n_nonzero;
        }
      }

    ++it;
    }

  result.mem_resize(new_n_nonzero);

  typename SpProxy<T1>::const_iterator_type it2 = pa.begin();
  uword cur_pos = 0;
  while(it2 != pa.end())
    {
    if(Proxy<T2>::prefer_at_accessor == false)
      {
      const typename T1::elem_type val = (*it2) / pb[(it2.col() * pb.get_n_rows()) + it2.row()];
      if(val != 0)
        {
        access::rw(result.values[cur_pos]) = val;
        access::rw(result.row_indices[cur_pos]) = it2.row();
        ++access::rw(result.col_ptrs[it2.col() + 1]);
        ++cur_pos;
        }
      }
    else
      {
      const typename T1::elem_type val = (*it2) / pb.at(it2.row(), it2.col());
      if(val != 0)
        {
        access::rw(result.values[cur_pos]) = val;
        access::rw(result.row_indices[cur_pos]) = it2.row();
        ++access::rw(result.col_ptrs[it2.col() + 1]);
        ++cur_pos;
        }
      }

    ++it2;
    }

  // Fix column pointers
  for(uword col = 1; col <= result.n_cols; ++col)
    {
    access::rw(result.col_ptrs[col]) += result.col_ptrs[col - 1];
    }

  return result;
  }



//! element-wise division of one dense and one sparse object
template<typename T1, typename T2>
inline
typename
enable_if2
  <
  (is_arma_type<T1>::value && is_arma_sparse_type<T2>::value && is_same_type<typename T1::elem_type, typename T2::elem_type>::value),
  Mat<typename T1::elem_type>
  >::result
operator/
  (
  const   Base<typename T1::elem_type, T1>& x,
  const SpBase<typename T2::elem_type, T2>& y
  )
  {
  arma_extra_debug_sigprint();

  const   Proxy<T1> pa(x.get_ref());
  const SpProxy<T2> pb(y.get_ref());

  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "element-wise division");

  Mat<typename T1::elem_type> result(pa.get_n_rows(), pa.get_n_cols());

  result.fill(Datum<typename T1::elem_type>::inf);

  // Now divide each element
  typename SpProxy<T2>::const_iterator_type it = pb.begin();

  while(it != pb.end())
    {
    if(Proxy<T1>::prefer_at_accessor == false)
      {
      const uword index = (it.col() * result.n_rows) + it.row();
      result[index] = pa[index] / (*it);
      }
    else
      {
      result.at(it.row(), it.col()) = pa.at(it.row(), it.col()) / (*it);
      }

    ++it;
    }

  return result;
  }



//! @}
