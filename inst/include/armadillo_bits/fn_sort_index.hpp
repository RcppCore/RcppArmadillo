// Copyright (C) 2009-2013 Conrad Sanderson
// Copyright (C) 2009-2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup fn_sort_index
//! @{




template<typename T1, typename T2>
struct arma_sort_index_packet
  {
  T1 val;
  T2 index;
  };



class arma_sort_index_helper_ascend
  {
  public:
  
  template<typename T1, typename T2>
  arma_inline
  bool
  operator() (const arma_sort_index_packet<T1,T2>& A, const arma_sort_index_packet<T1,T2>& B) const
    {
    return (A.val < B.val);
    }
  };



class arma_sort_index_helper_descend
  {
  public:
  
  template<typename T1, typename T2>
  arma_inline
  bool
  operator() (const arma_sort_index_packet<T1,T2>& A, const arma_sort_index_packet<T1,T2>& B) const
    {
    return (A.val > B.val);
    }
  };



template<typename umat_elem_type, typename eT, const uword sort_type, const uword sort_stable>
void
inline
sort_index_helper(umat_elem_type* out_mem, const eT* in_mem, const uword n_elem)
  {
  arma_extra_debug_sigprint();
  
  std::vector< arma_sort_index_packet<eT, umat_elem_type> > packet_vec(n_elem);
  
  for(uword i=0; i<n_elem; ++i)
    {
    packet_vec[i].val   = in_mem[i];
    packet_vec[i].index = i;
    }
  
  
  if(sort_type == 0)
    {
    // ascend
    
    arma_sort_index_helper_ascend comparator;
    
    if(sort_stable == 0)
      {
      std::sort( packet_vec.begin(), packet_vec.end(), comparator );
      }
    else
      {
      std::stable_sort( packet_vec.begin(), packet_vec.end(), comparator );
      }
    }
  else
    {
    // descend
    
    arma_sort_index_helper_descend comparator;
    
    if(sort_stable == 0)
      {
      std::sort( packet_vec.begin(), packet_vec.end(), comparator );
      }
    else
      {
      std::stable_sort( packet_vec.begin(), packet_vec.end(), comparator );
      }
    }
  
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = packet_vec[i].index;
    }
  }



template<typename T1>
inline
umat
sort_index
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword sort_type = 0,
  const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  if(A.is_empty() == true)
    {
    return umat();
    }
  
  arma_debug_check( (A.is_vec() == false), "sort_index(): currently only handles vectors");
  
  typedef typename umat::elem_type out_elem_type;
  
  umat out(A.n_rows, A.n_cols);
  
  if(sort_type == 0)
    {
    sort_index_helper<out_elem_type, eT, 0, 0>(out.memptr(), A.mem, A.n_elem);
    }
  else
    {
    sort_index_helper<out_elem_type, eT, 1, 0>(out.memptr(), A.mem, A.n_elem);
    }
  
  return out;
  }



template<typename T1>
inline
umat
stable_sort_index
  (
  const Base<typename T1::elem_type,T1>& X,
  const uword sort_type = 0,
  const typename arma_not_cx<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  if(A.is_empty() == true)
    {
    return umat();
    }
  
  arma_debug_check( (A.is_vec() == false), "stable_sort_index(): currently only handles vectors");
  
  typedef typename umat::elem_type out_elem_type;
  
  umat out(A.n_rows, A.n_cols);
  
  if(sort_type == 0)
    {
    sort_index_helper<out_elem_type, eT, 0, 1>(out.memptr(), A.mem, A.n_elem);
    }
  else
    {
    sort_index_helper<out_elem_type, eT, 1, 1>(out.memptr(), A.mem, A.n_elem);
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  ( (is_arma_type<T1>::value == true) && (is_same_type<T2, char>::value == true) && (is_cx<typename T1::elem_type>::value == false) ),
  umat
  >::result
sort_index
  (
  const T1& X,
  const T2* sort_direction
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X);
  const Mat<eT>& A = tmp.M;
  
  if(A.is_empty() == true)
    {
    return umat();
    }
  
  arma_debug_check( (A.is_vec() == false), "sort_index(): currently only handles vectors");
  
  const char sig = sort_direction[0];
  
  arma_debug_check( ((sig != 'a') && (sig != 'd')), "sort_index(): unknown sort direction" );
  
  typedef typename umat::elem_type out_elem_type;
  
  umat out(A.n_rows, A.n_cols);
  
  if(sig == 'a')
    {
    sort_index_helper<out_elem_type, eT, 0, 0>(out.memptr(), A.mem, A.n_elem);
    }
  else
    {
    sort_index_helper<out_elem_type, eT, 1, 0>(out.memptr(), A.mem, A.n_elem);
    }
  
  return out;
  }



template<typename T1, typename T2>
inline
typename
enable_if2
  <
  ( (is_arma_type<T1>::value == true) && (is_same_type<T2, char>::value == true) && (is_cx<typename T1::elem_type>::value == false) ),
  umat
  >::result
stable_sort_index
  (
  const T1& X,
  const T2* sort_direction
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X);
  const Mat<eT>& A = tmp.M;
  
  if(A.is_empty() == true)
    {
    return umat();
    }
  
  arma_debug_check( (A.is_vec() == false), "stable_sort_index(): currently only handles vectors");
  
  const char sig = sort_direction[0];
  
  arma_debug_check( ((sig != 'a') && (sig != 'd')), "stable_sort_index(): unknown sort direction" );
  
  typedef typename umat::elem_type out_elem_type;
  
  umat out(A.n_rows, A.n_cols);
  
  if(sig == 'a')
    {
    sort_index_helper<out_elem_type, eT, 0, 1>(out.memptr(), A.mem, A.n_elem);
    }
  else
    {
    sort_index_helper<out_elem_type, eT, 1, 1>(out.memptr(), A.mem, A.n_elem);
    }
  
  return out;
  }



//! @}
