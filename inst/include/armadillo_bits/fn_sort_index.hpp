// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup fn_sort_index
//! @{




template<typename T1, typename T2>
struct arma_sort_index_packet_ascend
  {
  T1 val;
  T2 index;
  };



template<typename T1, typename T2>
struct arma_sort_index_packet_descend
  {
  T1 val;
  T2 index;
  };



template<typename T1, typename T2>
inline
bool
operator< (const arma_sort_index_packet_ascend<T1,T2>& A, const arma_sort_index_packet_ascend<T1,T2>& B)
  {
  return A.val < B.val;
  }



template<typename T1, typename T2>
inline
bool
operator< (const arma_sort_index_packet_descend<T1,T2>& A, const arma_sort_index_packet_descend<T1,T2>& B)
  {
  return A.val > B.val;
  }



template<typename umat_elem_type, typename packet_type, typename eT>
void
inline
sort_index_helper(umat_elem_type* out_mem, std::vector<packet_type>& packet_vec, const eT* in_mem)
  {
  arma_extra_debug_sigprint();
  
  const uword n_elem = packet_vec.size();
  
  for(uword i=0; i<n_elem; ++i)
    {
    packet_vec[i].val   = in_mem[i];
    packet_vec[i].index = i;
    }
  
  std::sort( packet_vec.begin(), packet_vec.end() );
  
  for(uword i=0; i<n_elem; ++i)
    {
    out_mem[i] = packet_vec[i].index;
    }
  }



template<typename T1>
inline
umat
sort_index(const Base<typename T1::elem_type,T1>& X, const uword sort_type = 0)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  arma_type_check(( is_complex<eT>::value == true ));
  
  const unwrap<T1> tmp(X.get_ref());
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
    std::vector< arma_sort_index_packet_ascend<eT,out_elem_type> > packet_vec(A.n_elem);
    
    sort_index_helper(out.memptr(), packet_vec, A.mem);
    }
  else
    {
    std::vector< arma_sort_index_packet_descend<eT,out_elem_type> > packet_vec(A.n_elem);
    
    sort_index_helper(out.memptr(), packet_vec, A.mem);
    }
  
  return out;
  }


//! @}
