// Copyright (C) 2012 Ryan Curtin
// Copyright (C) 2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup spglue_minus
//! @{



template<typename T1, typename T2>
inline
void
spglue_minus::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_minus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> pa(X.A);
  const SpProxy<T2> pb(X.B);
  
  const bool is_alias = pa.is_alias(out) || pb.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_minus::apply_noalias(out, pa, pb);
    }
  else
    {
    SpMat<eT> tmp;
    spglue_minus::apply_noalias(tmp, pa, pb);
    
    out.steal_mem(tmp);
    }
  }



template<typename eT, typename T1, typename T2>
arma_hot
inline
void
spglue_minus::apply_noalias(SpMat<eT>& result, const SpProxy<T1>& pa, const SpProxy<T2>& pb)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(pa.get_n_rows(), pa.get_n_cols(), pb.get_n_rows(), pb.get_n_cols(), "subtraction");

  result.set_size(pa.get_n_rows(), pa.get_n_cols());

  // Resize memory to correct size.
  result.mem_resize(n_unique(pa, pb, op_n_unique_sub()));

  // Now iterate across both matrices.
  typename SpProxy<T1>::const_iterator_type x_it = pa.begin();
  typename SpProxy<T2>::const_iterator_type y_it = pb.begin();

  uword cur_val = 0;
  while((x_it.pos() < pa.get_n_nonzero()) || (y_it.pos() < pb.get_n_nonzero()))
    {
    if(x_it == y_it)
      {
      const typename T1::elem_type val = (*x_it) - (*y_it);
      if (val != 0)
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
      if((x_it.col() < y_it.col()) || ((x_it.col() == y_it.col()) && (x_it.row() < y_it.row()))) // if y is closer to the end
        {
        access::rw(result.values[cur_val]) = (*x_it);
        access::rw(result.row_indices[cur_val]) = x_it.row();
        ++access::rw(result.col_ptrs[x_it.col() + 1]);
        ++cur_val;
        ++x_it;
        }
      else
        {
        access::rw(result.values[cur_val]) = -(*y_it);
        access::rw(result.row_indices[cur_val]) = y_it.row();
        ++access::rw(result.col_ptrs[y_it.col() + 1]);
        ++cur_val;
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



//
//
// spglue_minus2: scalar*(A - B)



template<typename T1, typename T2>
inline
void
spglue_minus2::apply(SpMat<typename T1::elem_type>& out, const SpGlue<T1,T2,spglue_minus2>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const SpProxy<T1> pa(X.A);
  const SpProxy<T2> pb(X.B);
  
  const bool is_alias = pa.is_alias(out) || pb.is_alias(out);
  
  if(is_alias == false)
    {
    spglue_minus::apply_noalias(out, pa, pb);
    }
  else
    {
    SpMat<eT> tmp;
    spglue_minus::apply_noalias(tmp, pa, pb);
    
    out.steal_mem(tmp);
    }
  
  out *= X.aux;
  }



//! @}
