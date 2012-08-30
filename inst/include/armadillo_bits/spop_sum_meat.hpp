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


//! \addtogroup spop_sum
//! @{



template<typename T1>
arma_hot
inline
void
spop_sum::apply(SpMat<typename T1::elem_type>& out, const SpOp<T1,spop_sum>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const uword dim = in.aux_uword_a;
  
  const SpProxy<T1> p(in.m);
  
  if(p.is_alias(out) == false)
    {
    if(dim == 0) // find the sum in each column
      {
      out.zeros(1, p.get_n_cols());
      
      typename SpProxy<T1>::const_iterator_type it     = p.begin();
      
      while(it.pos() < p.get_n_nonzero())
        {
        out.at(0, it.col()) += (*it);
        ++it;
        }
      }
    else // find the sum in each row
      {
      out.zeros(p.get_n_rows(), 1);
      
      typename SpProxy<T1>::const_iterator_type it     = p.begin();
      
      while(it.pos() < p.get_n_nonzero())
        {
        out.at(it.row(), 0) += (*it);
        ++it;
        }
      }
    }
  else
    {
    if(dim == 0) // find the sum in each column
      {
      SpMat<eT> result(1, p.get_n_cols());
      
      typename SpProxy<T1>::const_iterator_type it     = p.begin();
      
      while(it.pos() < p.get_n_nonzero())
        {
        result.at(0, it.col()) += (*it);
        ++it;
        }
      
      out.steal_mem(result);
      }
    else // find the sum in each row
      {
      SpMat<eT> result(p.get_n_rows(), 1);
      
      typename SpProxy<T1>::const_iterator_type it     = p.begin();
      
      while(it.pos() < p.get_n_nonzero())
        {
        result.at(it.row(), 0) += (*it);
        ++it;
        }
      
      out.steal_mem(result);
      }
    }
  }




//! @}
