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


//! \addtogroup fn_n_unique
//! @{


//! \brief
//! Get the number of unique nonzero elements in two sparse matrices.
//! This is very useful for determining the amount of memory necessary before
//! a sparse matrix operation on two matrices.

template<typename T1, typename T2, typename op_n_unique_type>
inline
uword
n_unique
  (
  const SpBase<typename T1::elem_type, T1>& x,
  const SpBase<typename T2::elem_type, T2>& y,
  const op_n_unique_type junk
  )
  {
  arma_extra_debug_sigprint();
  
  const SpProxy<T1> pa(x.get_ref());
  const SpProxy<T2> pb(y.get_ref());
  
  return n_unique(pa,pb,junk);
  }



template<typename T1, typename T2, typename op_n_unique_type>
arma_hot
inline
uword
n_unique
  (
  const SpProxy<T1>& pa,
  const SpProxy<T2>& pb,
  const op_n_unique_type junk
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  // Use iterators.
  typename SpProxy<T1>::const_iterator_type x_it = pa.begin();
  typename SpProxy<T2>::const_iterator_type y_it = pb.begin();

  uword total_n_nonzero = 0;

  while((x_it != pa.end()) || (y_it != pb.end()))
    {
    if(x_it == y_it)
      {
      if(op_n_unique_type::eval((*x_it), (*y_it)) != typename T1::elem_type(0))
        {
        ++total_n_nonzero;
        }

      ++x_it;
      ++y_it;
      }
    else
      {
      if((x_it.col() < y_it.col()) || ((x_it.col() == y_it.col()) && (x_it.row() < y_it.row()))) // if y is closer to the end
        {
        if(op_n_unique_type::eval((*x_it), typename T1::elem_type(0)) != typename T1::elem_type(0))
          {
          ++total_n_nonzero;
          }

        ++x_it;
        }
      else // x is closer to the end
        {
        if(op_n_unique_type::eval(typename T1::elem_type(0), (*y_it)) != typename T1::elem_type(0))
          {
          ++total_n_nonzero;
          }

        ++y_it;
        }
      }
    }

  return total_n_nonzero;
  }


// Simple operators.
struct op_n_unique_add
  {
  template<typename eT> inline static eT eval(const eT& l, const eT& r) { return (l + r); }
  };

struct op_n_unique_sub
  {
  template<typename eT> inline static eT eval(const eT& l, const eT& r) { return (l - r); }
  };

struct op_n_unique_mul
  {
  template<typename eT> inline static eT eval(const eT& l, const eT& r) { return (l * r); }
  };

struct op_n_unique_count
  {
  template<typename eT> inline static eT eval(const eT& l, const eT& r) { return 1; }
  };



//! @}
