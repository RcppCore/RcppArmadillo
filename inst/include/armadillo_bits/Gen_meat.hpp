// Copyright (C) 2011 NICTA (www.nicta.com.au)
// Copyright (C) 2011 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup Gen
//! @{



template<typename eT, typename gen_type>
arma_inline
Gen<eT, gen_type>::Gen(const uword in_n_rows, const uword in_n_cols)
  : n_rows(in_n_rows)
  , n_cols(in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT, typename gen_type>
arma_inline
Gen<eT, gen_type>::~Gen()
  {
  arma_extra_debug_sigprint();
  }



template<typename eT, typename gen_type>
arma_inline
eT
Gen<eT, gen_type>::generate()
  {
       if(is_same_type<gen_type, gen_ones_full>::value == true) { return eT(1);                   }
  else if(is_same_type<gen_type, gen_zeros    >::value == true) { return eT(0);                   }
  else if(is_same_type<gen_type, gen_randu    >::value == true) { return eT(eop_aux_randu<eT>()); }
  else if(is_same_type<gen_type, gen_randn    >::value == true) { return eT(eop_aux_randn<eT>()); }
  else                                                          { return eT();                    }
  }



template<typename eT, typename gen_type>
arma_inline
eT
Gen<eT, gen_type>::operator[](const uword i) const
  {
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    return ((i % n_rows) == (i / n_rows)) ? eT(1) : eT(0);
    }
  else
    {
    return Gen<eT, gen_type>::generate();
    }
  }



template<typename eT, typename gen_type>
arma_inline
eT
Gen<eT, gen_type>::at(const uword row, const uword col) const
  {
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    return (row == col) ? eT(1) : eT(0);
    }
  else
    {
    return Gen<eT, gen_type>::generate();
    }
  }



template<typename eT, typename gen_type>
inline
void
Gen<eT, gen_type>::apply(Mat<eT>& out) const
  {
  arma_extra_debug_sigprint();
  
  // NOTE: we're assuming that the matrix has already been set to the correct size;
  // this is done by either the Mat contructor or operator=()
  
       if(is_same_type<gen_type, gen_ones_diag>::value == true) { out.eye();   }
  else if(is_same_type<gen_type, gen_ones_full>::value == true) { out.ones();  }
  else if(is_same_type<gen_type, gen_zeros    >::value == true) { out.zeros(); }
  else if(is_same_type<gen_type, gen_randu    >::value == true) { out.randu(); }
  else if(is_same_type<gen_type, gen_randn    >::value == true) { out.randn(); }
  }



template<typename eT, typename gen_type>
inline
void
Gen<eT, gen_type>::apply_inplace_plus(Mat<eT>& out) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "addition");
  
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword i=0; i<N; ++i)
      {
      out.at(i,i) += eT(1);
      }
    }
  else
    {
          eT*   out_mem = out.memptr();
    const uword n_elem  = out.n_elem;
    
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = Gen<eT, gen_type>::generate();
      const eT tmp_j = Gen<eT, gen_type>::generate();
      
      out_mem[i] += tmp_i;
      out_mem[j] += tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] += Gen<eT, gen_type>::generate();
      }
    }
  
  }




template<typename eT, typename gen_type>
inline
void
Gen<eT, gen_type>::apply_inplace_minus(Mat<eT>& out) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "subtraction");
  
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword i=0; i<N; ++i)
      {
      out.at(i,i) -= eT(1);
      }
    }
  else
    {
          eT*   out_mem = out.memptr();
    const uword n_elem  = out.n_elem;
    
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = Gen<eT, gen_type>::generate();
      const eT tmp_j = Gen<eT, gen_type>::generate();
      
      out_mem[i] -= tmp_i;
      out_mem[j] -= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] -= Gen<eT, gen_type>::generate();
      }
    }
  
  }




template<typename eT, typename gen_type>
inline
void
Gen<eT, gen_type>::apply_inplace_schur(Mat<eT>& out) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "element-wise multiplication");
  
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword i=0; i<N; ++i)
      {
      for(uword row=0;   row<i;      ++row) { out.at(row,i) = eT(0); }
      for(uword row=i+1; row<n_rows; ++row) { out.at(row,i) = eT(0); }
      }
    }
  else
    {
          eT*   out_mem = out.memptr();
    const uword n_elem  = out.n_elem;
    
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = Gen<eT, gen_type>::generate();
      const eT tmp_j = Gen<eT, gen_type>::generate();
      
      out_mem[i] *= tmp_i;
      out_mem[j] *= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] *= Gen<eT, gen_type>::generate();
      }
    }
  
  }




template<typename eT, typename gen_type>
inline
void
Gen<eT, gen_type>::apply_inplace_div(Mat<eT>& out) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "element-wise division");
  
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword i=0; i<N; ++i)
      {
      const eT zero = eT(0);
      
      for(uword row=0;   row<i;      ++row) { out.at(row,i) /= zero; }
      for(uword row=i+1; row<n_rows; ++row) { out.at(row,i) /= zero; }
      }
    }
  else
    {
          eT*   out_mem = out.memptr();
    const uword n_elem  = out.n_elem;
    
    uword i,j;
    
    for(i=0, j=1; j<n_elem; i+=2, j+=2)
      {
      const eT tmp_i = Gen<eT, gen_type>::generate();
      const eT tmp_j = Gen<eT, gen_type>::generate();
      
      out_mem[i] /= tmp_i;
      out_mem[j] /= tmp_j;
      }
    
    if(i < n_elem)
      {
      out_mem[i] /= Gen<eT, gen_type>::generate();
      }
    }
  
  }




//! @}
