// Copyright (C) 2011-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2011-2012 Conrad Sanderson
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



template<typename T1, typename gen_type>
arma_inline
Gen<T1, gen_type>::Gen(const uword in_n_rows, const uword in_n_cols)
  : n_rows(in_n_rows)
  , n_cols(in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename gen_type>
arma_inline
Gen<T1, gen_type>::~Gen()
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename gen_type>
arma_inline
typename T1::elem_type
Gen<T1, gen_type>::generate()
  {
  typedef typename T1::elem_type eT;
  
       if(is_same_type<gen_type, gen_ones_full>::value == true) { return eT(1);                   }
  else if(is_same_type<gen_type, gen_zeros    >::value == true) { return eT(0);                   }
  else if(is_same_type<gen_type, gen_randu    >::value == true) { return eT(eop_aux_randu<eT>()); }
  else if(is_same_type<gen_type, gen_randn    >::value == true) { return eT(eop_aux_randn<eT>()); }
  else                                                          { return eT();                    }
  }



template<typename T1, typename gen_type>
arma_inline
typename T1::elem_type
Gen<T1, gen_type>::operator[](const uword ii) const
  {
  typedef typename T1::elem_type eT;
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    return ((ii % n_rows) == (ii / n_rows)) ? eT(1) : eT(0);
    }
  else
    {
    return Gen<T1, gen_type>::generate();
    }
  }



template<typename T1, typename gen_type>
arma_inline
typename T1::elem_type
Gen<T1, gen_type>::at(const uword row, const uword col) const
  {
  typedef typename T1::elem_type eT;
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    return (row == col) ? eT(1) : eT(0);
    }
  else
    {
    return Gen<T1, gen_type>::generate();
    }
  }



template<typename T1, typename gen_type>
inline
void
Gen<T1, gen_type>::apply(Mat<typename T1::elem_type>& out) const
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



template<typename T1, typename gen_type>
inline
void
Gen<T1, gen_type>::apply_inplace_plus(Mat<typename T1::elem_type>& out) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "addition");
  
  typedef typename T1::elem_type eT;
  
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword iq=0; iq < N; ++iq)
      {
      out.at(iq,iq) += eT(1);
      }
    }
  else
    {
          eT*   out_mem = out.memptr();
    const uword n_elem  = out.n_elem;
    
    uword iq,jq;
    for(iq=0, jq=1; jq < n_elem; iq+=2, jq+=2)
      {
      const eT tmp_i = Gen<T1, gen_type>::generate();
      const eT tmp_j = Gen<T1, gen_type>::generate();
      
      out_mem[iq] += tmp_i;
      out_mem[jq] += tmp_j;
      }
    
    if(iq < n_elem)
      {
      out_mem[iq] += Gen<T1, gen_type>::generate();
      }
    }
  
  }




template<typename T1, typename gen_type>
inline
void
Gen<T1, gen_type>::apply_inplace_minus(Mat<typename T1::elem_type>& out) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "subtraction");
  
  typedef typename T1::elem_type eT;
  
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword iq=0; iq < N; ++iq)
      {
      out.at(iq,iq) -= eT(1);
      }
    }
  else
    {
          eT*   out_mem = out.memptr();
    const uword n_elem  = out.n_elem;
    
    uword iq,jq;
    for(iq=0, jq=1; jq < n_elem; iq+=2, jq+=2)
      {
      const eT tmp_i = Gen<T1, gen_type>::generate();
      const eT tmp_j = Gen<T1, gen_type>::generate();
      
      out_mem[iq] -= tmp_i;
      out_mem[jq] -= tmp_j;
      }
    
    if(iq < n_elem)
      {
      out_mem[iq] -= Gen<T1, gen_type>::generate();
      }
    }
  
  }




template<typename T1, typename gen_type>
inline
void
Gen<T1, gen_type>::apply_inplace_schur(Mat<typename T1::elem_type>& out) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "element-wise multiplication");
  
  typedef typename T1::elem_type eT;
  
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword iq=0; iq < N; ++iq)
      {
      for(uword row=0;    row < iq;     ++row) { out.at(row,iq) = eT(0); }
      for(uword row=iq+1; row < n_rows; ++row) { out.at(row,iq) = eT(0); }
      }
    }
  else
    {
          eT*   out_mem = out.memptr();
    const uword n_elem  = out.n_elem;
    
    uword iq,jq;
    for(iq=0, jq=1; jq < n_elem; iq+=2, jq+=2)
      {
      const eT tmp_i = Gen<T1, gen_type>::generate();
      const eT tmp_j = Gen<T1, gen_type>::generate();
      
      out_mem[iq] *= tmp_i;
      out_mem[jq] *= tmp_j;
      }
    
    if(iq < n_elem)
      {
      out_mem[iq] *= Gen<T1, gen_type>::generate();
      }
    }
  
  }




template<typename T1, typename gen_type>
inline
void
Gen<T1, gen_type>::apply_inplace_div(Mat<typename T1::elem_type>& out) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "element-wise division");
  
  typedef typename T1::elem_type eT;
  
  
  if(is_same_type<gen_type, gen_ones_diag>::value == true)
    {
    const uword N = (std::min)(n_rows, n_cols);
    
    for(uword iq=0; iq < N; ++iq)
      {
      const eT zero = eT(0);
      
      for(uword row=0;    row < iq;     ++row) { out.at(row,iq) /= zero; }
      for(uword row=iq+1; row < n_rows; ++row) { out.at(row,iq) /= zero; }
      }
    }
  else
    {
          eT*   out_mem = out.memptr();
    const uword n_elem  = out.n_elem;
    
    uword iq,jq;
    for(iq=0, jq=1; jq < n_elem; iq+=2, jq+=2)
      {
      const eT tmp_i = Gen<T1, gen_type>::generate();
      const eT tmp_j = Gen<T1, gen_type>::generate();
      
      out_mem[iq] /= tmp_i;
      out_mem[jq] /= tmp_j;
      }
    
    if(iq < n_elem)
      {
      out_mem[iq] /= Gen<T1, gen_type>::generate();
      }
    }
  
  }




//! @}

