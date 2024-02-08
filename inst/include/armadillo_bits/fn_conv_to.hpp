// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup fn_conv_to
//! @{



//! conversion from Armadillo Base and BaseCube objects to scalars
//! NOTE: use as_scalar() instead; this functionality is kept only for compatibility with old user code
template<typename out_eT>
class conv_to
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static out_eT from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);

  template<typename in_eT, typename T1>
  inline static out_eT from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static out_eT from(const BaseCube<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static out_eT from(const BaseCube<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
out_eT
conv_to<out_eT>::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_type_check(( is_supported_elem_type<out_eT>::value == false ));
  
  const Proxy<T1> P(in.get_ref());
  
  arma_debug_check( (P.get_n_elem() != 1), "conv_to(): given object does not have exactly one element" );
  
  return out_eT(Proxy<T1>::use_at ? P.at(0,0) : P[0]);
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
out_eT
conv_to<out_eT>::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_type_check(( is_supported_elem_type<out_eT>::value == false ));
  
  const Proxy<T1> P(in.get_ref());
  
  arma_debug_check( (P.get_n_elem() != 1), "conv_to(): given object does not have exactly one element" );
  
  out_eT out;
  
  arrayops::convert_cx_scalar(out, (Proxy<T1>::use_at ? P.at(0,0) : P[0]));
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
out_eT
conv_to<out_eT>::from(const BaseCube<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_type_check(( is_supported_elem_type<out_eT>::value == false ));
  
  const ProxyCube<T1> P(in.get_ref());
  
  arma_debug_check( (P.get_n_elem() != 1), "conv_to(): given object does not have exactly one element" );
  
  return out_eT(ProxyCube<T1>::use_at ? P.at(0,0,0) : P[0]);
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
out_eT
conv_to<out_eT>::from(const BaseCube<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_type_check(( is_supported_elem_type<out_eT>::value == false ));
  
  const ProxyCube<T1> P(in.get_ref());
  
  arma_debug_check( (P.get_n_elem() != 1), "conv_to(): given object does not have exactly one element" );
  
  out_eT out;
  
  arrayops::convert_cx_scalar(out, (ProxyCube<T1>::use_at ? P.at(0,0,0) : P[0]));
  
  return out;
  }



//! conversion to Armadillo matrices from Armadillo Base objects, as well as from std::vector
template<typename out_eT>
class conv_to< Mat<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static Mat<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static Mat<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  
  //
  
  template<typename in_eT, typename T1>
  inline static Mat<out_eT> from(const SpBase<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static Mat<out_eT> from(const SpBase<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  
  //
  
  template<typename in_eT>
  inline static Mat<out_eT> from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT>
  inline static Mat<out_eT> from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  const Mat<in_eT>& X  = tmp.M;
  
  Mat<out_eT> out(X.n_rows, X.n_cols, arma_nozeros_indicator());
  
  arrayops::convert( out.memptr(), X.memptr(), X.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  const Mat<in_eT>& X  = tmp.M;
  
  Mat<out_eT> out(X.n_rows, X.n_cols, arma_nozeros_indicator());
  
  arrayops::convert_cx( out.memptr(), X.memptr(), X.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const SpBase<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap_spmat<T1> U(in.get_ref());
  const SpMat<in_eT>&    X = U.M;
  
  Mat<out_eT> out(X.n_rows, X.n_cols, arma_zeros_indicator());
  
  podarray<out_eT> tmp(X.n_nonzero);
  
  arrayops::convert( tmp.memptr(), X.values, X.n_nonzero );
  
  typename SpMat<in_eT>::const_iterator it     = X.begin();
  typename SpMat<in_eT>::const_iterator it_end = X.end();
  
  for(uword count=0; it != it_end; ++it, ++count)  { out.at(it.row(), it.col()) = tmp[count]; }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const SpBase<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap_spmat<T1> U(in.get_ref());
  const SpMat<in_eT>&    X = U.M;
  
  Mat<out_eT> out(X.n_rows, X.n_cols, arma_zeros_indicator());
  
  podarray<out_eT> tmp(X.n_nonzero);
  
  arrayops::convert_cx( tmp.memptr(), X.values, X.n_nonzero );
  
  typename SpMat<in_eT>::const_iterator it     = X.begin();
  typename SpMat<in_eT>::const_iterator it_end = X.end();
  
  for(uword count=0; it != it_end; ++it, ++count)  { out.at(it.row(), it.col()) = tmp[count]; }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
arma_warn_unused
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Mat<out_eT> out(N, 1, arma_nozeros_indicator());
  
  if(N > 0)
    {
    arrayops::convert( out.memptr(), &(in[0]), N );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
arma_warn_unused
inline
Mat<out_eT>
conv_to< Mat<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Mat<out_eT> out(N, 1, arma_nozeros_indicator());
  
  if(N > 0)
    {
    arrayops::convert_cx( out.memptr(), &(in[0]), N );
    }
  
  return out;
  }



//! conversion to Armadillo row vectors from Armadillo Base objects, as well as from std::vector
template<typename out_eT>
class conv_to< Row<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static Row<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static Row<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  
  //
  
  template<typename in_eT>
  inline static Row<out_eT> from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT>
  inline static Row<out_eT> from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  const Mat<in_eT>& X  = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object cannot be interpreted as a vector" );
  
  Row<out_eT> out(X.n_elem, arma_nozeros_indicator());
  
  arrayops::convert( out.memptr(), X.memptr(), X.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  const Mat<in_eT>& X  = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object cannot be interpreted as a vector" );
  
  Row<out_eT> out(X.n_rows, X.n_cols, arma_nozeros_indicator());
  
  arrayops::convert_cx( out.memptr(), X.memptr(), X.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
arma_warn_unused
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Row<out_eT> out(N, arma_nozeros_indicator());
  
  if(N > 0)
    {
    arrayops::convert( out.memptr(), &(in[0]), N );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
arma_warn_unused
inline
Row<out_eT>
conv_to< Row<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Row<out_eT> out(N, arma_nozeros_indicator());
  
  if(N > 0)
    {
    arrayops::convert_cx( out.memptr(), &(in[0]), N );
    }
  
  return out;
  }



//! conversion to Armadillo column vectors from Armadillo Base objects, as well as from std::vector
template<typename out_eT>
class conv_to< Col<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static Col<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static Col<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  
  //
  
  template<typename in_eT>
  inline static Col<out_eT> from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT>
  inline static Col<out_eT> from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  const Mat<in_eT>& X  = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object cannot be interpreted as a vector" );
  
  Col<out_eT> out(X.n_elem, arma_nozeros_indicator());
  
  arrayops::convert( out.memptr(), X.memptr(), X.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  const Mat<in_eT>& X  = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object cannot be interpreted as a vector" );
  
  Col<out_eT> out(X.n_rows, X.n_cols, arma_nozeros_indicator());
  
  arrayops::convert_cx( out.memptr(), X.memptr(), X.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
arma_warn_unused
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Col<out_eT> out(N, arma_nozeros_indicator());
  
  if(N > 0)
    {
    arrayops::convert( out.memptr(), &(in[0]), N );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
arma_warn_unused
inline
Col<out_eT>
conv_to< Col<out_eT> >::from(const std::vector<in_eT>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const uword N = uword( in.size() );
  
  Col<out_eT> out(N, arma_nozeros_indicator());
  
  if(N > 0)
    {
    arrayops::convert_cx( out.memptr(), &(in[0]), N );
    }
  
  return out;
  }



//! convert between SpMat types
template<typename out_eT>
class conv_to< SpMat<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static SpMat<out_eT> from(const SpBase<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static SpMat<out_eT> from(const SpBase<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  
  //
  
  template<typename in_eT, typename T1>
  inline static SpMat<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static SpMat<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
SpMat<out_eT>
conv_to< SpMat<out_eT> >::from(const SpBase<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap_spmat<T1>  tmp(in.get_ref());
  const SpMat<in_eT>& X = tmp.M;
  
  SpMat<out_eT> out(arma_layout_indicator(), X);
  
  arrayops::convert( access::rwp(out.values), X.values, X.n_nonzero );
  
  out.remove_zeros();
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
SpMat<out_eT>
conv_to< SpMat<out_eT> >::from(const SpBase<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap_spmat<T1>  tmp(in.get_ref());
  const SpMat<in_eT>& X = tmp.M;
  
  SpMat<out_eT> out(arma_layout_indicator(), X);
  
  arrayops::convert_cx( access::rwp(out.values), X.values, X.n_nonzero );
  
  out.remove_zeros();
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
SpMat<out_eT>
conv_to< SpMat<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  SpMat<out_eT> out;
  
  const quasi_unwrap<T1> U(in.get_ref());
  const Mat<in_eT>&  X = U.M;
  
  if(is_same_type<out_eT,in_eT>::yes)
    {
    const Mat<out_eT>& Y = reinterpret_cast<const Mat<out_eT>&>(X);
    
    SpMat<out_eT> tmp(Y);
    
    out.steal_mem(tmp);
    }
  else
    {
    const uword X_n_rows = X.n_rows;
    const uword X_n_cols = X.n_cols;
    const uword X_n_elem = X.n_elem;
  
    const in_eT* X_mem = X.memptr();
    
    uword X_nnz = 0;
    
    for(uword i=0; i < X_n_elem; ++i)  { X_nnz += (X_mem[i] != in_eT(0)) ? uword(1) : uword(0); }
    
    podarray< in_eT> X_nonzeros(X_nnz);
    podarray<out_eT> Y_nonzeros(X_nnz);
    
    for(uword i=0,count=0; i < X_n_elem; ++i)
      {
      const in_eT X_val = X_mem[i];
      
      if(X_val != in_eT(0))  { X_nonzeros[count] = X_val; ++count; }
      }
    
    arrayops::convert( Y_nonzeros.memptr(), X_nonzeros.memptr(), X_nnz );
    
    if(X_nnz == 0)
      {
      out.set_size(X_n_rows, X.n_cols);
      }
    else
      {
      SpMat<out_eT> tmp(arma_reserve_indicator(), X_n_rows, X_n_cols, X_nnz);
      
      uword count = 0;
      
      for(uword c=0; c < X_n_cols; ++c)
      for(uword r=0; r < X_n_rows; ++r)
        {
        const in_eT X_val = (*X_mem);  ++X_mem;
        
        if(X_val != in_eT(0))
          {
          access::rw(tmp.values[count])      = Y_nonzeros[count];
          access::rw(tmp.row_indices[count]) = r;
          access::rw(tmp.col_ptrs[c + 1])++;
          ++count;
          }
        }
      
      // Sum column counts to be column pointers.
      for(uword c=1; c <= tmp.n_cols; ++c)
        {
        access::rw(tmp.col_ptrs[c]) += tmp.col_ptrs[c - 1];
        }
      
      tmp.remove_zeros();  // in case conversion resulted in an element equal to zero
      
      out.steal_mem(tmp);
      }
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
SpMat<out_eT>
conv_to< SpMat<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  SpMat<out_eT> out;
  
  const quasi_unwrap<T1> U(in.get_ref());
  const Mat<in_eT>&  X = U.M;
  
  if(is_same_type<out_eT,in_eT>::yes)
    {
    const Mat<out_eT>& Y = reinterpret_cast<const Mat<out_eT>&>(X);
    
    SpMat<out_eT> tmp(Y);
    
    out.steal_mem(tmp);
    }
  else
    {
    const uword X_n_rows = X.n_rows;
    const uword X_n_cols = X.n_cols;
    const uword X_n_elem = X.n_elem;
  
    const in_eT* X_mem = X.memptr();
    
    uword X_nnz = 0;
    
    for(uword i=0; i < X_n_elem; ++i)  { X_nnz += (X_mem[i] != in_eT(0)) ? uword(1) : uword(0); }
    
    podarray< in_eT> X_nonzeros(X_nnz);
    podarray<out_eT> Y_nonzeros(X_nnz);
    
    for(uword i=0,count=0; i < X_n_elem; ++i)
      {
      const in_eT X_val = X_mem[i];
      
      if(X_val != in_eT(0))  { X_nonzeros[count] = X_val; ++count; }
      }
    
    arrayops::convert_cx( Y_nonzeros.memptr(), X_nonzeros.memptr(), X_nnz );
    
    if(X_nnz == 0)
      {
      out.set_size(X_n_rows, X.n_cols);
      }
    else
      {
      SpMat<out_eT> tmp(arma_reserve_indicator(), X_n_rows, X_n_cols, X_nnz);
      
      uword count = 0;
      
      for(uword c=0; c < X_n_cols; ++c)
      for(uword r=0; r < X_n_rows; ++r)
        {
        const in_eT X_val = (*X_mem);  ++X_mem;
        
        if(X_val != in_eT(0))
          {
          access::rw(tmp.values[count])      = Y_nonzeros[count];
          access::rw(tmp.row_indices[count]) = r;
          access::rw(tmp.col_ptrs[c + 1])++;
          ++count;
          }
        }
      
      // Sum column counts to be column pointers.
      for(uword c=1; c <= tmp.n_cols; ++c)
        {
        access::rw(tmp.col_ptrs[c]) += tmp.col_ptrs[c - 1];
        }
      
      tmp.remove_zeros();  // in case conversion resulted in an element equal to zero
      
      out.steal_mem(tmp);
      }
    }
  
  return out;
  }



//! conversion to Armadillo cubes from Armadillo BaseCube objects
template<typename out_eT>
class conv_to< Cube<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static Cube<out_eT> from(const BaseCube<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static Cube<out_eT> from(const BaseCube<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Cube<out_eT>
conv_to< Cube<out_eT> >::from(const BaseCube<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap_cube<T1>  tmp( in.get_ref() );
  const Cube<in_eT>& X = tmp.M;
  
  Cube<out_eT> out(X.n_rows, X.n_cols, X.n_slices, arma_nozeros_indicator());
  
  arrayops::convert( out.memptr(), X.memptr(), X.n_elem );
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
Cube<out_eT>
conv_to< Cube<out_eT> >::from(const BaseCube<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const unwrap_cube<T1>  tmp( in.get_ref() );
  const Cube<in_eT>& X = tmp.M;
  
  Cube<out_eT> out(X.n_rows, X.n_cols, X.n_slices, arma_nozeros_indicator());
  
  arrayops::convert_cx( out.memptr(), X.memptr(), X.n_elem );
  
  return out;
  }



//! conversion to std::vector from Armadillo Base objects
template<typename out_eT>
class conv_to< std::vector<out_eT> >
  {
  public:
  
  template<typename in_eT, typename T1>
  inline static std::vector<out_eT> from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk = nullptr);
  
  template<typename in_eT, typename T1>
  inline static std::vector<out_eT> from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk = nullptr);
  };



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
std::vector<out_eT>
conv_to< std::vector<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_not_cx<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  const Mat<in_eT>& X  = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object cannot be interpreted as a vector" );
  
  const uword N = X.n_elem;
  
  std::vector<out_eT> out(N);
  
  if(N > 0)
    {
    arrayops::convert( &(out[0]), X.memptr(), N );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
arma_warn_unused
inline
std::vector<out_eT>
conv_to< std::vector<out_eT> >::from(const Base<in_eT, T1>& in, const typename arma_cx_only<in_eT>::result* junk)
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  const quasi_unwrap<T1> tmp(in.get_ref());
  const Mat<in_eT>& X  = tmp.M;
  
  arma_debug_check( ( (X.is_vec() == false) && (X.is_empty() == false) ), "conv_to(): given object cannot be interpreted as a vector" );
  
  const uword N = X.n_elem;
  
  std::vector<out_eT> out(N);
  
  if(N > 0)
    {
    arrayops::convert_cx( &(out[0]), X.memptr(), N );
    }
  
  return out;
  }



//! @}
