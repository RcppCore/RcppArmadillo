// Copyright (C) 2010 NICTA (www.nicta.com.au)
// Copyright (C) 2010 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup eglue_core
//! @{



class eglue_plus : public eglue_core<eglue_plus>
  {
  public:
  
  inline static const char* text() { return "addition"; }
  };



class eglue_minus : public eglue_core<eglue_minus>
  {
  public:
  
  inline static const char* text() { return "subtraction"; }
  };



class eglue_div : public eglue_core<eglue_div>
  {
  public:
  
  inline static const char* text() { return "element-wise division"; }
  };



class eglue_schur : public eglue_core<eglue_schur>
  {
  public:
  
  inline static const char* text() { return "element-wise multiplication"; }
  };



//
// matrices



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(x.get_n_rows(), x.get_n_cols());
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = out.n_elem;
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] = tmp_i;\
      out_mem[j] = tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] = P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else
    {
    arma_stop("eglue_core::apply_proxy(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_plus(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, x.P1, "addition");
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = out.n_elem;
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] += tmp_i;\
      out_mem[j] += tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] += P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else
    {
    arma_stop("eglue_core::apply_inplace_plus(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_minus(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, x.P1, "subtraction");
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = out.n_elem;
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] -= tmp_i;\
      out_mem[j] -= tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] -= P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else 
    {
    arma_stop("eglue_core::apply_inplace_minus(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_schur(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, x.P1, "element-wise multiplication");
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = out.n_elem;
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] *= tmp_i;\
      out_mem[j] *= tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] *= P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else
    {
    arma_stop("eglue_core::apply_inplace_schur(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_div(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, x.P1, "element-wise division");
  
  typedef typename T1::elem_type      eT;
  typedef typename Proxy<T1>::ea_type ea_type1;
  typedef typename Proxy<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = out.n_elem;
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] /= tmp_i;\
      out_mem[j] /= tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] /= P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else
    {
    arma_stop("eglue_core::apply_inplace_div(): unhandled eglue_type");
    }
  }



//
// cubes



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(x.get_n_rows(), x.get_n_cols(), x.get_n_slices());
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type1;
  typedef typename ProxyCube<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
        eT* out_mem = out.memptr();
  const u32 n_elem  = out.n_elem;
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] = tmp_i;\
      out_mem[j] = tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] = P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else
    {
    arma_stop("eglue_core::apply(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_plus(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, x.get_n_rows(), x.get_n_cols(), x.get_n_slices(), "addition");
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type1;
  typedef typename ProxyCube<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
  const u32 n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] += tmp_i;\
      out_mem[j] += tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] += P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else
    {
    arma_stop("eglue_core::apply_inplace_plus(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_minus(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, x.get_n_rows(), x.get_n_cols(), x.get_n_slices(), "subtraction");
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type1;
  typedef typename ProxyCube<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
  const u32 n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] -= tmp_i;\
      out_mem[j] -= tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] -= P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else 
    {
    arma_stop("eglue_core::apply_inplace_minus(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_schur(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, x.get_n_rows(), x.get_n_cols(), x.get_n_slices(), "element-wise multiplication");
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type1;
  typedef typename ProxyCube<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
  const u32 n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] *= tmp_i;\
      out_mem[j] *= tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] *= P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else
    {
    arma_stop("eglue_core::apply_inplace_schur(): unhandled eglue_type");
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_div(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, x.get_n_rows(), x.get_n_cols(), x.get_n_slices(), "element-wise division");
  
  typedef typename T1::elem_type          eT;
  typedef typename ProxyCube<T1>::ea_type ea_type1;
  typedef typename ProxyCube<T2>::ea_type ea_type2;
  
  ea_type1 P1 = x.P1.get_ea();
  ea_type2 P2 = x.P2.get_ea();
  
  const u32 n_elem  = out.n_elem;
        eT* out_mem = out.memptr();
  
  #undef  arma_applier
  #define arma_applier(operator) \
    {\
    u32 i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operator##= P2[i];\
      tmp_j operator##= P2[j];\
      \
      out_mem[i] /= tmp_i;\
      out_mem[j] /= tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] /= P1[i] operator P2[i];\
      }\
    }
  
       if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier(+); }
  else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier(-); }
  else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier(/); }
  else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier(*); }
  else
    {
    arma_stop("eglue_core::apply_inplace_div(): unhandled eglue_type");
    }
  }



//! @}
