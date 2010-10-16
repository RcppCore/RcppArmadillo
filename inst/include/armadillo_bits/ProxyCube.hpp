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


//! \addtogroup ProxyCube
//! @{



template<typename T1>
class ProxyCube
  {
  public:
  inline ProxyCube(const T1& A)
    {
    arma_type_check< is_arma_type<T1>::value == false >::apply();
    }
  };



// ea_type is the "element accessor" type,
// which can provide access to elements via operator[]

template<typename eT>
class ProxyCube< Cube<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<eT>                                 stored_type;
  typedef const eT*                                ea_type;
  
  arma_aligned const Cube<eT>& Q;
  
  inline explicit ProxyCube(const Cube<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows()       const { return Q.n_rows;       }
  arma_inline u32 get_n_cols()       const { return Q.n_cols;       }
  arma_inline u32 get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline u32 get_n_slices()     const { return Q.n_slices;     }
  arma_inline u32 get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename T1, typename op_type>
class ProxyCube< OpCube<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<elem_type>                          stored_type;
  typedef const elem_type*                         ea_type;
  
  arma_aligned const Cube<elem_type> Q;
  
  inline explicit ProxyCube(const OpCube<T1, op_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows()       const { return Q.n_rows;       }
  arma_inline u32 get_n_cols()       const { return Q.n_cols;       }
  arma_inline u32 get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline u32 get_n_slices()     const { return Q.n_slices;     }
  arma_inline u32 get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename T1, typename T2, typename glue_type>
class ProxyCube< GlueCube<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Cube<elem_type>                          stored_type;
  typedef const elem_type*                         ea_type;
  
  arma_aligned const Cube<elem_type> Q;
  
  inline explicit ProxyCube(const GlueCube<T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }

  arma_inline u32 get_n_rows()       const { return Q.n_rows;       }
  arma_inline u32 get_n_cols()       const { return Q.n_cols;       }
  arma_inline u32 get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline u32 get_n_slices()     const { return Q.n_slices;     }
  arma_inline u32 get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename eT>
class ProxyCube< subview_cube<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_cube<eT>                         stored_type;
  typedef const subview_cube<eT>&                  ea_type;
  
  arma_aligned const subview_cube<eT>& Q;
  
  inline explicit ProxyCube(const subview_cube<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows()       const { return Q.n_rows;       }
  arma_inline u32 get_n_cols()       const { return Q.n_cols;       }
  arma_inline u32 get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline u32 get_n_slices()     const { return Q.n_slices;     }
  arma_inline u32 get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  
  arma_inline ea_type get_ea() const { return Q; }
  };



template<typename T1, typename eop_type>
class ProxyCube< eOpCube<T1, eop_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eOpCube<T1, eop_type>                    stored_type;
  typedef const eOpCube<T1, eop_type>&             ea_type;
  
  arma_aligned const eOpCube<T1, eop_type>& Q;
  
  inline explicit ProxyCube(const eOpCube<T1, eop_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows()       const { return Q.get_n_rows();       }
  arma_inline u32 get_n_cols()       const { return Q.get_n_cols();       }
  arma_inline u32 get_n_elem_slice() const { return Q.get_n_elem_slice(); }
  arma_inline u32 get_n_slices()     const { return Q.get_n_slices();     }
  arma_inline u32 get_n_elem()       const { return Q.get_n_elem();       }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  
  arma_inline ea_type get_ea() const { return Q; }
  };



template<typename T1, typename T2, typename eglue_type>
class ProxyCube< eGlueCube<T1, T2, eglue_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eGlueCube<T1, T2, eglue_type>            stored_type;
  typedef const eGlueCube<T1, T2, eglue_type>&     ea_type;
  
  arma_aligned const eGlueCube<T1, T2, eglue_type>& Q;
  
  inline explicit ProxyCube(const eGlueCube<T1, T2, eglue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows()       const { return Q.get_n_rows();       }
  arma_inline u32 get_n_cols()       const { return Q.get_n_cols();       }
  arma_inline u32 get_n_elem_slice() const { return Q.get_n_elem_slice(); }
  arma_inline u32 get_n_slices()     const { return Q.get_n_slices();     }
  arma_inline u32 get_n_elem()       const { return Q.get_n_elem();       }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  
  arma_inline ea_type get_ea() const { return Q; }
  };



template<typename out_eT, typename T1, typename op_type>
class ProxyCube< mtOpCube<out_eT, T1, op_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Cube<out_eT>                 stored_type;
  typedef          const elem_type*             ea_type;
  
  arma_aligned const Cube<out_eT> Q;
  
  inline explicit ProxyCube(const mtOpCube<out_eT, T1, op_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows()       const { return Q.n_rows;       }
  arma_inline u32 get_n_cols()       const { return Q.n_cols;       }
  arma_inline u32 get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline u32 get_n_slices()     const { return Q.n_slices;     }
  arma_inline u32 get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



template<typename out_eT, typename T1, typename T2, typename glue_type>
class ProxyCube< mtGlueCube<out_eT, T1, T2, glue_type > >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Cube<out_eT>                 stored_type;
  typedef          const elem_type*             ea_type;
  
  arma_aligned const Cube<out_eT> Q;
  
  inline explicit ProxyCube(const mtGlueCube<out_eT, T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline u32 get_n_rows()       const { return Q.n_rows;       }
  arma_inline u32 get_n_cols()       const { return Q.n_cols;       }
  arma_inline u32 get_n_elem_slice() const { return Q.n_elem_slice; }
  arma_inline u32 get_n_slices()     const { return Q.n_slices;     }
  arma_inline u32 get_n_elem()       const { return Q.n_elem;       }
  
  arma_inline elem_type operator[] (const u32 i)                                   const { return Q[i];                  }
  arma_inline elem_type at         (const u32 row, const u32 col, const u32 slice) const { return Q.at(row, col, slice); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  };



//! @}
