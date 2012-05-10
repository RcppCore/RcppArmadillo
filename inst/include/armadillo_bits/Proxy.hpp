// Copyright (C) 2010-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup Proxy
//! @{


// ea_type is the "element accessor" type,
// which can provide access to elements via operator[]



template<typename T1>
struct Proxy_default
  {
  inline Proxy_default(const T1&)
    {
    arma_type_check(( is_arma_type<T1>::value == false ));
    }
  };



template<typename T1>
struct Proxy_fixed
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef T1                                       stored_type;
  typedef const elem_type*                         ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = true;
  static const bool fake_mat           = false;
  
  static const bool is_row = T1::is_row;
  static const bool is_col = T1::is_col;
  
  arma_aligned const T1& Q;
  
  inline explicit Proxy_fixed(const T1& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline static uword get_n_rows() { return T1::n_rows; }
  arma_inline static uword get_n_cols() { return T1::n_cols; }
  arma_inline static uword get_n_elem() { return T1::n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&Q) == void_ptr(&X)); }
  };



template<typename T1, bool condition>
struct Proxy_redirect {};

template<typename T1>
struct Proxy_redirect<T1, false> { typedef Proxy_default<T1> result; };

template<typename T1>
struct Proxy_redirect<T1, true>  { typedef Proxy_fixed<T1>   result; };



template<typename T1>
class Proxy : public Proxy_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  public:
  inline Proxy(const T1& A)
    : Proxy_redirect< T1, is_Mat_fixed<T1>::value >::result(A)
    {
    }
  };



template<typename eT>
class Proxy< Mat<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = false;
  static const bool is_col = false;
  
  arma_aligned const Mat<eT>& Q;
  
  inline explicit Proxy(const Mat<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&Q) == void_ptr(&X)); }
  };



template<typename eT>
class Proxy< Col<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Col<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = false;
  static const bool is_col = true;
  
  arma_aligned const Col<eT>& Q;
  
  inline explicit Proxy(const Col<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];   }
  arma_inline elem_type at         (const uword row, const uword) const { return Q[row]; }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&Q) == void_ptr(&X)); }
  };



template<typename eT>
class Proxy< Row<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Row<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = true;
  static const bool is_col = false;
  
  arma_aligned const Row<eT>& Q;
  
  inline explicit Proxy(const Row<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return 1;        }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];   }
  arma_inline elem_type at         (const uword, const uword col) const { return Q[col]; }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&Q) == void_ptr(&X)); }
  };



template<typename T1, typename gen_type>
class Proxy< Gen<T1, gen_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Gen<T1, gen_type>                        stored_type;
  typedef const Gen<T1, gen_type>&                 ea_type;
  
  static const bool prefer_at_accessor = Gen<T1, gen_type>::prefer_at_accessor;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = Gen<T1, gen_type>::is_row;
  static const bool is_col = Gen<T1, gen_type>::is_col;
  
  arma_aligned const Gen<T1, gen_type>& Q;
  
  inline explicit Proxy(const Gen<T1, gen_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return (is_row ? 1 : Q.n_rows);                           }
  arma_inline uword get_n_cols() const { return (is_col ? 1 : Q.n_cols);                           }
  arma_inline uword get_n_elem() const { return (is_row ? 1 : Q.n_rows) * (is_col ? 1 : Q.n_cols); }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>&) const { return false; }
  };



template<typename T1, typename op_type>
class Proxy< Op<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = Op<T1, op_type>::is_row;
  static const bool is_col = Op<T1, op_type>::is_col;
  
  arma_aligned const Mat<elem_type> Q;
  
  inline explicit Proxy(const Op<T1, op_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>&) const { return false; }
  };



template<typename T1>
struct Proxy_xtrans_default
  {
  typedef typename T1::elem_type eT;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  arma_aligned const Mat<eT> Q;
  
  arma_hot
  inline Proxy_xtrans_default(const T1& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>&) const { return false; }
  };



template<typename T1>
struct Proxy_xtrans_vector
  {
  inline Proxy_xtrans_vector(const T1&) {}
  };



template<typename T1>
struct Proxy_xtrans_vector< Op<T1, op_htrans> >
  {
  typedef typename T1::elem_type eT;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = quasi_unwrap<T1>::has_subview;
  static const bool is_fixed           = false;
  static const bool fake_mat           = true;
  
  arma_aligned const quasi_unwrap<T1> U; // avoid copy if T1 is a Row, Col or subview_col
  arma_aligned const Mat<eT>          Q;
  
  inline Proxy_xtrans_vector(const Op<T1, op_htrans>& A)
    : U(A.m)
    , Q(const_cast<eT*>(U.M.memptr()), U.M.n_cols, U.M.n_rows, false, false)
    {
    arma_extra_debug_sigprint();
    }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return U.is_alias(X); }
  };



template<typename T1>
struct Proxy_xtrans_vector< Op<T1, op_strans> >
  {
  typedef typename T1::elem_type eT;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = quasi_unwrap<T1>::has_subview;
  static const bool is_fixed           = false;
  static const bool fake_mat           = true;
  
  arma_aligned const quasi_unwrap<T1> U; // avoid copy if T1 is a Row, Col or subview_col
  arma_aligned const Mat<eT>          Q;
  
  inline Proxy_xtrans_vector(const Op<T1, op_strans>& A)
    : U(A.m)
    , Q(const_cast<eT*>(U.M.memptr()), U.M.n_cols, U.M.n_rows, false, false)
    {
    arma_extra_debug_sigprint();
    }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return U.is_alias(X); }
  };



template<typename T1, bool condition>
struct Proxy_xtrans_redirect {};

template<typename T1>
struct Proxy_xtrans_redirect<T1, false> { typedef Proxy_xtrans_default<T1> result; };

template<typename T1>
struct Proxy_xtrans_redirect<T1, true>  { typedef Proxy_xtrans_vector<T1>  result; };



template<typename T1>
class Proxy< Op<T1, op_htrans> >
  : public
    Proxy_xtrans_redirect
      <
      Op<T1, op_htrans>,
      ((is_complex<typename T1::elem_type>::value == false) && ((Op<T1, op_htrans>::is_row) || (Op<T1, op_htrans>::is_col)) )
      >::result
  {
  public:
  
  typedef
  typename
  Proxy_xtrans_redirect
    <
    Op<T1, op_htrans>,
    ((is_complex<typename T1::elem_type>::value == false) && ((Op<T1, op_htrans>::is_row) || (Op<T1, op_htrans>::is_col)) )
    >::result
  Proxy_xtrans;
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  
  static const bool prefer_at_accessor = Proxy_xtrans::prefer_at_accessor;
  static const bool has_subview        = Proxy_xtrans::has_subview;
  static const bool is_fixed           = Proxy_xtrans::is_fixed;
  static const bool fake_mat           = Proxy_xtrans::fake_mat;
  
  // NOTE: the Op class takes care of swapping row and col for op_htrans
  static const bool is_row = Op<T1, op_htrans>::is_row;
  static const bool is_col = Op<T1, op_htrans>::is_col;
  
  using Proxy_xtrans::Q;
  
  inline explicit Proxy(const Op<T1, op_htrans>& A)
    : Proxy_xtrans(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return Proxy_xtrans::is_alias(X); }
  };



template<typename T1>
class Proxy< Op<T1, op_strans> >
  : public
    Proxy_xtrans_redirect
      <
      Op<T1, op_strans>,
      ( (Op<T1, op_strans>::is_row) || (Op<T1, op_strans>::is_col) )
      >::result
  {
  public:
  
  typedef
  typename
  Proxy_xtrans_redirect
    <
    Op<T1, op_strans>,
    ( (Op<T1, op_strans>::is_row) || (Op<T1, op_strans>::is_col) )
    >::result
  Proxy_xtrans;
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  
  static const bool prefer_at_accessor = Proxy_xtrans::prefer_at_accessor;
  static const bool has_subview        = Proxy_xtrans::has_subview;
  static const bool is_fixed           = Proxy_xtrans::is_fixed;
  static const bool fake_mat           = Proxy_xtrans::fake_mat;
  
  // NOTE: the Op class takes care of swapping row and col for op_strans
  static const bool is_row = Op<T1, op_strans>::is_row;
  static const bool is_col = Op<T1, op_strans>::is_col;
  
  using Proxy_xtrans::Q;
  
  inline explicit Proxy(const Op<T1, op_strans>& A)
    : Proxy_xtrans(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return Proxy_xtrans::is_alias(X); }
  };



template<typename T1, typename T2, typename glue_type>
class Proxy< Glue<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = Glue<T1, T2, glue_type>::is_row;
  static const bool is_col = Glue<T1, T2, glue_type>::is_col;
  
  arma_aligned const Mat<elem_type> Q;
  
  inline explicit Proxy(const Glue<T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }

  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>&) const { return false; }
  };



template<typename eT>
class Proxy< subview<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview<eT>                              stored_type;
  typedef const subview<eT>&                       ea_type;
  
  static const bool prefer_at_accessor = true;
  static const bool has_subview        = true;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = false;
  static const bool is_col = false;
  
  arma_aligned const subview<eT>& Q;
  
  inline explicit Proxy(const subview<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.m)) == void_ptr(&X)); }
  };



template<typename eT>
class Proxy< subview_col<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_col<eT>                          stored_type;
  typedef const eT*                                ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = true;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = false;
  static const bool is_col = true;
  
  arma_aligned const subview_col<eT>& Q;
  
  inline explicit Proxy(const subview_col<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];   }
  arma_inline elem_type at         (const uword row, const uword) const { return Q[row]; }
  
  arma_inline ea_type get_ea() const { return Q.colptr(0); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.m)) == void_ptr(&X)); }
  };



template<typename eT>
class Proxy< subview_row<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_row<eT>                          stored_type;
  typedef const subview_row<eT>&                   ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = true;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = true;
  static const bool is_col = false;
  
  arma_aligned const subview_row<eT>& Q;
  
  inline explicit Proxy(const subview_row<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return 1;        }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];   }
  arma_inline elem_type at         (const uword, const uword col) const { return Q[col]; }
  
  arma_inline ea_type get_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.m)) == void_ptr(&X)); }
  };



template<typename eT, typename T1>
class Proxy< subview_elem1<eT,T1> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = false;
  static const bool is_col = true;
  
  arma_aligned const Mat<eT> Q;
  
  inline explicit Proxy(const subview_elem1<eT,T1>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];   }
  arma_inline elem_type at         (const uword row, const uword) const { return Q[row]; }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>&) const { return false; }
  };



template<typename eT, typename T1, typename T2>
class Proxy< subview_elem2<eT,T1,T2> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = false;
  static const bool is_col = false;
  
  arma_aligned const Mat<eT> Q;
  
  inline explicit Proxy(const subview_elem2<eT,T1,T2>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>&) const { return false; }
  };



template<typename eT>
class Proxy< diagview<eT> >
  {
  public:
  
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef diagview<eT>                             stored_type;
  typedef const diagview<eT>&                      ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = true;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = false;
  static const bool is_col = true;
  
  arma_aligned const diagview<eT>& Q;
  
  inline explicit Proxy(const diagview<eT>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];         }
  arma_inline elem_type at         (const uword row, const uword) const { return Q.at(row, 0); }
  
  arma_inline ea_type get_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.m)) == void_ptr(&X)); }
  };




template<typename T1, typename eop_type>
class Proxy< eOp<T1, eop_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eOp<T1, eop_type>                        stored_type;
  typedef const eOp<T1, eop_type>&                 ea_type;
  
  static const bool prefer_at_accessor = eOp<T1, eop_type>::prefer_at_accessor;
  static const bool has_subview        = eOp<T1, eop_type>::has_subview;
  static const bool is_fixed           = eOp<T1, eop_type>::is_fixed;
  static const bool fake_mat           = eOp<T1, eop_type>::fake_mat;
  
  static const bool is_row = eOp<T1, eop_type>::is_row;
  static const bool is_col = eOp<T1, eop_type>::is_col;
  
  arma_aligned const eOp<T1, eop_type>& Q;
  
  inline explicit Proxy(const eOp<T1, eop_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.get_n_rows(); }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.get_n_cols(); }
  arma_inline uword get_n_elem() const { return Q.get_n_elem();              }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return Q.P.is_alias(X); }
  };



template<typename T1, typename T2, typename eglue_type>
class Proxy< eGlue<T1, T2, eglue_type > >
  {
  public:
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eGlue<T1, T2, eglue_type>                stored_type;
  typedef const eGlue<T1, T2, eglue_type>&         ea_type;
  
  static const bool prefer_at_accessor = eGlue<T1, T2, eglue_type>::prefer_at_accessor;
  static const bool has_subview        = eGlue<T1, T2, eglue_type>::has_subview;
  static const bool is_fixed           = eGlue<T1, T2, eglue_type>::is_fixed;
  static const bool fake_mat           = eGlue<T1, T2, eglue_type>::fake_mat;
  
  static const bool is_row = eGlue<T1, T2, eglue_type>::is_row;
  static const bool is_col = eGlue<T1, T2, eglue_type>::is_col;
  
  arma_aligned const eGlue<T1, T2, eglue_type>& Q;
  
  inline explicit Proxy(const eGlue<T1, T2, eglue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.get_n_rows(); }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.get_n_cols(); }
  arma_inline uword get_n_elem() const { return Q.get_n_elem();              }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];           }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row, col); }
  
  arma_inline ea_type get_ea()                          const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (Q.P1.is_alias(X) || Q.P2.is_alias(X)); }
  };



template<typename out_eT, typename T1, typename op_type>
class Proxy< mtOp<out_eT, T1, op_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Mat<out_eT>                  stored_type;
  typedef          const elem_type*             ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = mtOp<out_eT, T1, op_type>::is_row;
  static const bool is_col = mtOp<out_eT, T1, op_type>::is_col;
  
  arma_aligned const Mat<out_eT> Q;
  
  inline explicit Proxy(const mtOp<out_eT, T1, op_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];          }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row,col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>&) const { return false; }
  };



template<typename out_eT, typename T1, typename T2, typename glue_type>
class Proxy< mtGlue<out_eT, T1, T2, glue_type > >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Mat<out_eT>                  stored_type;
  typedef          const elem_type*             ea_type;
  
  static const bool prefer_at_accessor = false;
  static const bool has_subview        = false;
  static const bool is_fixed           = false;
  static const bool fake_mat           = false;
  
  static const bool is_row = mtGlue<out_eT, T1, T2, glue_type>::is_row;
  static const bool is_col = mtGlue<out_eT, T1, T2, glue_type>::is_col;
  
  arma_aligned const Mat<out_eT> Q;
  
  inline explicit Proxy(const mtGlue<out_eT, T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                    const { return Q[i];          }
  arma_inline elem_type at         (const uword row, const uword col) const { return Q.at(row,col); }
  
  arma_inline ea_type get_ea() const { return Q.memptr(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>&) const { return false; }
  };



//! @}
