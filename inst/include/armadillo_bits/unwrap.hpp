// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup unwrap
//! @{



template<typename T1>
struct unwrap_default
  {
  typedef typename T1::elem_type eT;
  
  inline unwrap_default(const T1& A)   // TODO: change this to Base ?
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Mat<eT> M;
  };



template<typename T1>
struct unwrap_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit unwrap_Mat_fixed(const T1& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct unwrap_redirect {};

template<typename T1>
struct unwrap_redirect<T1, false> { typedef unwrap_default<T1>   result; };

template<typename T1>
struct unwrap_redirect<T1, true>  { typedef unwrap_Mat_fixed<T1> result; };


template<typename T1>
class unwrap : public unwrap_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  public:
  
  inline unwrap(const T1& A)
    : unwrap_redirect< T1, is_Mat_fixed<T1>::value >::result(A)
    {
    }
  };



template<typename eT>
class unwrap< Mat<eT> >
  {
  public:
  
  inline unwrap(const Mat<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Mat<eT>& M;
  };



template<typename eT>
class unwrap< Row<eT> >
  {
  public:
  
  inline unwrap(const Row<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Row<eT>& M;
  };



template<typename eT>
class unwrap< Col<eT> >
  {
  public:
  
  inline unwrap(const Col<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Col<eT>& M;
  };



template<typename out_eT, typename T1, typename T2, typename glue_type>
class unwrap< mtGlue<out_eT, T1, T2, glue_type> >
  {
  public:
  
  inline unwrap(const mtGlue<out_eT, T1, T2, glue_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Mat<out_eT> M;
  };


template<typename out_eT, typename T1, typename op_type>
class unwrap< mtOp<out_eT, T1, op_type> >
  {
  public:
  
  inline unwrap(const mtOp<out_eT, T1, op_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Mat<out_eT> M;
  };



//
//
//


template<typename T1>
class unwrap_check
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_check(const T1& A, const Mat<eT>&)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Mat<eT> M;
  };



template<typename eT>
class unwrap_check< Mat<eT> >
  {
  public:

  inline
  unwrap_check(const Mat<eT>& A, const Mat<eT>& B)
    : M_local( (&A == &B) ? new Mat<eT>(A) : 0 )
    , M      ( (&A == &B) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class unwrap_check< Row<eT> >
  {
  public:
  
  inline
  unwrap_check(const Row<eT>& A, const Mat<eT>& B)
    : M_local( (&A == reinterpret_cast<const Row<eT>*>(&B)) ? new Row<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Row<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  
  // the order below is important
  const Row<eT>* M_local;
  const Row<eT>& M;
  };



template<typename eT>
class unwrap_check< Col<eT> >
  {
  public:

  inline
  unwrap_check(const Col<eT>& A, const Mat<eT>& B)
    : M_local( (&A == reinterpret_cast<const Col<eT>*>(&B)) ? new Col<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Col<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  
  // the order below is important
  const Col<eT>* M_local;
  const Col<eT>& M;
  
  };



//
//
//



template<typename T1>
class unwrap_check_mixed
  {
  public:
  
  typedef typename T1::elem_type eT1;
  
  template<typename eT2>
  inline
  unwrap_check_mixed(const T1& A, const Mat<eT2>&)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Mat<eT1> M;
  };



template<typename eT1>
class unwrap_check_mixed< Mat<eT1> >
  {
  public:
  
  template<typename eT2>
  inline
  unwrap_check_mixed(const Mat<eT1>& A, const Mat<eT2>& B)
    : M_local( (void_ptr(&A) == void_ptr(&B)) ? new Mat<eT1>(A) : 0 )
    , M      ( (void_ptr(&A) == void_ptr(&B)) ? (*M_local)      : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check_mixed()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  
  // the order below is important
  const Mat<eT1>* M_local;
  const Mat<eT1>& M;
  };



template<typename eT1>
class unwrap_check_mixed< Row<eT1> >
  {
  public:
  
  template<typename eT2>
  inline
  unwrap_check_mixed(const Row<eT1>& A, const Mat<eT2>& B)
    : M_local( (void_ptr(&A) == void_ptr(&B)) ? new Row<eT1>(A) : 0 )
    , M      ( (void_ptr(&A) == void_ptr(&B)) ? (*M_local)      : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check_mixed()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  
  // the order below is important
  const Row<eT1>* M_local;
  const Row<eT1>& M;
  };



template<typename eT1>
class unwrap_check_mixed< Col<eT1> >
  {
  public:
  
  template<typename eT2>
  inline
  unwrap_check_mixed(const Col<eT1>& A, const Mat<eT2>& B)
    : M_local( (void_ptr(&A) == void_ptr(&B)) ? new Col<eT1>(A) : 0 )
    , M      ( (void_ptr(&A) == void_ptr(&B)) ? (*M_local)      : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check_mixed()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  
  // the order below is important
  const Col<eT1>* M_local;
  const Col<eT1>& M;
  };



//



template<typename T1>
struct partial_unwrap_default
  {
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_default(const T1& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };


template<typename T1>
struct partial_unwrap_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit partial_unwrap_Mat_fixed(const T1& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct partial_unwrap_redirect {};

template<typename T1>
struct partial_unwrap_redirect<T1, false> { typedef partial_unwrap_default<T1>   result; };

template<typename T1>
struct partial_unwrap_redirect<T1, true>  { typedef partial_unwrap_Mat_fixed<T1> result; };

template<typename T1>
class partial_unwrap : public partial_unwrap_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  public:
  
  inline partial_unwrap(const T1& A)
    : partial_unwrap_redirect< T1, is_Mat_fixed<T1>::value >::result(A)
    {
    }
  };



template<typename eT>
class partial_unwrap< Mat<eT> >
  {
  public:
  
  inline
  partial_unwrap(const Mat<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< Row<eT> >
  {
  public:
  
  inline
  partial_unwrap(const Row<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< Col<eT> >
  {
  public:
  
  inline
  partial_unwrap(const Col<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< subview_col<eT> >
  {
  public:
  
  inline
  partial_unwrap(const subview_col<eT>& A)
    : M( const_cast<eT*>( A.colptr(0) ), A.n_rows, 1, false, false )
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_htrans_default
  {
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_htrans_default(const Op<T1, op_htrans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };


template<typename T1>
struct partial_unwrap_htrans_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit partial_unwrap_htrans_Mat_fixed(const Op<T1, op_htrans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct partial_unwrap_htrans_redirect {};

template<typename T1>
struct partial_unwrap_htrans_redirect<T1, false> { typedef partial_unwrap_htrans_default<T1>   result; };

template<typename T1>
struct partial_unwrap_htrans_redirect<T1, true>  { typedef partial_unwrap_htrans_Mat_fixed<T1> result; };

template<typename T1>
class partial_unwrap< Op<T1, op_htrans> > : public partial_unwrap_htrans_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  public:
  
  inline partial_unwrap(const Op<T1, op_htrans>& A)
    : partial_unwrap_htrans_redirect< T1, is_Mat_fixed<T1>::value >::result(A)
    {
    }
  };



template<typename eT>
class partial_unwrap< Op< Mat<eT>, op_htrans> >
  {
  public:
  
  inline
  partial_unwrap(const Op< Mat<eT>, op_htrans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< Op< Row<eT>, op_htrans> >
  {
  public:
  
  inline
  partial_unwrap(const Op< Row<eT>, op_htrans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< Op< Col<eT>, op_htrans> >
  {
  public:
  
  inline
  partial_unwrap(const Op< Col<eT>, op_htrans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< Op< subview_col<eT>, op_htrans> >
  {
  public:
  
  inline
  partial_unwrap(const Op< subview_col<eT>, op_htrans>& A)
    : M( const_cast<eT*>( A.m.colptr(0) ), A.m.n_rows, 1, false, false )
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_htrans2_default
  {
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_htrans2_default(const Op<T1, op_htrans2>& A)
    : val(A.aux)
    , M  (A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };


template<typename T1>
struct partial_unwrap_htrans2_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit partial_unwrap_htrans2_Mat_fixed(const Op<T1, op_htrans2>& A)
    : val(A.aux)
    , M  (A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct partial_unwrap_htrans2_redirect {};

template<typename T1>
struct partial_unwrap_htrans2_redirect<T1, false> { typedef partial_unwrap_htrans2_default<T1>   result; };

template<typename T1>
struct partial_unwrap_htrans2_redirect<T1, true>  { typedef partial_unwrap_htrans2_Mat_fixed<T1> result; };

template<typename T1>
class partial_unwrap< Op<T1, op_htrans2> > : public partial_unwrap_htrans2_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  public:
  
  inline partial_unwrap(const Op<T1, op_htrans2>& A)
    : partial_unwrap_htrans2_redirect< T1, is_Mat_fixed<T1>::value >::result(A)
    {
    }
  };



template<typename eT>
class partial_unwrap< Op< Mat<eT>, op_htrans2> >
  {
  public:
  
  inline
  partial_unwrap(const Op< Mat<eT>, op_htrans2>& A)
    : val(A.aux)
    , M  (A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< Op< Row<eT>, op_htrans2> >
  {
  public:
  
  inline
  partial_unwrap(const Op< Row<eT>, op_htrans2>& A)
    : val(A.aux)
    , M  (A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< Op< Col<eT>, op_htrans2> >
  {
  public:
  
  inline
  partial_unwrap(const Op< Col<eT>, op_htrans2>& A)
    : val(A.aux)
    , M  (A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< Op< subview_col<eT>, op_htrans2> >
  {
  public:
  
  inline
  partial_unwrap(const Op< subview_col<eT>, op_htrans2>& A)
    : val( A.aux )
    , M  ( const_cast<eT*>( A.m.colptr(0) ), A.m.n_rows, 1, false, false )
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_scalar_times_default
  {
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_scalar_times_default(const eOp<T1, eop_scalar_times>& A)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_scalar_times_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit partial_unwrap_scalar_times_Mat_fixed(const eOp<T1, eop_scalar_times>& A)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct partial_unwrap_scalar_times_redirect {};

template<typename T1>
struct partial_unwrap_scalar_times_redirect<T1, false> { typedef partial_unwrap_scalar_times_default<T1>   result; };

template<typename T1>
struct partial_unwrap_scalar_times_redirect<T1, true>  { typedef partial_unwrap_scalar_times_Mat_fixed<T1> result; };


template<typename T1>
class partial_unwrap< eOp<T1, eop_scalar_times> > : public partial_unwrap_scalar_times_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  typedef typename T1::elem_type eT;
  
  public:
  inline partial_unwrap(const eOp<T1, eop_scalar_times>& A)
    : partial_unwrap_scalar_times_redirect< T1, is_Mat_fixed<T1>::value >::result(A)
    {
    }
  };



template<typename eT>
class partial_unwrap< eOp<Mat<eT>, eop_scalar_times> >
  {
  public:
  
  inline
  partial_unwrap(const eOp<Mat<eT>,eop_scalar_times>& A)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< eOp<Row<eT>, eop_scalar_times> >
  {
  public:
  
  inline
  partial_unwrap(const eOp<Row<eT>,eop_scalar_times>& A)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap< eOp<Col<eT>, eop_scalar_times> >
  {
  public:
  
  inline
  partial_unwrap(const eOp<Col<eT>,eop_scalar_times>& A)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    }
  
  inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



//



template<typename T1>
struct partial_unwrap_check_default
  {
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_check_default(const T1& A, const Mat<eT>&)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };


template<typename T1>
struct partial_unwrap_check_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit partial_unwrap_check_Mat_fixed(const T1& A, const Mat<eT>& B)
    : M_local( (&A == &B) ? new Mat<eT>(A) : 0 )
    , M      ( (&A == &B) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check_Mat_fixed()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct partial_unwrap_check_redirect {};

template<typename T1>
struct partial_unwrap_check_redirect<T1, false> { typedef partial_unwrap_check_default<T1>   result; };

template<typename T1>
struct partial_unwrap_check_redirect<T1, true>  { typedef partial_unwrap_check_Mat_fixed<T1> result; };

template<typename T1>
class partial_unwrap_check : public partial_unwrap_check_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  typedef typename T1::elem_type eT;
  
  public:
  inline partial_unwrap_check(const T1& A, const Mat<eT>& B)
    : partial_unwrap_check_redirect< T1, is_Mat_fixed<T1>::value >::result(A, B)
    {
    }
  };



template<typename eT>
class partial_unwrap_check< Mat<eT> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Mat<eT>& A, const Mat<eT>& B)
    : M_local ( (&A == &B) ? new Mat<eT>(A) : 0 )
    , M       ( (&A == &B) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< Row<eT> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Row<eT>& A, const Mat<eT>& B)
    : M_local ( (&A == &B) ? new Mat<eT>(A) : 0 )
    , M       ( (&A == &B) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< Col<eT> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Col<eT>& A, const Mat<eT>& B)
    : M_local ( (&A == &B) ? new Mat<eT>(A) : 0 )
    , M       ( (&A == &B) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< subview_col<eT> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const subview_col<eT>& A, const Mat<eT>& B)
    : M( const_cast<eT*>( A.colptr(0) ), A.n_rows, 1, (&(A.m) == &B), false )
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_check_htrans_default
  {
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_check_htrans_default(const Op<T1, op_htrans>& A, const Mat<eT>&)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };


template<typename T1>
struct partial_unwrap_check_htrans_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit partial_unwrap_check_htrans_Mat_fixed(const Op<T1, op_htrans>& A, const Mat<eT>& B)
    : M_local( (&(A.m) == &B) ? new Mat<eT>(A.m) : 0   )
    , M      ( (&(A.m) == &B) ? (*M_local)       : A.m )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check_htrans_Mat_fixed()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct partial_unwrap_check_htrans_redirect {};

template<typename T1>
struct partial_unwrap_check_htrans_redirect<T1, false> { typedef partial_unwrap_check_htrans_default<T1>   result; };

template<typename T1>
struct partial_unwrap_check_htrans_redirect<T1, true>  { typedef partial_unwrap_check_htrans_Mat_fixed<T1> result; };


template<typename T1>
class partial_unwrap_check< Op<T1, op_htrans> > : public partial_unwrap_check_htrans_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  typedef typename T1::elem_type eT;
  
  public:
  inline partial_unwrap_check(const Op<T1, op_htrans>& A, const Mat<eT>& B)
    : partial_unwrap_check_htrans_redirect< T1, is_Mat_fixed<T1>::value >::result(A, B)
    {
    }
  };



template<typename eT>
class partial_unwrap_check< Op< Mat<eT>, op_htrans> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Op< Mat<eT>, op_htrans>& A, const Mat<eT>& B)
    : M_local ( (&A.m == &B) ? new Mat<eT>(A.m) : 0   )
    , M       ( (&A.m == &B) ? (*M_local)       : A.m )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< Op< Row<eT>, op_htrans> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Op< Row<eT>, op_htrans>& A, const Mat<eT>& B)
    : M_local ( (&A.m == &B) ? new Mat<eT>(A.m) : 0   )
    , M       ( (&A.m == &B) ? (*M_local)       : A.m )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< Op< Col<eT>, op_htrans> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Op< Col<eT>, op_htrans>& A, const Mat<eT>& B)
    : M_local ( (&A.m == &B) ? new Mat<eT>(A.m) : 0   )
    , M       ( (&A.m == &B) ? (*M_local)       : A.m )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< Op< subview_col<eT>, op_htrans> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Op< subview_col<eT>, op_htrans>& A, const Mat<eT>& B)
    : M( const_cast<eT*>( A.m.colptr(0) ), A.m.n_rows, 1, (&(A.m.m) == &B), false )
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return eT(1); }
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_check_htrans2_default
  {
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_check_htrans2_default(const Op<T1, op_htrans2>& A, const Mat<eT>&)
    : val(A.aux)
    , M  (A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_check_htrans2_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit partial_unwrap_check_htrans2_Mat_fixed(const Op<T1, op_htrans2>& A, const Mat<eT>& B)
    : val    (A.aux)
    , M_local( (&(A.m) == &B) ? new Mat<eT>(A.m) : 0   )
    , M      ( (&(A.m) == &B) ? (*M_local)       : A.m )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check_htrans2_Mat_fixed()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct partial_unwrap_check_htrans2_redirect {};

template<typename T1>
struct partial_unwrap_check_htrans2_redirect<T1, false> { typedef partial_unwrap_check_htrans2_default<T1>   result; };

template<typename T1>
struct partial_unwrap_check_htrans2_redirect<T1, true>  { typedef partial_unwrap_check_htrans2_Mat_fixed<T1> result; };


template<typename T1>
class partial_unwrap_check< Op<T1, op_htrans2> > : public partial_unwrap_check_htrans2_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  typedef typename T1::elem_type eT;
  
  public:
  inline partial_unwrap_check(const Op<T1, op_htrans2>& A, const Mat<eT>& B)
    : partial_unwrap_check_htrans2_redirect< T1, is_Mat_fixed<T1>::value >::result(A, B)
    {
    }
  };



template<typename eT>
class partial_unwrap_check< Op< Mat<eT>, op_htrans2> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Op< Mat<eT>, op_htrans2>& A, const Mat<eT>& B)
    : val     (A.aux)
    , M_local ( (&A.m == &B) ? new Mat<eT>(A.m) : 0   )
    , M       ( (&A.m == &B) ? (*M_local)       : A.m )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  // the order below is important
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< Op< Row<eT>, op_htrans2> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Op< Row<eT>, op_htrans2>& A, const Mat<eT>& B)
    : val     (A.aux)
    , M_local ( (&A.m == &B) ? new Mat<eT>(A.m) : 0   )
    , M       ( (&A.m == &B) ? (*M_local)       : A.m )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  // the order below is important
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< Op< Col<eT>, op_htrans2> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Op< Col<eT>, op_htrans2>& A, const Mat<eT>& B)
    : val     (A.aux)
    , M_local ( (&A.m == &B) ? new Mat<eT>(A.m) : 0   )
    , M       ( (&A.m == &B) ? (*M_local)       : A.m )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  // the order below is important
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< Op< subview_col<eT>, op_htrans2> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const Op< subview_col<eT>, op_htrans2>& A, const Mat<eT>& B)
    : val( A.aux )
    , M  ( const_cast<eT*>( A.m.colptr(0) ), A.m.n_rows, 1, (&(A.m.m) == &B), false )
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_check_scalar_times_default
  {
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_check_scalar_times_default(const eOp<T1, eop_scalar_times>& A, const Mat<eT>&)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };



template<typename T1>
struct partial_unwrap_check_scalar_times_Mat_fixed
  {
  typedef typename T1::elem_type eT;
  
  inline explicit partial_unwrap_check_scalar_times_Mat_fixed(const eOp<T1, eop_scalar_times>& A, const Mat<eT>& B)
    : val    ( A.aux )
    , M_local( (&(A.P.Q) == &B) ? new Mat<eT>(A.P.Q) : 0     )
    , M      ( (&(A.P.Q) == &B) ? (*M_local)         : A.P.Q )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check_scalar_times_Mat_fixed()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename T1, bool condition>
struct partial_unwrap_check_scalar_times_redirect {};

template<typename T1>
struct partial_unwrap_check_scalar_times_redirect<T1, false> { typedef partial_unwrap_check_scalar_times_default<T1>   result; };

template<typename T1>
struct partial_unwrap_check_scalar_times_redirect<T1, true>  { typedef partial_unwrap_check_scalar_times_Mat_fixed<T1> result; };


template<typename T1>
class partial_unwrap_check< eOp<T1, eop_scalar_times> > : public partial_unwrap_check_scalar_times_redirect<T1, is_Mat_fixed<T1>::value >::result
  {
  typedef typename T1::elem_type eT;
  
  public:
  inline partial_unwrap_check(const eOp<T1, eop_scalar_times>& A, const Mat<eT>& B)
    : partial_unwrap_check_scalar_times_redirect< T1, is_Mat_fixed<T1>::value >::result(A, B)
    {
    }
  };



template<typename eT>
class partial_unwrap_check< eOp<Mat<eT>, eop_scalar_times> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const eOp<Mat<eT>,eop_scalar_times>& A, const Mat<eT>& B)
    : val    (A.aux)
    , M_local( (&(A.P.Q) == &B) ? new Mat<eT>(A.P.Q) : 0     )
    , M      ( (&(A.P.Q) == &B) ? *M_local           : A.P.Q )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< eOp<Row<eT>, eop_scalar_times> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const eOp<Row<eT>,eop_scalar_times>& A, const Mat<eT>& B)
    : val(A.aux)
    , M_local( (&(A.P.Q) == &B) ? new Mat<eT>(A.P.Q) : 0     )
    , M      ( (&(A.P.Q) == &B) ? *M_local           : A.P.Q )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< eOp<Col<eT>, eop_scalar_times> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const eOp<Col<eT>,eop_scalar_times>& A, const Mat<eT>& B)
    : val    ( A.aux )
    , M_local( (&(A.P.Q) == &B) ? new Mat<eT>(A.P.Q) : 0     )
    , M      ( (&(A.P.Q) == &B) ? *M_local           : A.P.Q )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local) { delete M_local; }
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename eT>
class partial_unwrap_check< eOp<subview_col<eT>, eop_scalar_times> >
  {
  public:
  
  arma_hot inline
  partial_unwrap_check(const eOp<subview_col<eT>,eop_scalar_times>& A, const Mat<eT>& B)
    : val( A.aux )
    , M  ( const_cast<eT*>( A.P.Q.colptr(0) ), A.P.Q.n_rows, 1, (&(A.P.Q.m) == &B), false )
    {
    arma_extra_debug_sigprint();
    }
  
  arma_hot arma_inline eT get_val() const { return val; }
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };



//! @}
