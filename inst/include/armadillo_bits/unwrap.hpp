// Copyright (C) 2008-2011 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2011 Conrad Sanderson
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
class unwrap
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline unwrap(const T1& A)   // TODO: change this to Base ?
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  const Mat<eT> M;
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
  unwrap_check(const T1& A, const Mat<eT>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    arma_ignore(B);
    }
  
  inline
  ~unwrap_check()
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
    
    if(M_local)
      {
      delete M_local;
      }
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
    
    if(M_local)
      {
      delete M_local;
      }
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
    
    if(M_local)
      {
      delete M_local;
      }
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
  unwrap_check_mixed(const T1& A, const Mat<eT2>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    arma_ignore(B);
    }
  
  inline
  ~unwrap_check_mixed()
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
    
    if(M_local)
      {
      delete M_local;
      }
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
    
    if(M_local)
      {
      delete M_local;
      }
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  // the order below is important
  const Col<eT1>* M_local;
  const Col<eT1>& M;
  };



//



template<typename T1>
class partial_unwrap
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap(const T1& A)  // TODO: change this to Base ?
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT> M;
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
  
  
  inline
  ~partial_unwrap()
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
  
  
  inline
  ~partial_unwrap()
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
  
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename T1>
class partial_unwrap< Op<T1, op_htrans> >
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline
  partial_unwrap(const Op<T1,op_htrans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT> M;
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
  
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
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
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
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
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT>& M;
  };



template<typename T1>
class partial_unwrap< Op<T1, op_htrans2> >
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline
  partial_unwrap(const Op<T1,op_htrans2>& A)
    : val(A.aux)
    , M  (A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return val; }
  
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
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
  
  inline
  ~partial_unwrap()
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
  
  inline
  ~partial_unwrap()
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
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return val; }
  
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



template<typename T1>
class partial_unwrap< eOp<T1, eop_scalar_times> >
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline
  partial_unwrap(const eOp<T1,eop_scalar_times>& A)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return val; }
  
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
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
  
  inline
  ~partial_unwrap()
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
  
  inline
  ~partial_unwrap()
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
  
  inline
  ~partial_unwrap()
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
class partial_unwrap_check
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline partial_unwrap_check(const T1& A, const Mat<eT>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    arma_ignore(B);
    }
  
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };



template<typename eT>
class partial_unwrap_check< Mat<eT> >
  {
  public:
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
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
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
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
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
  static const bool do_trans = false;
  static const bool do_times = false;
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename T1>
class partial_unwrap_check< Op<T1, op_htrans> >
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline
  partial_unwrap_check(const Op<T1,op_htrans>& A, const Mat<eT>& B)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    arma_ignore(B);
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  const Mat<eT> M;
  };



template<typename eT>
class partial_unwrap_check< Op< Mat<eT>, op_htrans> >
  {
  public:
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
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
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
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
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return eT(1); }
  
  
  static const bool do_trans = true;
  static const bool do_times = false;
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename T1>
class partial_unwrap_check< Op<T1, op_htrans2> >
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline
  partial_unwrap_check(const Op<T1,op_htrans2>& A, const Mat<eT>& B)
    : val(A.aux)
    , M  (A.m)
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return val; }
  
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };



template<typename eT>
class partial_unwrap_check< Op< Mat<eT>, op_htrans2> >
  {
  public:
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return val; }
  
  
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
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return val; }
  
  
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
  
  inline
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
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  inline eT get_val() const { return val; }
  
  
  static const bool do_trans = true;
  static const bool do_times = true;
  
  // the order below is important
  const eT       val;
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  };



template<typename T1>
class partial_unwrap_check< eOp<T1, eop_scalar_times> >
  {
  public:
  
  typedef typename T1::elem_type eT;
  
  inline
  partial_unwrap_check(const eOp<T1,eop_scalar_times>& A, const Mat<eT>& B)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    arma_ignore(B);
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return val; }
  
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT      val;
  const Mat<eT> M;
  };



template<typename eT>
class partial_unwrap_check< eOp<Mat<eT>, eop_scalar_times> >
  {
  public:
  
  inline
  partial_unwrap_check(const eOp<Mat<eT>,eop_scalar_times>& A, const Mat<eT>& B)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    arma_ignore(B);
    }
  
  inline
  ~partial_unwrap_check()
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
class partial_unwrap_check< eOp<Row<eT>, eop_scalar_times> >
  {
  public:
  
  inline
  partial_unwrap_check(const eOp<Row<eT>,eop_scalar_times>& A, const Mat<eT>& B)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    arma_ignore(B);
    }
  
  inline
  ~partial_unwrap_check()
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
class partial_unwrap_check< eOp<Col<eT>, eop_scalar_times> >
  {
  public:
  
  inline
  partial_unwrap_check(const eOp<Col<eT>,eop_scalar_times>& A, const Mat<eT>& B)
    : val(A.aux)
    , M  (A.P.Q)
    {
    arma_extra_debug_sigprint();
    arma_ignore(B);
    }
  
  inline
  ~partial_unwrap_check()
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline eT get_val() const { return val; }
  
  
  static const bool do_trans = false;
  static const bool do_times = true;
  
  const eT       val;
  const Mat<eT>& M;
  };



//! @}
