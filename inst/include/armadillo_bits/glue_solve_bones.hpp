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



//! \addtogroup glue_solve
//! @{



class glue_solve_gen
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = false;
    static constexpr bool is_col  = T2::is_col;
    static constexpr bool is_xvec = false;
    };
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_gen>& X);
  
  template<typename eT, typename T1, typename T2> inline static bool apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags);
  };



class glue_solve_tri_default
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = false;
    static constexpr bool is_col  = T2::is_col;
    static constexpr bool is_xvec = false;
    };
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_tri_default>& X);
  
  template<typename eT, typename T1, typename T2> inline static bool apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags);
  };



class glue_solve_tri
  {
  public:
  
  template<typename T1, typename T2>
  struct traits
    {
    static constexpr bool is_row  = false;
    static constexpr bool is_col  = T2::is_col;
    static constexpr bool is_xvec = false;
    };
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_tri>& X);
  
  template<typename eT, typename T1, typename T2> inline static bool apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags);
  };



namespace solve_opts
  {
  struct opts
    {
    const uword flags;
    
    inline explicit opts(const uword in_flags);
    
    inline const opts operator+(const opts& rhs) const;
    };
  
  inline
  opts::opts(const uword in_flags)
    : flags(in_flags)
    {}
  
  inline
  const opts
  opts::operator+(const opts& rhs) const
    {
    const opts result( flags | rhs.flags );
    
    return result;
    }
  
  // The values below (eg. 1u << 1) are for internal Armadillo use only.
  // The values can change without notice.
  
  static constexpr uword flag_none         = uword(0       );
  static constexpr uword flag_fast         = uword(1u <<  0);
  static constexpr uword flag_equilibrate  = uword(1u <<  1);
  static constexpr uword flag_no_approx    = uword(1u <<  2);
  static constexpr uword flag_triu         = uword(1u <<  3);
  static constexpr uword flag_tril         = uword(1u <<  4);
  static constexpr uword flag_no_band      = uword(1u <<  5);
  static constexpr uword flag_no_sympd     = uword(1u <<  6);
  static constexpr uword flag_allow_ugly   = uword(1u <<  7);
  static constexpr uword flag_likely_sympd = uword(1u <<  8);
  static constexpr uword flag_refine       = uword(1u <<  9);
  static constexpr uword flag_no_trimat    = uword(1u << 10);
  
  struct opts_none         : public opts { inline opts_none()         : opts(flag_none        ) {} };
  struct opts_fast         : public opts { inline opts_fast()         : opts(flag_fast        ) {} };
  struct opts_equilibrate  : public opts { inline opts_equilibrate()  : opts(flag_equilibrate ) {} };
  struct opts_no_approx    : public opts { inline opts_no_approx()    : opts(flag_no_approx   ) {} };
  struct opts_triu         : public opts { inline opts_triu()         : opts(flag_triu        ) {} };
  struct opts_tril         : public opts { inline opts_tril()         : opts(flag_tril        ) {} };
  struct opts_no_band      : public opts { inline opts_no_band()      : opts(flag_no_band     ) {} };
  struct opts_no_sympd     : public opts { inline opts_no_sympd()     : opts(flag_no_sympd    ) {} };
  struct opts_allow_ugly   : public opts { inline opts_allow_ugly()   : opts(flag_allow_ugly  ) {} };
  struct opts_likely_sympd : public opts { inline opts_likely_sympd() : opts(flag_likely_sympd) {} };
  struct opts_refine       : public opts { inline opts_refine()       : opts(flag_refine      ) {} };
  struct opts_no_trimat    : public opts { inline opts_no_trimat()    : opts(flag_no_trimat   ) {} };
  
  static const opts_none         none;
  static const opts_fast         fast;
  static const opts_equilibrate  equilibrate;
  static const opts_no_approx    no_approx;
  static const opts_triu         triu;
  static const opts_tril         tril;
  static const opts_no_band      no_band;
  static const opts_no_sympd     no_sympd;
  static const opts_allow_ugly   allow_ugly;
  static const opts_likely_sympd likely_sympd;
  static const opts_refine       refine;
  static const opts_no_trimat    no_trimat;
  }



//! @}
