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


//! \addtogroup arma_rng
//! @{


#if defined(ARMA_RNG_ALT)
  #undef ARMA_USE_EXTERN_RNG
#endif


// NOTE: mt19937_64_instance_warmup is used as a workaround
// NOTE: for thread_local issue on macOS 11 and/or AppleClang 12.0
// NOTE: see https://gitlab.com/conradsnicta/armadillo-code/-/issues/173
// NOTE: if this workaround causes problems, please report it and 
// NOTE: disable the workaround by uncommenting the code block below:

// #if defined(__APPLE__) || defined(__apple_build_version__)
//   #if !defined(ARMA_DONT_DISABLE_EXTERN_RNG)
//     #undef ARMA_USE_EXTERN_RNG
//   #endif
// #endif


// NOTE: workaround for another thread_local issue on macOS
// NOTE: where GCC (not Clang) may not have support for thread_local

#if (defined(__APPLE__) && defined(__GNUG__) && !defined(__clang__))
  #if !defined(ARMA_DONT_DISABLE_EXTERN_RNG)
    #undef ARMA_USE_EXTERN_RNG
  #endif
#endif



#if defined(ARMA_USE_EXTERN_RNG)
  extern thread_local std::mt19937_64 mt19937_64_instance;
  
  #if defined(__APPLE__) || defined(__apple_build_version__)
    namespace
      {
      struct mt19937_64_instance_warmup
        {
        inline mt19937_64_instance_warmup()
          {
          typename std::mt19937_64::result_type junk = mt19937_64_instance();
          arma_ignore(junk);
          }
        };
      
      static mt19937_64_instance_warmup mt19937_64_instance_warmup_run;
      }
  #endif
#endif



class arma_rng
  {
  public:
  
  #if   defined(ARMA_RNG_ALT)
    typedef arma_rng_alt::seed_type      seed_type;
  #elif defined(ARMA_USE_EXTERN_RNG)
    typedef std::mt19937_64::result_type seed_type;
  #else
    typedef arma_rng_cxx03::seed_type    seed_type;
  #endif
  
  #if   defined(ARMA_RNG_ALT)
    static constexpr int rng_method = 2;
  #elif defined(ARMA_USE_EXTERN_RNG)
    static constexpr int rng_method = 1;
  #else
    static constexpr int rng_method = 0;
  #endif
  
  inline static void set_seed(const seed_type val);
  inline static void set_seed_random();
  
  template<typename eT> struct randi;
  template<typename eT> struct randu;
  template<typename eT> struct randn;
  template<typename eT> struct randg;
  };



inline
void
arma_rng::set_seed(const arma_rng::seed_type val)
  {
  #if   defined(ARMA_RNG_ALT)
    {
    arma_rng_alt::set_seed(val);
    }
  #elif defined(ARMA_USE_EXTERN_RNG)
    {
    mt19937_64_instance.seed(val);
    }
  #else
    {
    arma_rng_cxx03::set_seed(val);
    }
  #endif
  }



arma_cold
inline
void
arma_rng::set_seed_random()
  {
  seed_type seed1 = seed_type(0);
  seed_type seed2 = seed_type(0);
  seed_type seed3 = seed_type(0);
  seed_type seed4 = seed_type(0);
  
  bool have_seed = false;
  
  try
    {
    std::random_device rd;
    
    if(rd.entropy() > double(0))  { seed1 = static_cast<seed_type>( rd() ); }
    
    if(seed1 != seed_type(0))  { have_seed = true; }
    }
  catch(...) {}
  
  
  if(have_seed == false)
    {
    try
      {
      union
        {
        seed_type     a;
        unsigned char b[sizeof(seed_type)];
        } tmp;
      
      tmp.a = seed_type(0);
      
      std::ifstream f("/dev/urandom", std::ifstream::binary);
      
      if(f.good())  { f.read((char*)(&(tmp.b[0])), sizeof(seed_type)); }
      
      if(f.good())
        {
        seed2 = tmp.a;
        
        if(seed2 != seed_type(0))  { have_seed = true; }
        }
      }
    catch(...) {}
    }
  
  
  if(have_seed == false)
    {
    // get better-than-nothing seeds in case reading /dev/urandom failed 
    
    const std::chrono::system_clock::time_point tp_now = std::chrono::system_clock::now();
    
    auto since_epoch_usec = std::chrono::duration_cast<std::chrono::microseconds>(tp_now.time_since_epoch()).count();
    
    seed3 = static_cast<seed_type>( since_epoch_usec & 0xFFFF );
    
    union
      {
      uword*        a;
      unsigned char b[sizeof(uword*)];
      } tmp;
    
    tmp.a = (uword*)malloc(sizeof(uword));
    
    if(tmp.a != nullptr)
      {
      for(size_t i=0; i<sizeof(uword*); ++i)  { seed4 += seed_type(tmp.b[i]); }
      
      free(tmp.a);
      }
    }
  
  arma_rng::set_seed( seed1 + seed2 + seed3 + seed4 );
  }



//



template<typename eT>
struct arma_rng::randi
  {
  inline
  operator eT ()
    {
    #if   defined(ARMA_RNG_ALT)
      {
      return eT( arma_rng_alt::randi_val() );
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      constexpr double scale = double(std::numeric_limits<int>::max()) / double(std::mt19937_64::max());
      
      return eT( double(mt19937_64_instance()) * scale );
      }
    #else
      {
      return eT( arma_rng_cxx03::randi_val() );
      }
    #endif
    }
  
  
  inline
  static
  int
  max_val()
    {
    #if   defined(ARMA_RNG_ALT)
      {
      return arma_rng_alt::randi_max_val();
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      return std::numeric_limits<int>::max();
      }
    #else
      {
      return arma_rng_cxx03::randi_max_val();
      }
    #endif
    }
  
  
  inline
  static
  void
  fill(eT* mem, const uword N, const int a, const int b)
    {
    #if   defined(ARMA_RNG_ALT)
      {
      arma_rng_alt::randi_fill(mem, N, a, b);
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::uniform_int_distribution<int> local_i_distr(a, b);
      
      for(uword i=0; i<N; ++i)  { mem[i] = eT(local_i_distr(mt19937_64_instance)); }
      }
    #else
      {
      if(N == uword(1))  { arma_rng_cxx03::randi_fill(mem, uword(1), a, b); return; }
      
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                    local_engine;
      std::uniform_int_distribution<int> local_i_distr(a, b);
      
      local_engine.seed( local_seed_type(std::rand()) );
      
      for(uword i=0; i<N; ++i)  { mem[i] = eT(local_i_distr(local_engine)); }
      }
    #endif
    }
  };



//



template<typename eT>
struct arma_rng::randu
  {
  inline
  operator eT ()
    {
    #if   defined(ARMA_RNG_ALT)
      {
      return eT( arma_rng_alt::randu_val() );
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      constexpr double scale = double(1.0) / double(std::mt19937_64::max());
      
      return eT( double(mt19937_64_instance()) * scale );
      }
    #else
      {
      return eT( arma_rng_cxx03::randu_val() );
      }
    #endif
    }
  
  
  inline
  static
  void
  fill(eT* mem, const uword N)
    {
    #if defined(ARMA_RNG_ALT)
      {
      for(uword i=0; i < N; ++i)  { mem[i] = eT( arma_rng_alt::randu_val() ); }
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::uniform_real_distribution<double> local_u_distr;
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( local_u_distr(mt19937_64_instance) ); }
      }
    #else
      {
      if(N == uword(1))  { mem[0] = eT( arma_rng_cxx03::randu_val() ); return; }
      
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                        local_engine;
      std::uniform_real_distribution<double> local_u_distr;
      
      local_engine.seed( local_seed_type(std::rand()) );
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( local_u_distr(local_engine) ); }
      }
    #endif
    }
  
  
  inline
  static
  void
  fill(eT* mem, const uword N, const double a, const double b)
    {
    #if defined(ARMA_RNG_ALT)
      {
      const double r = b - a;
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( arma_rng_alt::randu_val() * r + a ); }
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::uniform_real_distribution<double> local_u_distr(a,b);
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( local_u_distr(mt19937_64_instance) ); }
      }
    #else
      {
      if(N == uword(1))  { mem[0] = eT( arma_rng_cxx03::randu_val() * (b - a) + a ); return; }
      
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                        local_engine;
      std::uniform_real_distribution<double> local_u_distr(a,b);
      
      local_engine.seed( local_seed_type(std::rand()) );
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( local_u_distr(local_engine) ); }
      }
    #endif
    }
  };



template<typename T>
struct arma_rng::randu< std::complex<T> >
  {
  arma_inline
  operator std::complex<T> ()
    {
    #if defined(ARMA_RNG_ALT)
      {
      const T a = T( arma_rng_alt::randu_val() );
      const T b = T( arma_rng_alt::randu_val() );
      
      return std::complex<T>(a, b);
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::uniform_real_distribution<double> local_u_distr;
      
      const T a = T( local_u_distr(mt19937_64_instance) );
      const T b = T( local_u_distr(mt19937_64_instance) );
      
      return std::complex<T>(a, b);
      }
    #else
      {
      const T a = T( arma_rng_cxx03::randu_val() );
      const T b = T( arma_rng_cxx03::randu_val() );
      
      return std::complex<T>(a, b);
      }
    #endif
    }
  
  
  inline
  static
  void
  fill(std::complex<T>* mem, const uword N)
    {
    #if defined(ARMA_RNG_ALT)
      {
      for(uword i=0; i < N; ++i)
        {
        const T a = T( arma_rng_alt::randu_val() );
        const T b = T( arma_rng_alt::randu_val() );
        
        mem[i] = std::complex<T>(a, b);
        }
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::uniform_real_distribution<double> local_u_distr;
      
      for(uword i=0; i < N; ++i)
        {
        const T a = T( local_u_distr(mt19937_64_instance) );
        const T b = T( local_u_distr(mt19937_64_instance) );
        
        mem[i] = std::complex<T>(a, b);
        }
      }
    #else
      {
      if(N == uword(1))
        {
        const T a = T( arma_rng_cxx03::randu_val() );
        const T b = T( arma_rng_cxx03::randu_val() );
        
        mem[0] = std::complex<T>(a, b);
        
        return;
        }
      
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                        local_engine;
      std::uniform_real_distribution<double> local_u_distr;
      
      local_engine.seed( local_seed_type(std::rand()) );
      
      for(uword i=0; i < N; ++i)
        {
        const T a = T( local_u_distr(local_engine) );
        const T b = T( local_u_distr(local_engine) );
        
        mem[i] = std::complex<T>(a, b);
        }
      }
    #endif
    }
  
  
  inline
  static
  void
  fill(std::complex<T>* mem, const uword N, const double a, const double b)
    {
    #if defined(ARMA_RNG_ALT)
      {
      const double r = b - a;
      
      for(uword i=0; i < N; ++i)
        {
        const T tmp1 = T( arma_rng_alt::randu_val() * r + a );
        const T tmp2 = T( arma_rng_alt::randu_val() * r + a );
        
        mem[i] = std::complex<T>(tmp1, tmp2);
        }
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::uniform_real_distribution<double> local_u_distr(a,b);
      
      for(uword i=0; i < N; ++i)
        {
        const T tmp1 = T( local_u_distr(mt19937_64_instance) );
        const T tmp2 = T( local_u_distr(mt19937_64_instance) );
        
        mem[i] = std::complex<T>(tmp1, tmp2);
        }
      }
    #else
      {
      if(N == uword(1))
        {
        const double r = b - a;
        
        const T tmp1 = T( arma_rng_cxx03::randu_val() * r + a);
        const T tmp2 = T( arma_rng_cxx03::randu_val() * r + a);
        
        mem[0] = std::complex<T>(tmp1, tmp2);
        
        return;
        }
      
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                        local_engine;
      std::uniform_real_distribution<double> local_u_distr(a,b);
      
      local_engine.seed( local_seed_type(std::rand()) );
      
      for(uword i=0; i < N; ++i)
        {
        const T tmp1 = T( local_u_distr(local_engine) );
        const T tmp2 = T( local_u_distr(local_engine) );
        
        mem[i] = std::complex<T>(tmp1, tmp2);
        }
      }
    #endif
    }
  };



//



template<typename eT>
struct arma_rng::randn
  {
  inline
  operator eT () const
    {
    #if   defined(ARMA_RNG_ALT)
      {
      return eT( arma_rng_alt::randn_val() );
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::normal_distribution<double> local_n_distr;
      
      return eT( local_n_distr(mt19937_64_instance) );
      }
    #else
      {
      return eT( arma_rng_cxx03::randn_val() );
      }
    #endif
    }
  
  
  inline
  static
  void
  dual_val(eT& out1, eT& out2)
    {
    #if   defined(ARMA_RNG_ALT)
      {
      arma_rng_alt::randn_dual_val(out1, out2);
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::normal_distribution<double> local_n_distr;
      
      out1 = eT( local_n_distr(mt19937_64_instance) );
      out2 = eT( local_n_distr(mt19937_64_instance) );
      }
    #else
      {
      arma_rng_cxx03::randn_dual_val(out1, out2);
      }
    #endif
    }
  
  
  inline
  static
  void
  fill(eT* mem, const uword N)
    {
    #if defined(ARMA_RNG_ALT)
      {
      // NOTE: old method to avoid regressions in user code that assumes specific sequence
      
      uword i, j;
      
      for(i=0, j=1; j < N; i+=2, j+=2)  { arma_rng_alt::randn_dual_val( mem[i], mem[j] ); }
      
      if(i < N)  { mem[i] = eT( arma_rng_alt::randn_val() ); }
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::normal_distribution<double> local_n_distr;
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( local_n_distr(mt19937_64_instance) ); }
      }
    #else
      {
      if(N == uword(1))  { mem[0] = eT( arma_rng_cxx03::randn_val() ); return; }
      
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                  local_engine;
      std::normal_distribution<double> local_n_distr;
      
      local_engine.seed( local_seed_type(std::rand()) );
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( local_n_distr(local_engine) ); }
      }
    #endif
    }
  
  
  inline
  static
  void
  fill(eT* mem, const uword N, const double mu, const double sd)
    {
    #if defined(ARMA_RNG_ALT)
      {
      // NOTE: old method to avoid regressions in user code that assumes specific sequence
      
      uword i, j;
      
      for(i=0, j=1; j < N; i+=2, j+=2)
        {
        eT val_i = eT(0);
        eT val_j = eT(0);
        
        arma_rng_alt::randn_dual_val( val_i, val_j );
        
        mem[i] = (val_i * sd) + mu;
        mem[j] = (val_j * sd) + mu;
        }
      
      if(i < N)
        {
        const eT val_i = eT( arma_rng_alt::randn_val() );
         
        mem[i] = (val_i * sd) + mu;
        }
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::normal_distribution<double> local_n_distr(mu, sd);
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( local_n_distr(mt19937_64_instance) ); }
      }
    #else
      {
      if(N == uword(1))
        {
        const eT val = eT( arma_rng_cxx03::randn_val() );
        
        mem[0] = (val * sd) + mu;
        
        return;
        }
      
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                  local_engine;
      std::normal_distribution<double> local_n_distr(mu, sd);
      
      local_engine.seed( local_seed_type(std::rand()) );
      
      for(uword i=0; i < N; ++i)  { mem[i] = eT( local_n_distr(local_engine) ); }
      }
    #endif
    }
  };



template<typename T>
struct arma_rng::randn< std::complex<T> >
  {
  inline
  operator std::complex<T> () const
    {
    #if defined(_MSC_VER)
      // attempt at workaround for MSVC bug
      // does MS even test their so-called compilers before release?
      T a;
      T b;
    #else
      T a(0);
      T b(0);
    #endif
    
    arma_rng::randn<T>::dual_val(a, b);
    
    return std::complex<T>(a, b);
    }
  
  
  inline
  static
  void
  dual_val(std::complex<T>& out1, std::complex<T>& out2)
    {
    #if defined(_MSC_VER)
      T a;
      T b;
    #else
      T a(0);
      T b(0);
    #endif
    
    arma_rng::randn<T>::dual_val(a,b);
    out1 = std::complex<T>(a,b);
    
    arma_rng::randn<T>::dual_val(a,b);
    out2 = std::complex<T>(a,b);
    }
  
  
  inline
  static
  void
  fill(std::complex<T>* mem, const uword N)
    {
    #if defined(ARMA_RNG_ALT)
      {
      for(uword i=0; i < N; ++i)  { mem[i] = std::complex<T>( arma_rng::randn< std::complex<T> >() ); }
      }
    #elif defined(ARMA_USE_EXTERN_RNG)
      {
      std::normal_distribution<double> local_n_distr;
      
      for(uword i=0; i < N; ++i)
        {
        const T a = T( local_n_distr(mt19937_64_instance) );
        const T b = T( local_n_distr(mt19937_64_instance) );
        
        mem[i] = std::complex<T>(a,b);
        }
      }
    #else
      {
      if(N == uword(1))
        {
        T a = T(0);
        T b = T(0);
        
        arma_rng_cxx03::randn_dual_val(a,b);
        
        mem[0] = std::complex<T>(a,b);
        
        return;
        }
      
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                  local_engine;
      std::normal_distribution<double> local_n_distr;
      
      local_engine.seed( local_seed_type(std::rand()) );
      
      for(uword i=0; i < N; ++i)
        {
        const T a = T( local_n_distr(local_engine) );
        const T b = T( local_n_distr(local_engine) );
        
        mem[i] = std::complex<T>(a,b);
        }
      }
    #endif
    }
  
  
  inline
  static
  void
  fill(std::complex<T>* mem, const uword N, const double mu, const double sd)
    {
    arma_rng::randn< std::complex<T> >::fill(mem, N);
    
    if( (mu == double(0)) && (sd == double(1)) )  { return; }
    
    for(uword i=0; i<N; ++i)
      {
      const std::complex<T>& val = mem[i];
      
      mem[i] = std::complex<T>( ((val.real() * sd) + mu), ((val.imag() * sd) + mu) );
      }
    }
  };



//



template<typename eT>
struct arma_rng::randg
  {
  inline
  static
  void
  fill(eT* mem, const uword N, const double a, const double b)
    {
    #if defined(ARMA_USE_EXTERN_RNG)
      {
      std::gamma_distribution<double> local_g_distr(a,b);
      
      for(uword i=0; i<N; ++i)  { mem[i] = eT(local_g_distr(mt19937_64_instance)); }
      }
    #else
      {
      typedef typename std::mt19937_64::result_type local_seed_type;
      
      std::mt19937_64                 local_engine;
      std::gamma_distribution<double> local_g_distr(a,b);
      
      local_engine.seed( local_seed_type(arma_rng::randi<local_seed_type>()) );
      
      for(uword i=0; i<N; ++i)  { mem[i] = eT(local_g_distr(local_engine)); }
      }
    #endif
    }
  };



//! @}
