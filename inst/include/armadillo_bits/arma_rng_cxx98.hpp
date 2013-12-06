// Copyright (C) 2013 Conrad Sanderson
// Copyright (C) 2013 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup arma_rng_cxx98
//! @{



class arma_rng_cxx98
  {
  public:
  
  typedef unsigned int seed_type;
  
  inline static void set_seed(const seed_type val);
  
  arma_inline static int    randi_val();
  arma_inline static double randu_val();
       inline static double randn_val();
  
  template<typename eT>
  inline static void randn_dual_val(eT& out1, eT& out2);
  
  template<typename eT>
  inline static void randi_fill(eT* mem, const uword N, const int a, const int b);
  
  inline static int randi_max_val();
  };



inline
void
arma_rng_cxx98::set_seed(const arma_rng_cxx98::seed_type val)
  {
  std::srand(val);
  }



arma_inline
int
arma_rng_cxx98::randi_val()
  {
  return std::rand();
  }



arma_inline
double
arma_rng_cxx98::randu_val()
  {
  return double( double(std::rand()) * ( double(1) / double(RAND_MAX) ) );
  }



inline
double
arma_rng_cxx98::randn_val()
  {
  // polar form of the Box-Muller transformation:
  // http://en.wikipedia.org/wiki/Box-Muller_transformation
  // http://en.wikipedia.org/wiki/Marsaglia_polar_method
  
  double tmp1;
  double tmp2;
  double w;
  
  do
    {
    tmp1 = double(2) * double(std::rand()) * (double(1) / double(RAND_MAX)) - double(1);
    tmp2 = double(2) * double(std::rand()) * (double(1) / double(RAND_MAX)) - double(1);
    
    w = tmp1*tmp1 + tmp2*tmp2;
    }
  while ( w >= double(1) );
  
  return double( tmp1 * std::sqrt( (double(-2) * std::log(w)) / w) );
  }



template<typename eT>
inline
void
arma_rng_cxx98::randn_dual_val(eT& out1, eT& out2)
  {
  // make sure we are internally using at least floats
  typedef typename promote_type<eT,float>::result eTp;
  
  eTp tmp1;
  eTp tmp2;
  eTp w;
  
  do
    {
    tmp1 = eTp(2) * eTp(std::rand()) * (eTp(1) / eTp(RAND_MAX)) - eTp(1);
    tmp2 = eTp(2) * eTp(std::rand()) * (eTp(1) / eTp(RAND_MAX)) - eTp(1);
    
    w = tmp1*tmp1 + tmp2*tmp2;
    }
  while ( w >= eTp(1) );
  
  const eTp k = std::sqrt( (eTp(-2) * std::log(w)) / w);
  
  out1 = eT(tmp1*k);
  out2 = eT(tmp2*k);
  }



template<typename eT>
inline
void
arma_rng_cxx98::randi_fill(eT* mem, const uword N, const int a, const int b)
  {
  if( (a == 0) && (b == RAND_MAX) )
    {
    for(uword i=0; i<N; ++i)
      {
      mem[i] = std::rand();
      }
    }
  else
    {
    const uword length = b - a + 1;
    
    const double scale = double(length) / double(RAND_MAX);
    
    for(uword i=0; i<N; ++i)
      {
      mem[i] = int( double(std::rand()) * scale ) + a;
      }
    }
  }



inline
int
arma_rng_cxx98::randi_max_val()
  {
  return RAND_MAX;
  }



//! @}
