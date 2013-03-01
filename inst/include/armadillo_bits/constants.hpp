// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup constants
//! @{


namespace priv
  {
  class Datum_helper
    {
    public:
    
    template<typename eT>
    static
    typename arma_real_only<eT>::result
    nan(typename arma_real_only<eT>::result* junk = 0)
      {
      arma_ignore(junk);
      
      if(std::numeric_limits<eT>::has_quiet_NaN == true)
        {
        return std::numeric_limits<eT>::quiet_NaN();
        }
      else
        {
        return eT(0);
        }
      }
    
    
    template<typename eT>
    static
    typename arma_cx_only<eT>::result
    nan(typename arma_cx_only<eT>::result* junk = 0)
      {
      arma_ignore(junk);
      
      typedef typename get_pod_type<eT>::result T;
      
      return eT( Datum_helper::nan<T>(), Datum_helper::nan<T>() );
      }
    
    
    template<typename eT>
    static
    typename arma_integral_only<eT>::result
    nan(typename arma_integral_only<eT>::result* junk = 0)
      {
      arma_ignore(junk);
      
      return eT(0);
      }
    
    
    template<typename eT>
    static
    typename arma_real_only<eT>::result
    inf(typename arma_real_only<eT>::result* junk = 0)
      {
      arma_ignore(junk);
      
      if(std::numeric_limits<eT>::has_infinity == true)
        {
        return std::numeric_limits<eT>::infinity();
        }
      else
        {
        return std::numeric_limits<eT>::max();
        }
      }
    
    
    template<typename eT>
    static
    typename arma_cx_only<eT>::result
    inf(typename arma_cx_only<eT>::result* junk = 0)
      {
      arma_ignore(junk);
      
      typedef typename get_pod_type<eT>::result T;
      
      return eT( Datum_helper::inf<T>(), Datum_helper::inf<T>() );
      }
    

    template<typename eT>
    static
    typename arma_integral_only<eT>::result
    inf(typename arma_integral_only<eT>::result* junk = 0)
      {
      arma_ignore(junk);
      
      return std::numeric_limits<eT>::max();
      }
    
    };
  }



//! various constants.
//! Physical constants taken from NIST and WolframAlpha on 2009-06-23
//! http://physics.nist.gov/cuu/Constants
//! http://www.wolframalpha.com
//! See also http://en.wikipedia.org/wiki/Physical_constant


template<typename eT>
class Datum
  {
  public:
  
  static const eT pi;       //!< ratio of any circle's circumference to its diameter
  static const eT e;        //!< base of the natural logarithm
  static const eT euler;    //!< Euler's constant, aka Euler-Mascheroni constant
  static const eT gratio;   //!< golden ratio
  static const eT sqrt2;    //!< square root of 2
  static const eT eps;      //!< the difference between 1 and the least value greater than 1 that is representable
  static const eT log_min;  //!< log of the minimum representable value
  static const eT log_max;  //!< log of the maximum representable value
  static const eT nan;      //!< "not a number"
  static const eT inf;      //!< infinity 

  // 
  
  static const eT m_u;       //!< atomic mass constant (in kg)
  static const eT N_A;       //!< Avogadro constant
  static const eT k;         //!< Boltzmann constant (in joules per kelvin)
  static const eT k_evk;     //!< Boltzmann constant (in eV/K)
  static const eT a_0;       //!< Bohr radius (in meters)
  static const eT mu_B;      //!< Bohr magneton
  static const eT Z_0;       //!< characteristic impedance of vacuum (in ohms)
  static const eT G_0;       //!< conductance quantum (in siemens)
  static const eT k_e;       //!< Coulomb's constant (in meters per farad)
  static const eT eps_0;     //!< electric constant (in farads per meter)
  static const eT m_e;       //!< electron mass (in kg)
  static const eT eV;        //!< electron volt (in joules)
  static const eT ec;        //!< elementary charge (in coulombs)
  static const eT F;         //!< Faraday constant (in coulombs)
  static const eT alpha;     //!< fine-structure constant
  static const eT alpha_inv; //!< inverse fine-structure constant
  static const eT K_J;       //!< Josephson constant
  static const eT mu_0;      //!< magnetic constant (in henries per meter)
  static const eT phi_0;     //!< magnetic flux quantum (in webers)
  static const eT R;         //!< molar gas constant (in joules per mole kelvin)
  static const eT G;         //!< Newtonian constant of gravitation (in newton square meters per kilogram squared)
  static const eT h;         //!< Planck constant (in joule seconds)
  static const eT h_bar;     //!< Planck constant over 2 pi, aka reduced Planck constant (in joule seconds)
  static const eT m_p;       //!< proton mass (in kg)
  static const eT R_inf;     //!< Rydberg constant (in reciprocal meters)
  static const eT c_0;       //!< speed of light in vacuum (in meters per second)
  static const eT sigma;     //!< Stefan-Boltzmann constant
  static const eT R_k;       //!< von Klitzing constant (in ohms)
  static const eT b;         //!< Wien wavelength displacement law constant
  };


// the long lengths of the constants are for future support of "long double"
// and any smart compiler that does high-precision computation at compile-time
  
template<typename eT> const eT Datum<eT>::pi        = eT(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679);
template<typename eT> const eT Datum<eT>::e         = eT(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274);
template<typename eT> const eT Datum<eT>::euler     = eT(0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495);
template<typename eT> const eT Datum<eT>::gratio    = eT(1.6180339887498948482045868343656381177203091798057628621354486227052604628189024497072072041893911374);
template<typename eT> const eT Datum<eT>::sqrt2     = eT(1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727);
template<typename eT> const eT Datum<eT>::eps       = std::numeric_limits<eT>::epsilon();
template<typename eT> const eT Datum<eT>::log_min   = std::log(std::numeric_limits<eT>::min());
template<typename eT> const eT Datum<eT>::log_max   = std::log(std::numeric_limits<eT>::max());
template<typename eT> const eT Datum<eT>::nan       = priv::Datum_helper::nan<eT>();
template<typename eT> const eT Datum<eT>::inf       = priv::Datum_helper::inf<eT>();

template<typename eT> const eT Datum<eT>::m_u       = eT(1.660538782e-27);
template<typename eT> const eT Datum<eT>::N_A       = eT(6.02214179e23);
template<typename eT> const eT Datum<eT>::k         = eT(1.3806504e-23);
template<typename eT> const eT Datum<eT>::k_evk     = eT(8.617343e-5);
template<typename eT> const eT Datum<eT>::a_0       = eT(0.52917720859e-10);
template<typename eT> const eT Datum<eT>::mu_B      = eT(927.400915e-26);
template<typename eT> const eT Datum<eT>::Z_0       = eT(3.76730313461771e-2);
template<typename eT> const eT Datum<eT>::G_0       = eT(7.7480917004e-5);
template<typename eT> const eT Datum<eT>::k_e       = eT(8.9875517873681764e9);
template<typename eT> const eT Datum<eT>::eps_0     = eT(8.85418781762039e-12);
template<typename eT> const eT Datum<eT>::m_e       = eT(9.10938215e-31);
template<typename eT> const eT Datum<eT>::eV        = eT(1.602176487e-19);
template<typename eT> const eT Datum<eT>::ec        = eT(1.602176487e-19);
template<typename eT> const eT Datum<eT>::F         = eT(96485.3399);
template<typename eT> const eT Datum<eT>::alpha     = eT(7.2973525376e-3);
template<typename eT> const eT Datum<eT>::alpha_inv = eT(137.035999679);
template<typename eT> const eT Datum<eT>::K_J       = eT(483597.891e9);
template<typename eT> const eT Datum<eT>::mu_0      = eT(1.25663706143592e-06);
template<typename eT> const eT Datum<eT>::phi_0     = eT(2.067833667e-15);
template<typename eT> const eT Datum<eT>::R         = eT(8.314472);
template<typename eT> const eT Datum<eT>::G         = eT(6.67428e-11);
template<typename eT> const eT Datum<eT>::h         = eT(6.62606896e-34);
template<typename eT> const eT Datum<eT>::h_bar     = eT(1.054571628e-34);
template<typename eT> const eT Datum<eT>::m_p       = eT(1.672621637e-27);
template<typename eT> const eT Datum<eT>::R_inf     = eT(10973731.568527);
template<typename eT> const eT Datum<eT>::c_0       = eT(299792458.0);
template<typename eT> const eT Datum<eT>::sigma     = eT(5.670400e-8);
template<typename eT> const eT Datum<eT>::R_k       = eT(25812.807557);
template<typename eT> const eT Datum<eT>::b         = eT(2.8977685e-3);



typedef Datum<float>  fdatum;
typedef Datum<double> datum;




namespace priv
  {
  
  template<typename eT>
  static
  arma_inline
  arma_hot
  typename arma_real_only<eT>::result
  most_neg(typename arma_real_only<eT>::result* junk = 0)
    {
    arma_ignore(junk);
    
    if(std::numeric_limits<eT>::has_infinity == true)
      {
      return -(std::numeric_limits<eT>::infinity());
      }
    else
      {
      return -(std::numeric_limits<eT>::max());
      }
    }
  
  
  template<typename eT>
  static
  arma_inline
  arma_hot
  typename arma_integral_only<eT>::result
  most_neg(typename arma_integral_only<eT>::result* junk = 0)
    {
    arma_ignore(junk);
    
    return std::numeric_limits<eT>::min();
    }
  
  
  template<typename eT>
  static
  arma_inline
  arma_hot
  typename arma_real_only<eT>::result
  most_pos(typename arma_real_only<eT>::result* junk = 0)
    {
    arma_ignore(junk);
    
    if(std::numeric_limits<eT>::has_infinity == true)
      {
      return std::numeric_limits<eT>::infinity();
      }
    else
      {
      return std::numeric_limits<eT>::max();
      }
    }
  
  
  template<typename eT>
  static
  arma_inline
  arma_hot
  typename arma_integral_only<eT>::result
  most_pos(typename arma_integral_only<eT>::result* junk = 0)
    {
    arma_ignore(junk);
    
    return std::numeric_limits<eT>::max();
    }

  }



//! @}
