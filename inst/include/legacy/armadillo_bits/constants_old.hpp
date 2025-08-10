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


//! \addtogroup constants_old
//! @{


// DO NOT USE IN NEW CODE !!!
// the Math and Phy classes are kept for compatibility with old code;
// for new code, use the Datum class instead
// eg. instead of math::pi(), use datum::pi

template<typename eT>
class Math
  {
  public:
  
  arma_frown("use datum::pi instead")      static eT pi()      { return eT(Datum<eT>::pi);      }
  arma_frown("use datum::e instead")       static eT e()       { return eT(Datum<eT>::e);       }
  arma_frown("use datum::euler instead")   static eT euler()   { return eT(Datum<eT>::euler);   }
  arma_frown("use datum::gratio instead")  static eT gratio()  { return eT(Datum<eT>::gratio);  }
  arma_frown("use datum::sqrt2 instead")   static eT sqrt2()   { return eT(Datum<eT>::sqrt2);   }
  arma_frown("use datum::eps instead")     static eT eps()     { return eT(Datum<eT>::eps);     }
  arma_frown("use datum::log_min instead") static eT log_min() { return eT(Datum<eT>::log_min); }
  arma_frown("use datum::log_max instead") static eT log_max() { return eT(Datum<eT>::log_max); }
  arma_frown("use datum::nan instead")     static eT nan()     { return eT(Datum<eT>::nan);     }
  arma_frown("use datum::inf instead")     static eT inf()     { return eT(Datum<eT>::inf);     }
  };



template<typename eT>
class Phy
  {
  public:
  
  arma_deprecated static eT m_u()       { return eT(Datum<eT>::m_u);       }
  arma_deprecated static eT N_A()       { return eT(Datum<eT>::N_A);       }
  arma_deprecated static eT k()         { return eT(Datum<eT>::k);         }
  arma_deprecated static eT k_evk()     { return eT(Datum<eT>::k_evk);     }
  arma_deprecated static eT a_0()       { return eT(Datum<eT>::a_0);       }
  arma_deprecated static eT mu_B()      { return eT(Datum<eT>::mu_B);      }
  arma_deprecated static eT Z_0()       { return eT(Datum<eT>::Z_0);       }
  arma_deprecated static eT G_0()       { return eT(Datum<eT>::G_0);       }
  arma_deprecated static eT k_e()       { return eT(Datum<eT>::k_e);       }
  arma_deprecated static eT eps_0()     { return eT(Datum<eT>::eps_0);     }
  arma_deprecated static eT m_e()       { return eT(Datum<eT>::m_e);       }
  arma_deprecated static eT eV()        { return eT(Datum<eT>::eV);        }
  arma_deprecated static eT e()         { return eT(Datum<eT>::ec);        }
  arma_deprecated static eT F()         { return eT(Datum<eT>::F);         }
  arma_deprecated static eT alpha()     { return eT(Datum<eT>::alpha);     }
  arma_deprecated static eT alpha_inv() { return eT(Datum<eT>::alpha_inv); }
  arma_deprecated static eT K_J()       { return eT(Datum<eT>::K_J);       }
  arma_deprecated static eT mu_0()      { return eT(Datum<eT>::mu_0);      }
  arma_deprecated static eT phi_0()     { return eT(Datum<eT>::phi_0);     }
  arma_deprecated static eT R()         { return eT(Datum<eT>::R);         }
  arma_deprecated static eT G()         { return eT(Datum<eT>::G);         }
  arma_deprecated static eT h()         { return eT(Datum<eT>::h);         }
  arma_deprecated static eT h_bar()     { return eT(Datum<eT>::h_bar);     }
  arma_deprecated static eT m_p()       { return eT(Datum<eT>::m_p);       }
  arma_deprecated static eT R_inf()     { return eT(Datum<eT>::R_inf);     }
  arma_deprecated static eT c_0()       { return eT(Datum<eT>::c_0);       }
  arma_deprecated static eT sigma()     { return eT(Datum<eT>::sigma);     }
  arma_deprecated static eT R_k()       { return eT(Datum<eT>::R_k);       }
  arma_deprecated static eT b()         { return eT(Datum<eT>::b);         }
  };



typedef Math<float>  fmath;
typedef Math<double> math;

typedef Phy<float>   fphy;
typedef Phy<double>  phy;



//! @}
