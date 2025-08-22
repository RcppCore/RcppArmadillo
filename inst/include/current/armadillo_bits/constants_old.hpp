// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (https://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// https://www.apache.org/licenses/LICENSE-2.0
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
struct Math
  {
  [[deprecated("use datum::pi instead")]]      static eT pi()      { return eT(Datum<eT>::pi);      }
  [[deprecated("use datum::e instead")]]       static eT e()       { return eT(Datum<eT>::e);       }
  [[deprecated("use datum::euler instead")]]   static eT euler()   { return eT(Datum<eT>::euler);   }
  [[deprecated("use datum::gratio instead")]]  static eT gratio()  { return eT(Datum<eT>::gratio);  }
  [[deprecated("use datum::sqrt2 instead")]]   static eT sqrt2()   { return eT(Datum<eT>::sqrt2);   }
  [[deprecated("use datum::eps instead")]]     static eT eps()     { return eT(Datum<eT>::eps);     }
  [[deprecated("use datum::log_min instead")]] static eT log_min() { return eT(Datum<eT>::log_min); }
  [[deprecated("use datum::log_max instead")]] static eT log_max() { return eT(Datum<eT>::log_max); }
  [[deprecated("use datum::nan instead")]]     static eT nan()     { return eT(Datum<eT>::nan);     }
  [[deprecated("use datum::inf instead")]]     static eT inf()     { return eT(Datum<eT>::inf);     }
  };



template<typename eT>
struct Phy
  {
  [[deprecated]] static eT m_u()       { return eT(Datum<eT>::m_u);       }
  [[deprecated]] static eT N_A()       { return eT(Datum<eT>::N_A);       }
  [[deprecated]] static eT k()         { return eT(Datum<eT>::k);         }
  [[deprecated]] static eT k_evk()     { return eT(Datum<eT>::k_evk);     }
  [[deprecated]] static eT a_0()       { return eT(Datum<eT>::a_0);       }
  [[deprecated]] static eT mu_B()      { return eT(Datum<eT>::mu_B);      }
  [[deprecated]] static eT Z_0()       { return eT(Datum<eT>::Z_0);       }
  [[deprecated]] static eT G_0()       { return eT(Datum<eT>::G_0);       }
  [[deprecated]] static eT k_e()       { return eT(Datum<eT>::k_e);       }
  [[deprecated]] static eT eps_0()     { return eT(Datum<eT>::eps_0);     }
  [[deprecated]] static eT m_e()       { return eT(Datum<eT>::m_e);       }
  [[deprecated]] static eT eV()        { return eT(Datum<eT>::eV);        }
  [[deprecated]] static eT e()         { return eT(Datum<eT>::ec);        }
  [[deprecated]] static eT F()         { return eT(Datum<eT>::F);         }
  [[deprecated]] static eT alpha()     { return eT(Datum<eT>::alpha);     }
  [[deprecated]] static eT alpha_inv() { return eT(Datum<eT>::alpha_inv); }
  [[deprecated]] static eT K_J()       { return eT(Datum<eT>::K_J);       }
  [[deprecated]] static eT mu_0()      { return eT(Datum<eT>::mu_0);      }
  [[deprecated]] static eT phi_0()     { return eT(Datum<eT>::phi_0);     }
  [[deprecated]] static eT R()         { return eT(Datum<eT>::R);         }
  [[deprecated]] static eT G()         { return eT(Datum<eT>::G);         }
  [[deprecated]] static eT h()         { return eT(Datum<eT>::h);         }
  [[deprecated]] static eT h_bar()     { return eT(Datum<eT>::h_bar);     }
  [[deprecated]] static eT m_p()       { return eT(Datum<eT>::m_p);       }
  [[deprecated]] static eT R_inf()     { return eT(Datum<eT>::R_inf);     }
  [[deprecated]] static eT c_0()       { return eT(Datum<eT>::c_0);       }
  [[deprecated]] static eT sigma()     { return eT(Datum<eT>::sigma);     }
  [[deprecated]] static eT R_k()       { return eT(Datum<eT>::R_k);       }
  [[deprecated]] static eT b()         { return eT(Datum<eT>::b);         }
  };



typedef Math<float>  fmath;
typedef Math<double> math;

typedef Phy<float>   fphy;
typedef Phy<double>  phy;



//! @}
