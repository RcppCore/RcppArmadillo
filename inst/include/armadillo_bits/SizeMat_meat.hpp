// Copyright (C) 2013-2015 Conrad Sanderson
// Copyright (C) 2013-2015 NICTA (www.nicta.com.au)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup SizeMat
//! @{



inline
SizeMat::SizeMat(const uword in_n_rows, const uword in_n_cols)
  : n_rows(in_n_rows)
  , n_cols(in_n_cols)
  {
  arma_extra_debug_sigprint();
  }



inline
uword
SizeMat::operator[](const uword dim) const
  {
  if(dim == 0)  { return n_rows; }
  if(dim == 1)  { return n_cols; }
  
  return uword(0);
  }



inline
uword
SizeMat::operator()(const uword dim) const
  {
  if(dim == 0)  { return n_rows; }
  if(dim == 1)  { return n_cols; }
  
  arma_debug_check(true, "size(): dimension out of bounds");
  
  return uword(0);
  }



inline
bool
SizeMat::operator==(const SizeMat& s) const
  {
  if(n_rows != s.n_rows)  { return false; }
  
  if(n_cols != s.n_cols)  { return false; }
  
  return true;
  }



inline
bool
SizeMat::operator!=(const SizeMat& s) const
  {
  if(n_rows != s.n_rows)  { return true; }
  
  if(n_cols != s.n_cols)  { return true; }
  
  return false;
  }



inline
bool
SizeMat::operator==(const SizeCube& s) const
  {
  if(n_rows   != s.n_rows  )  { return false; }
  
  if(n_cols   != s.n_cols  )  { return false; }
  
  if(uword(1) != s.n_slices)  { return false; }
  
  return true;
  }



inline
bool
SizeMat::operator!=(const SizeCube& s) const
  {
  if(n_rows   != s.n_rows  )  { return true; }
  
  if(n_cols   != s.n_cols  )  { return true; }
  
  if(uword(1) != s.n_slices)  { return true; }
  
  return false;
  }



//! @}
