// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup restrictors
//! @{



template<typename T> struct arma_scalar_only { };

template<> struct arma_scalar_only<char>   { typedef char   result; };
template<> struct arma_scalar_only<short>  { typedef short  result; };
template<> struct arma_scalar_only<int>    { typedef int    result; };
template<> struct arma_scalar_only<long>   { typedef long   result; };
template<> struct arma_scalar_only<float>  { typedef float  result; };
template<> struct arma_scalar_only<double> { typedef double result; };

template<> struct arma_scalar_only<unsigned char>  { typedef unsigned char  result; };
template<> struct arma_scalar_only<unsigned short> { typedef unsigned short result; };
template<> struct arma_scalar_only<unsigned int>   { typedef unsigned int   result; };
template<> struct arma_scalar_only<unsigned long>  { typedef unsigned long  result; };

template<typename T> struct arma_scalar_only< std::complex<T> > { typedef std::complex<T> result; };



template<typename T> struct arma_integral_only { };

template<> struct arma_integral_only<char>   { typedef char   result; };
template<> struct arma_integral_only<short>  { typedef short  result; };
template<> struct arma_integral_only<int>    { typedef int    result; };
template<> struct arma_integral_only<long>   { typedef long   result; };

template<> struct arma_integral_only<unsigned char>  { typedef unsigned char  result; };
template<> struct arma_integral_only<unsigned short> { typedef unsigned short result; };
template<> struct arma_integral_only<unsigned int>   { typedef unsigned int   result; };
template<> struct arma_integral_only<unsigned long>  { typedef unsigned long  result; };



template<typename T> struct arma_unsigned_integral_only { };

template<> struct arma_unsigned_integral_only<unsigned char>  { typedef unsigned char  result; };
template<> struct arma_unsigned_integral_only<unsigned short> { typedef unsigned short result; };
template<> struct arma_unsigned_integral_only<unsigned int>   { typedef unsigned int   result; };
template<> struct arma_unsigned_integral_only<unsigned long>  { typedef unsigned long  result; };



template<typename T> struct arma_signed_integral_only { };

template<> struct arma_signed_integral_only<char>   { typedef char   result; };
template<> struct arma_signed_integral_only<short>  { typedef short  result; };
template<> struct arma_signed_integral_only<int>    { typedef int    result; };
template<> struct arma_signed_integral_only<long>   { typedef long   result; };



template<typename T> struct arma_float_only { };

template<> struct arma_float_only<float>  { typedef float  result; };
template<> struct arma_float_only<double> { typedef double result; };



//! @}
