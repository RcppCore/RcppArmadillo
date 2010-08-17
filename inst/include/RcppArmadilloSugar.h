// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppArmadilloSugar.h: Rcpp/Armadillo glue
//
// Copyright (C)  2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
//
// This file is part of RcppArmadillo.
//
// RcppArmadillo is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppArmadillo is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RcppArmadillo__RcppArmadilloSugar__h
#define RcppArmadillo__RcppArmadilloSugar__h

namespace Rcpp{
namespace RcppArmadillo{

	template <int RTYPE, bool NA, typename T>
	class r_forward {
		public:
			
		typedef typename Rcpp::traits::storage_type<RTYPE>::type T1 ;
		typedef typename Rcpp::VectorBase<RTYPE,NA,T> SUGAR ;
		typedef typename arma::Mat<T1>::iterator mat_iterator ;
		
		inline static void apply(
			arma::Mat<T1>& out, 
			const arma::Op< Rcpp::VectorBase<RTYPE,NA,T> , r_forward<RTYPE,NA,T> >& in
		){
			const SUGAR& m = in.m ;
			int n = m.size() ;
			// deal with dimensions
			out.set_size( n, 1 ) ;
			mat_iterator first = out.begin(), last = out.end();
			// perhaps we should just use std::copy
			for( int i=0; first != last ; ++i){
				*first++ = m[i];
			}
		}
	} ;

	template <bool NA, typename T>
	class Complex_Imposter {
		public:
			typedef typename Rcpp::VectorBase<CPLXSXP,NA,T> SUGAR_EXP ;
			typedef std::complex<double> elem_type ;
		
			Complex_Imposter( const SUGAR_EXP& vec_) : vec(vec_){}
			inline int size() const { return vec.size() ; }
			inline std::complex<double> operator[]( int i ) const {
				Rcomplex x = vec[i] ;
				return std::complex<double>( x.r, x.i ) ;
			}
			
		private:
			const SUGAR_EXP& vec ;
	} ;

	
	
	template <bool NA, typename T>
	class r_complex_forward {
		public:
			typedef std::complex<double> T1 ;
			typedef Complex_Imposter<NA,T> SUGAR ;
			typedef arma::Mat<T1>::iterator mat_iterator ;
		
		inline static void apply(
			arma::Mat<T1>& out, 
			const arma::Op< Complex_Imposter<NA,T> , r_complex_forward<NA,T> >& in
		){
			const SUGAR& m = in.m ;
			int n = m.size() ;
			// deal with dimensions
			out.set_size( n, 1 ) ;
			mat_iterator first = out.begin(), last = out.end();
			for( int i=0; first != last ; ++i){
				*first++ = m[i] ;
			}
		}
			
	} ;
	
}


template <int RTYPE, bool NA, typename T>
inline arma::Op< 
	VectorBase<RTYPE,NA,T> , 
	RcppArmadillo::r_forward<RTYPE,NA,T>
>
forward( const VectorBase<RTYPE,NA,T>& x ) {
	return arma::Op< VectorBase<RTYPE,NA,T> , RcppArmadillo::r_forward<RTYPE,NA,T> >( x ) ;
}


template <bool NA, typename T>
inline arma::Op< 
	RcppArmadillo::Complex_Imposter<NA,T> , 
	RcppArmadillo::r_complex_forward<NA,T>
>
forward( const VectorBase<CPLXSXP,NA,T>& x ) {
	return arma::Op< 
		RcppArmadillo::Complex_Imposter<NA,T> ,
		RcppArmadillo::r_complex_forward<NA,T> >( x ) ;
}



} // Rcpp

#endif

