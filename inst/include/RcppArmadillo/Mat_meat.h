// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Mat_meat.h: Rcpp/Armadillo glue
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

#ifndef RCPPARMADILLO_MAT_MEAT_H
#define RCPPARMADILLO_MAT_MEAT_H

template <typename eT>
template <int RTYPE, bool NA, typename VECTOR>
inline Mat<eT>::Mat( const Rcpp::VectorBase<RTYPE,NA,VECTOR>& X ) 
	: n_rows( 0 )
	, n_cols( 0 )
	, n_elem( 0 )
	, use_aux_mem(false)
	, mem(mem)
{
	
	arma_extra_debug_sigprint_this(this);
	
	// TODO : deal with complex expressions because 
	// std::complex<double> != Rcomplex
	isnt_same_type<eT, typename Rcpp::traits::storage_type<RTYPE>::type >::check();
  	
	set_size( X.size(), 1 ) ;
		
	eT* ptr = memptr() ;
	for( u32 i=0; i<n_elem; ++i){
		ptr[i] = X[i] ;
	}
}

template <typename eT>
template <int RTYPE, bool NA, typename MATRIX>
inline Mat<eT>::Mat( const Rcpp::MatrixBase<RTYPE,NA,MATRIX>& X ) 
	: n_rows( 0 )
	, n_cols( 0 )
	, n_elem( 0 )
	, use_aux_mem(false)
	, mem(mem)
{
	
	arma_extra_debug_sigprint_this(this);
	
	// TODO : deal with complex expressions because 
	// std::complex<double> != Rcomplex
	isnt_same_type<eT, typename Rcpp::traits::storage_type<RTYPE>::type >::check();
  	
	u32 nr = X.nrow(), nc = X.ncol(), i_col, i_row, k ;
	set_size( nr, nc ) ;
		
	eT* ptr = memptr() ;
	for( i_col=0, k=0 ; i_col < nc; ++i_col){
		for( i_row = 0; i_row < nr ; ++i_row, ++k ){
			ptr[k] = X(i_row,i_col) ;
		}
	}
}

#endif
