// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppArmadillo.h: Rcpp/Armadillo glue
//
// Copyright (C)  2010 Dirk Eddelbuettel and Romain Francois
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

#ifndef RcppArmadillo__RcppArmadillo__h
#define RcppArmadillo__RcppArmadillo__h

#include <RcppCommon.h>
#include <armadillo>

#define RCPPARMA_FORWARD(TYPE) \
	template <> arma::Mat<TYPE> as< arma::Mat<TYPE> >( SEXP ) ; \
	template <> arma::Col<TYPE> as< arma::Col<TYPE> >( SEXP ) ; \
	template <> arma::Row<TYPE> as< arma::Row<TYPE> >( SEXP ) ; 

/* forward declarations */
namespace Rcpp{
	template <typename T> SEXP wrap ( const arma::Mat<T>& ) ;
	template <typename T> SEXP wrap ( const arma::Row<T>& ) ;
	template <typename T> SEXP wrap ( const arma::Col<T>& ) ;
#ifdef HAS_CUBE
	template <typename T> SEXP wrap ( const arma::Cube<T>& ) ;
#endif
	RCPPARMA_FORWARD(int)
	RCPPARMA_FORWARD(double)
	RCPPARMA_FORWARD(float)
	RCPPARMA_FORWARD(arma::s32)
}

#include <Rcpp.h>

RcppExport SEXP RcppArmadilloExample() ;
RcppExport SEXP RcppArmadilloExample_as_Mat( SEXP );
RcppExport SEXP RcppArmadilloExample_as_Col( SEXP );
RcppExport SEXP RcppArmadilloExample_as_Row( SEXP );

#define RCPPARMA_COPY(in,out,n) for( int i=0; i<n; i++) { out[i] = in[i] ; } 
#define RCPPARMA_COPY_CAST(in,out,n,TYPE) for( int i=0; i<n; i++) { out[i] = static_cast<TYPE>(in[i]) ; } 

namespace Rcpp{
namespace RcppArmadillo{

/* needs cast : float, unsigned int */
template <typename T> SEXP wrap_Mat__dispatch( const arma::Mat<T>& mat, ::Rcpp::traits::true_type ){
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > vec( 
		::Rcpp::Dimension( mat.n_rows, mat.n_cols ) );
	typedef typename ::Rcpp::traits::storage_type< ::Rcpp::traits::r_sexptype_traits<T>::rtype >::type STORAGE_TYPE ;  
	int n = mat.n_elem ;
	RCPPARMA_COPY_CAST(mat,vec,n,STORAGE_TYPE);
	return vec ;
}

/* no cast : double, int */
template <typename T> SEXP wrap_Mat__dispatch( const arma::Mat<T>& mat, ::Rcpp::traits::false_type ){
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > vec( 
		::Rcpp::Dimension( mat.n_rows, mat.n_cols ) );
	int n = mat.n_elem ;
	RCPPARMA_COPY(mat,vec,n) ;
	return vec ;
}


template <typename T> SEXP wrap_Col__dispatch( const arma::Col<T>& column, ::Rcpp::traits::true_type ){
	int n = column.n_rows ;
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > vec( 
		::Rcpp::Dimension( n, 1 ) );
	typedef typename ::Rcpp::traits::storage_type< ::Rcpp::traits::r_sexptype_traits<T>::rtype >::type STORAGE_TYPE ;  
	RCPPARMA_COPY_CAST(column,vec,n,STORAGE_TYPE) ;
	return vec ;
}

template <typename T> SEXP wrap_Col__dispatch( const arma::Col<T>& column, ::Rcpp::traits::false_type ){
	int n = column.n_rows ;
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > vec( 
		::Rcpp::Dimension( n, 1 ) );
	RCPPARMA_COPY(column,vec,n) ;
	return vec ;
}

template <typename T> SEXP wrap_Row__dispatch( const arma::Row<T>& row, ::Rcpp::traits::true_type){
	int n = row.n_cols ;
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > vec( 
		::Rcpp::Dimension( 1, n ) );
	typedef typename ::Rcpp::traits::storage_type< ::Rcpp::traits::r_sexptype_traits<T>::rtype >::type STORAGE_TYPE ;  
	RCPPARMA_COPY_CAST(row,vec,n,STORAGE_TYPE) ;
	return vec ;
}

template <typename T> SEXP wrap_Row__dispatch( const arma::Row<T>& row, ::Rcpp::traits::false_type){
	int n = row.n_cols ;
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > vec( 
		::Rcpp::Dimension( 1, n ) );
	RCPPARMA_COPY(row,vec,n) ;
	return vec ;
}

} /* namespace RcppArmadillo */
	
/* wrap */

template <typename T> SEXP wrap ( const arma::Mat<T>& mat ){
	return RcppArmadillo::wrap_Mat__dispatch( mat, typename ::Rcpp::traits::r_sexptype_needscast<T>() ) ;
} ;

template <typename T> SEXP wrap( const arma::Col<T>& column ){
	return RcppArmadillo::wrap_Col__dispatch( column, typename ::Rcpp::traits::r_sexptype_needscast<T>() ) ;
}

template <typename T> SEXP wrap( const arma::Row<T>& row ){
	return RcppArmadillo::wrap_Row__dispatch( row, typename ::Rcpp::traits::r_sexptype_needscast<T>() ) ;
}

#ifdef HAS_CUBE
template <typename T> SEXP wrap( const arma::Cube<T>& cube ){
	int n = cube.n_elem;
	SimpleVector< traits::r_sexptype_traits<T>::rtype > vec( 
		Dimension( cube.n_rows, cube.n_cols, cube.n_slices ) );
	RCPPARMA_COPY(cube,vec,n);
	return vec ;
}
#endif

} // namespace Rcpp


namespace Rcpp{
	
namespace RcppArmadillo{

/* when a cast is needed */
template <typename T> ::arma::Mat<T> convert_Mat__dispatch( SEXP x, T, ::Rcpp::traits::true_type ){
	::Rcpp::SimpleMatrix< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	::Rcpp::IntegerVector dim = input.attr("dim") ;
	::arma::Mat<T> out( dim[0], dim[1] ) ;
	int n = dim[0] * dim[1] ;
	RCPPARMA_COPY_CAST(input, out, n, T) ;
	return out;
}

/* when no cast is needed */
template <typename T> ::arma::Mat<T> convert_Mat__dispatch( SEXP x, T, ::Rcpp::traits::false_type ){
	::Rcpp::SimpleMatrix< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	::Rcpp::IntegerVector dim = input.attr("dim") ;
	::arma::Mat<T> out( dim[0], dim[1] ) ;
	int n = dim[0] * dim[1] ;
	RCPPARMA_COPY(input, out, n) ;
	return out;
}

/* when a cast is needed */
template <typename T> ::arma::Col<T> convert_Col__dispatch( SEXP x, T, ::Rcpp::traits::true_type ){
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	int n = input.size() ;
	::arma::Col<T> out( n ) ;
	RCPPARMA_COPY_CAST(input, out, n, T) ;
	return out;
}

/* when no cast is needed */
template <typename T> ::arma::Col<T> convert_Col__dispatch( SEXP x, T, ::Rcpp::traits::false_type ){
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	int n = input.size() ;
	::arma::Col<T> out( n ) ;
	RCPPARMA_COPY(input, out, n) ;
	return out;
}

/* when a cast is needed */
template <typename T> ::arma::Row<T> convert_Row__dispatch( SEXP x, T, ::Rcpp::traits::true_type ){
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	int n = input.size() ;
	::arma::Row<T> out( n ) ;
	RCPPARMA_COPY_CAST(input, out, n, T) ;
	return out;
}

/* when no cast is needed */
template <typename T> ::arma::Row<T> convert_Row__dispatch( SEXP x, T, ::Rcpp::traits::false_type ){
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	int n = input.size() ;
	::arma::Row<T> out( n ) ;
	RCPPARMA_COPY(input, out, n) ;
	return out;
}

/* dispatch depending on whether the type of data in the R vector is the same as T */
template <typename T> ::arma::Mat<T> convert_Mat( SEXP x, T t){
	return convert_Mat__dispatch( x, t, typename ::Rcpp::traits::r_sexptype_needscast<T>() ) ;
}
template <typename T> ::arma::Col<T> convert_Col( SEXP x, T t){
	return convert_Col__dispatch( x, t, typename ::Rcpp::traits::r_sexptype_needscast<T>() ) ;
}
template <typename T> ::arma::Mat<T> convert_Row( SEXP x, T t){
	return convert_Row__dispatch( x, t, typename ::Rcpp::traits::r_sexptype_needscast<T>() ) ;
}

} /* namespace RcppArmadillo */

/* as */

#define GENERATE_CONVERTERS(TYPE)  \
	template<> arma::Mat<TYPE> as< arma::Mat<TYPE> >(SEXP x){ return RcppArmadillo::convert_Mat<TYPE>(x, TYPE()) ; } ; \
	template<> arma::Col<TYPE> as< arma::Col<TYPE> >(SEXP x){ return RcppArmadillo::convert_Col<TYPE>(x, TYPE()) ; } ; \
	template<> arma::Row<TYPE> as< arma::Row<TYPE> >(SEXP x){ return RcppArmadillo::convert_Row<TYPE>(x, TYPE()) ; } ;

GENERATE_CONVERTERS(int)
GENERATE_CONVERTERS(arma::u32)
GENERATE_CONVERTERS(double)
GENERATE_CONVERTERS(float)

#undef GENERATE_CONVERTER
#undef RCPPARMA_COPY
#undef RCPPARMA_COPY_CAST
#undef RCPPARMA_FORWARD

}

#endif

