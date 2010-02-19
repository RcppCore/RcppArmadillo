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
	template <typename T> SEXP wrap ( const arma::field<T>& ) ;
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

template <typename T>
SEXP arma_wrap( const T& object, const ::Rcpp::Dimension& dim){
	::Rcpp::RObject x = ::Rcpp::wrap( object.memptr() , object.memptr() + object.n_elem ) ;
	x.attr( "dim" ) = dim ;
	return x; 
}

} /* namespace RcppArmadillo */
	
/* wrap */

template <typename T> SEXP wrap ( const arma::Mat<T>& data ){
	return RcppArmadillo::arma_wrap( data, Dimension( data.n_rows, data.n_cols ) ) ;
} ;

template <typename T> SEXP wrap( const arma::Col<T>& data ){
	return RcppArmadillo::arma_wrap( data, Dimension( data.n_elem, 1) ) ;
} ;

template <typename T> SEXP wrap( const arma::Row<T>& data ){
	return RcppArmadillo::arma_wrap(data, Dimension( 1, data.n_elem ) ) ;
} ;

#ifdef HAS_CUBE
template <typename T> SEXP wrap( const arma::Cube<T>& cube ){
	return RcppArmadillo::arma_wrap(data, Dimension(  data.n_rows, data.n_cols, data.n_slices ) ) ;
}
#endif

namespace RcppArmadillo {
	
/* Importer class for field<T> */
template <typename T> class FieldImporter {
public:
	typedef T r_import_type ;
	FieldImporter( const arma::field<T>& data_ ) : data(data_){}
	inline int size() const { return data.n_elem ; }
	inline T get(int i) const { return data[i] ; }
	SEXP wrap( int i) const { return ::Rcpp::wrap( data[i] ) ; } 
private:
	const arma::field<T>& data ;
} ;

} // namespace RcppArmadillo

template <typename T> SEXP wrap( const arma::field<T>& data){
	RObject x = wrap( RcppArmadillo::FieldImporter<T>( data ) ) ;
	x.attr("dim" ) = Dimension( data.n_rows, data.n_cols ) ;
	return x ;
}




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

