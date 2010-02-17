#ifndef RcppArmadillo__RcppArmadillo__h
#define RcppArmadillo__RcppArmadillo__h

#include <RcppCommon.h>
#include <armadillo>

/* forward declarations */
namespace Rcpp{
	template <typename T> SEXP wrap ( const arma::Mat<T>& ) ;
	template <typename T> SEXP wrap ( const arma::Row<T>& ) ;
	template <typename T> SEXP wrap ( const arma::Col<T>& ) ;
	
	template <> arma::Mat<int> as< arma::Mat<int> >( SEXP ) ;
	template <> arma::Mat<double> as< arma::Mat<double> >( SEXP ) ;
}

#include <Rcpp.h>

RcppExport SEXP RcppArmadilloExample() ;
RcppExport SEXP RcppArmadilloExample_as( SEXP );

namespace Rcpp{

/* wrap */

template <typename T> SEXP wrap ( const arma::Mat<T>& mat ){
	SimpleVector< traits::r_sexptype_traits<T>::rtype > vec( 
		Dimension( mat.n_rows, mat.n_cols ) );
	int n = mat.n_elem ;
	for( int i=0; i<n; i++){
		vec[i] = mat[i] ;
	}
	return vec ; 
} ;

template <typename T> SEXP wrap( const arma::Col<T>& column ){
	int n = column.n_rows ;
	SimpleVector< traits::r_sexptype_traits<T>::rtype > vec( 
		Dimension( n, 1 ) );
	for( int i=0; i<n; i++){
		vec[i] = column[i] ;
	}
	return vec ;
}

template <typename T> SEXP wrap( const arma::Row<T>& row ){
	int n = row.n_cols ;
	SimpleVector< traits::r_sexptype_traits<T>::rtype > vec( 
		Dimension( 1, n ) );
	for( int i=0; i<n; i++){
		vec[i] = row[i] ;
	}
	return vec ;
}

} // namespace Rcpp


namespace Rcpp{
	
namespace RcppArmadillo{

/* when a cast is needed */
template <typename T> ::arma::Mat<T> convert__dispatch( SEXP x, T, ::Rcpp::traits::true_type ){
	if( !Rf_isMatrix(x) ){
		throw std::range_error( "not a matrix" ) ;
	}
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	::Rcpp::IntegerVector dim = input.attr("dim") ;
	::arma::Mat<T> out( dim[0], dim[1] ) ;
	int n = dim[0] * dim[1] ;
	for( int i=0; i<n; i++){
		out[i] = static_cast<T>( input[i] ) ;
	}
	return out;
}

/* when no cast is needed */
template <typename T> ::arma::Mat<T> convert__dispatch( SEXP x, T, ::Rcpp::traits::false_type ){
	if( !Rf_isMatrix(x) ){
		throw std::range_error( "not a matrix" ) ;
	}
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	::Rcpp::IntegerVector dim = input.attr("dim") ;
	::arma::Mat<T> out( dim[0], dim[1] ) ;
	int n = dim[0] * dim[1] ;
	for( int i=0; i<n; i++){
		out[i] = input[i] ;
	}
	return out;
}

/* dispatch depending on whether the type of data in the R vector is the same as T */
template <typename T> ::arma::Mat<T> convert( SEXP x, T t){
	return convert__dispatch( x, t, typename ::Rcpp::traits::r_sexptype_needscast<T>() ) ;
}
	
}

/* as */

#define GENERATE_CONVERTERS(TYPE) template<> arma::Mat<TYPE> as< arma::Mat<TYPE> >(SEXP x){ return RcppArmadillo::convert<TYPE>(x, TYPE()) ; } ;

GENERATE_CONVERTERS(int)
GENERATE_CONVERTERS(arma::u32)
GENERATE_CONVERTERS(double)
GENERATE_CONVERTERS(float)

#undef GENERATE_CONVERTER

}

#endif

