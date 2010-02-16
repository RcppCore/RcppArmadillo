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
	
/* as */
template<> arma::Mat<int> as< arma::Mat<int> >(SEXP x){
	if( !Rf_isMatrix(x) ){
		throw std::range_error( "not a matrix" ) ;
	}
	SimpleVector< traits::r_sexptype_traits<int>::rtype > input(x);
	IntegerVector dim = input.attr("dim") ;
	arma::Mat<int> out( dim[0], dim[1] ) ;
	int n = dim[0] * dim[1] ;
	for( int i=0; i<n; i++){
		out[i] = input[i] ;
	}
	return out ;
}

/* as */
template<> arma::Mat<double> as< arma::Mat<double> >(SEXP x){
	if( !Rf_isMatrix(x) ){
		throw std::range_error( "not a matrix" ) ;
	}
	SimpleVector< traits::r_sexptype_traits<double>::rtype > input(x);
	IntegerVector dim = input.attr("dim") ;
	arma::Mat<double> out( dim[0], dim[1] ) ;
	int n = dim[0] * dim[1] ;
	for( int i=0; i<n; i++){
		out[i] = input[i] ;
	}
	return out ;
}

}

#endif

