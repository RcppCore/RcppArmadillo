#ifndef RcppArmadillo__RcppArmadillo__h
#define RcppArmadillo__RcppArmadillo__h

#include <RcppCommon.h>
#include <armadillo>

/* forward declarations */
namespace Rcpp{
	template <typename T> SEXP wrap ( const arma::Mat<T>& ) ;
	template <typename T> SEXP wrap ( const arma::Row<T>& ) ;
	template <typename T> SEXP wrap ( const arma::Col<T>& ) ;
}

#include <Rcpp.h>

RcppExport SEXP RcppArmadilloExample() ;

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

/* as */

} // namespace Rcpp

#endif

