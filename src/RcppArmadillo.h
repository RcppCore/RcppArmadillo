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

/* wrap */

template <typename T> SEXP wrap ( const arma::Mat<T>& mat ){
	SimpleVector< traits::r_sexptype_traits<T>::rtype > vec( 
		Dimension( mat.n_rows, mat.n_cols ) );
	int n = mat.n_elem ;
	RCPPARMA_COPY(mat,vec,n)
	return vec ; 
} ;

template <typename T> SEXP wrap( const arma::Col<T>& column ){
	int n = column.n_rows ;
	SimpleVector< traits::r_sexptype_traits<T>::rtype > vec( 
		Dimension( n, 1 ) );
	RCPPARMA_COPY(column,vec,n)
	return vec ;
}

template <typename T> SEXP wrap( const arma::Row<T>& row ){
	int n = row.n_cols ;
	SimpleVector< traits::r_sexptype_traits<T>::rtype > vec( 
		Dimension( 1, n ) );
	RCPPARMA_COPY(row,vec,n)
	return vec ;
}

#ifdef HAS_CUBE
template <typename T> SEXP wrap( const arma::Cube<T>& cube ){
	int n = cube.n_elem;
	SimpleVector< traits::r_sexptype_traits<T>::rtype > vec( 
		Dimension( cube.n_rows, cube.n_cols, cube.n_slices ) );
	RCPPARMA_COPY(cube,vec,n)
	return vec ;
}
#endif

} // namespace Rcpp


namespace Rcpp{
	
namespace RcppArmadillo{

/* when a cast is needed */
template <typename T> ::arma::Mat<T> convert_Mat__dispatch( SEXP x, T, ::Rcpp::traits::true_type ){
	if( !Rf_isMatrix(x) ){
		throw std::range_error( "not a matrix" ) ;
	}
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
	::Rcpp::IntegerVector dim = input.attr("dim") ;
	::arma::Mat<T> out( dim[0], dim[1] ) ;
	int n = dim[0] * dim[1] ;
	RCPPARMA_COPY_CAST(input, out, n, T) ;
	return out;
}

/* when no cast is needed */
template <typename T> ::arma::Mat<T> convert_Mat__dispatch( SEXP x, T, ::Rcpp::traits::false_type ){
	if( !Rf_isMatrix(x) ){
		throw std::range_error( "not a matrix" ) ;
	}
	::Rcpp::SimpleVector< ::Rcpp::traits::r_sexptype_traits<T>::rtype > input(x);
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


}

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

