#include "RcppArmadillo.h"

SEXP RcppArmadilloExample(){
	using namespace Rcpp ;
	using namespace arma ;
	
	List output(3); 
	List matrices(4) ;
	List rows(2) ;
	List cols(2) ;
	
	imat x1 = eye<imat>( 3,3 ) ;
	x1.swap_rows(0,1)  ;
	matrices[0] = x1 ;
	
	mat x2 = eye<mat>( 3,3 ) ;
	x2.swap_rows(0,1)  ;
	matrices[1] = x2 ;
	
	fmat x3 = eye<fmat>( 3, 3 );
	x3.swap_rows(0,1)  ;
	matrices[2] = x3 ;
	
	umat x4 = eye<umat>( 3, 3 );
	x4.swap_rows(0,1)  ;
	matrices[3] = x4 ;
	
	colvec r1 = zeros<mat>(5,1);
	cols[0] = r1  ;
	
	fcolvec r2 = zeros<fmat>(5,1);
	cols[1] = r2  ;
	
	rowvec c1 = zeros<mat>(1,5);
	rows[0] = c1  ;
	
	frowvec c2 = zeros<fmat>(1,5);
	rows[1] = c2  ;
	
	std::vector<std::string> names(4) ;
	names[0] = "Mat<int>" ;
	names[1] = "Mat<double>" ;
	names[2] = "Mat<float>" ;
	names[3] = "Mat<unsigned int>" ;
	matrices.names() = names ;
	
	names.resize(2) ;
	names[0] = "Row<double>" ;
	names[1] = "Row<float>" ;
	rows.names() = names ;
	
	names[0] = "Col<double>" ;
	names[1] = "Col<float>" ;
	cols.names() = names ;
	
	output[0] = matrices;
	output[1] = rows ;
	output[2] = cols ;
	
	names[0] = "matrices : Mat<T>" ; 
	names[1] = "rows : Row<T>" ;
	names.push_back( "columns : Col<T>" );
	output.names() = names ;
	
	return output ;
}

SEXP RcppArmadilloExample_as_Mat( SEXP input_ ){
	using namespace Rcpp ;
	
	List input(input_) ;
	arma::imat m1 = input[0] ; /* implicit as */
	arma::mat  m2 = input[1] ; /* implicit as */
	arma::umat m3 = input[0] ; /* implicit as */
	arma::fmat m4 = input[1] ; /* implicit as */
	
	List res(4) ;
	res[0] = arma::accu( m1 ) ;
	res[1] = arma::accu( m2 ) ;
	res[2] = arma::accu( m3 ) ;
	res[3] = arma::accu( m4 ) ;
	
	return res ;
}

SEXP RcppArmadilloExample_as_Col( SEXP input_ ){
	using namespace Rcpp ;
	
	List input(input_) ;
	arma::icolvec m1 = input[0] ; /* implicit as */
	arma::colvec  m2 = input[1] ; /* implicit as */
	arma::ucolvec m3 = input[0] ; /* implicit as */
	arma::fcolvec m4 = input[1] ; /* implicit as */
	
	List res(4) ;
	res[0] = arma::accu( m1 ) ;
	res[1] = arma::accu( m2 ) ;
	res[2] = arma::accu( m3 ) ;
	res[3] = arma::accu( m4 ) ;
	
	return res ;
}

SEXP RcppArmadilloExample_as_Row( SEXP input_ ){
	using namespace Rcpp ;
	
	List input(input_) ;
	arma::irowvec m1 = input[0] ; /* implicit as */
	arma::rowvec  m2 = input[1] ; /* implicit as */
	arma::urowvec m3 = input[0] ; /* implicit as */
	arma::frowvec m4 = input[1] ; /* implicit as */
	
	List res(4) ;
	res[0] = arma::accu( m1 ) ;
	res[1] = arma::accu( m2 ) ;
	res[2] = arma::accu( m3 ) ;
	res[3] = arma::accu( m4 ) ;
	
	return res ;
}
