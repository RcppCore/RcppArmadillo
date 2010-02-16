#include "RcppArmadillo.h"

SEXP RcppArmadilloExample(){
	using namespace Rcpp ;
	using namespace arma ;
	
	List output(4); 
	
	imat x1 = eye<imat>( 3,3 ) ;
	x1.swap_rows(0,1)  ;
	output[0] = x1 ;
	
	mat x2 = eye<mat>( 3,3 ) ;
	x2.swap_rows(0,1)  ;
	output[1] = x2 ;
	
	colvec y1 = zeros<colvec>(5,1);
	output[2] = y1  ;
	
	rowvec y2 = zeros<mat>(1,5);
	output[3] = y2  ;
	
	return output ;
}

SEXP RcppArmadilloExample_as( SEXP input_ ){
	using namespace Rcpp ;
	using namespace arma ;
	
	List input(input_) ;
	imat m1 = input[0] ; /* implicit as */
	mat  m2 = input[1] ; /* implicit as */
	
	List res(2) ;
	res[0] = accu( m1 ) ;
	res[1] = accu( m2 ) ;
	
	return res ;
}
