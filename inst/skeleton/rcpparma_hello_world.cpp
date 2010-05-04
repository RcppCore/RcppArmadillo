#include "rcpparma_hello_world.h"

using namespace Rcpp ;

SEXP rcpparma_hello_world(){
	
	arma::mat m1 = arma::eye<arma::mat>( 3, 3 ) ;
	arma::mat m2 = arma::eye<arma::mat>( 3, 3 ) ;
	                     
	List res ;                                  
	res["mat+mat"] = m1 + 3 * ( m1 + m2 );
	return res ;
}

