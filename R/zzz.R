
RcppArmadilloExample <- function(){
	.Call( "RcppArmadilloExample", PACKAGE = "RcppArmadillo" )
}

RcppArmadilloExample_as <- function(){
	integer_mat <- matrix( as.integer(diag(4)), ncol = 4, nrow = 4 )
	numeric_mat <- diag(5)
	.Call( "RcppArmadilloExample_as", 
		list( integer_mat, numeric_mat ), 
		PACKAGE = "RcppArmadillo" )
}
