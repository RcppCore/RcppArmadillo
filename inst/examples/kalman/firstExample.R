
library(inline)

g <- cxxfunction(signature(vs="numeric"), plugin="RcppArmadillo", body='
     arma::colvec v = Rcpp::as<arma::colvec>(vs);
     arma::mat op = v * v.t();
     double ip = arma::as_scalar(v.t() * v);
     return Rcpp::List::create(Rcpp::Named("outer")=op,
                               Rcpp::Named("inner")=ip);
')

g(7:11)
