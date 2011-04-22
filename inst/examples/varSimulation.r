#!/usr/bin/r
##
## varSimulation.r: Simulation of first-order vector autoregression data
##
## Copyright (C)  2011  Lance Bachmeier and Dirk Eddelbuettel
##
## This file is part of RcppArmadillo.
##
## RcppArmadillo is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppArmadillo is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

suppressMessages(require(inline))

code <- '
  arma::mat coeff = Rcpp::as<arma::mat>(a);
  arma::mat errors = Rcpp::as<arma::mat>(e);
  int m = errors.n_rows; int n = errors.n_cols;
  arma::mat simdata(m,n);
  simdata.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    simdata.row(row) = simdata.row(row-1)*trans(coeff)+errors.row(row);
  }
  return Rcpp::wrap(simdata);
'

rcppSim <- cxxfunction(signature(a="numeric",e="numeric"),
                       code,plugin="RcppArmadillo")

a <- matrix(c(0.5,0.1,0.1,0.5),nrow=2)
e <- matrix(rnorm(10000),ncol=2)
rcppData <- rcppSim(a,e)

rSim <- function(coeff,errors) {
  simdata <- matrix(0,nrow(errors),ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
}
rData <- rSim(a,e)
stopifnot(all.equal(rcppData, rData))

suppressMessages(require(compiler))
compRsim <- cmpfun(rSim)

compRData <- compRsim(a,e)
stopifnot(all.equal(rcppData, compRData))

library(rbenchmark)
res <- benchmark(rcppSim(a,e),
                 rSim(a,e),
                 compRsim(a,e),
                 columns=c("test", "replications", "elapsed",
                           "relative", "user.self", "sys.self"),
                 order="relative")
print(res)
