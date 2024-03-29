% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{armadillo_throttle_cores}
\alias{armadillo_throttle_cores}
\alias{armadillo_reset_cores}
\title{Throttle (or Reset) (Rcpp)Armadillo to Two Cores}
\usage{
armadillo_throttle_cores(n = 2)

armadillo_reset_cores()
}
\arguments{
\item{n}{Integer value of desired cores, default is two}
}
\description{
Helper functions to throttle use of cores by RcppArmadillo-internal
code on systems with OpenMP. On package load, the initial value is
saved and used to reset the value.
}
