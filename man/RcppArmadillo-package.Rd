\name{RcppArmadillo-package}
\alias{RcppArmadillo-package}
\alias{RcppArmadillo}
\alias{RcppArmadilloExample}
\docType{package}
\title{R and Armadillo Integration}
\description{The package brings the power of Armadillo to R.}
\section{Armadillo}{
  \code{Armadillo} is a C++ linear algebra library, aiming towards a good
  balance between speed and ease of use.

  It provides efficient classes for vectors, matrices and cubes, as well
  as many functions which operate on the classes (eg. contiguous and
  non-contiguous submatrix views).

  Various matrix decompositions are provided, and an automatic
  expression evaluator (via template meta-programming) combines several
  operations to increase efficiency.

  The syntax (API) is deliberately similar to Matlab. It is useful for
  algorithm development directly in C++, or quick conversion of research
  code into production environments.

  Armadillo has been primarily developed at NICTA (Australia) by Conrad
  Sanderson, with contributions from around the world.
}
\section{RcppArmadillo}{
  \code{RcppArmadillo} acts as a bridge between \code{Rcpp} and \code{Armadillo},
  allowing the programmer to write code using Armadillo classes that integrate
  seemlessly with \code{R} via \code{Rcpp}.
}
\section{Using RcppArmadillo}{
  The simplest way to get started is to create a skeleton of a package
  using \code{RcppArmadillo}. This can be done conveniently by the
  \code{\link{RcppArmadillo.package.skeleton}}
  function.
}
\section{Threading}{
  The Armadillo library can take advantage of OpenMP to execute computations in
  parallel via multi-threaded code. The number of cores uses can be set (or
  retrieved) explicitly via helper functions \code{armadillo_get_number_of_omp+threads()} and
  \code{armadillo_set_number_of_omp_threads()}. A default value is stored at package startup; it
  recognises R option value \code{Ncpus} and environment variable
  \code{OMP_THREAD_LIMIT}.  Additional helper functions
  \code{armadillo_throttle_cores()} and \code{armadillo_reset_cores()} are
  available to (temporarily) lower the number of cores uses and to reset to
  the package default value set at startup.
}
\author{
  For RcppArmadillo: Dirk Eddelbuettel, Romain Francois, Doug Bates and
  Binxiang Ni

  Maintainer: Dirk Eddelbuettel <edd@debian.org>

  For Armadillo: Conrad Sanderson
}
\references{
  Armadillo project: \url{https://arma.sourceforge.net/}

  Conrad Sanderson and Ryan Curtin.
  \href{https://arma.sourceforge.net/armadillo_joss_2016.pdf}{Armadillo:
  a template-based C++ library for linear algebra}.
  Journal of Open Source Software, Vol. 1, pp. 26, 2016.

  Dirk Eddelbuettel and Conrad Sanderson,
  "RcppArmadillo: Accelerating R with high-performance C++ linear algebra",
  Computational	Statistics and Data Analysis, 2014, 71, March, pages
  1054-1063, \doi{10.1016/j.csda.2013.02.005}.
}
\keyword{package}
\keyword{programming}
\keyword{interface}
