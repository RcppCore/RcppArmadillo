\name{armadillo_set_seed}
\alias{armadillo_set_seed}
\title{Set the Armadillo Random Number Generator to the given value}
\usage{
armadillo_set_seed(val)
}
\arguments{
  \item{val}{The seed used to initialize Armadillo's random
  number generator.}
}
\value{
The function is invoked for its side effect and has no
return value.
}
\description{
Set the Armadillo Random Number Generator to the given
value
}
\details{
Depending on whether RcppArmadillo was compiled for the C++98 standard
(currently the default) or for C++11 (optional), two different RNGs may be used.
This function resets either. For C++98, the R programming language's RNG is used.
For C++11, the RNG included in the \code{<random>} library is used only when
 \code{#define ARMA_USE_CXX11_RNG} is placed before \code{#include <RcppArmadillo.h>}.
 Otherwise, the R programming language's RNG will be used.
}
\note{
This has been found to not work as espected in
\pkg{RStudio} as its code also uses the system RNG library.
You may have to either not run within \pkg{RStudio} or
change your code to use a different RNG such as the one
from R.
}
\seealso{
The R documentation on its RNGs all of which are accessible
via \pkg{Rcpp}.
}

