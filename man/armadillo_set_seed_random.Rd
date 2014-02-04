\name{armadillo_set_seed_random}
\alias{armadillo_set_seed_random}
\title{Set the Armadillo Random Number Generator to a random value}
\usage{
armadillo_set_seed_random()
}
\value{
The function is invoked for its side effect and has no
return value.
}
\description{
Set the Armadillo Random Number Generator to a random value
}
\details{
Depending on whether RcppArmadillo was compiled for the
C++98 standard (currently the default) or for C++11
(optional), two different RNGs may be used. This function
resets either. For C++98, the system library RNG is used.
It may have poor performance, particularly on older
versions of Windows. For C++11, the RNG included in its
library is used.
}
\seealso{
The R documentation on its RNGs all of which are accessible
via \pkg{Rcpp}.
}

