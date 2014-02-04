\name{armadillo_version}
\alias{armadillo_version}
\title{Report the version of Armadillo}
\usage{
armadillo_version(single)
}
\arguments{
  \item{single}{A logical vector indicating whether a
  single return values is requested, or a named vector with
  three elements \code{major}, \code{minor} and
  \code{patch}.}
}
\value{
Depending on the value of \code{single}, either a single
number describing the Armadillo version or a named vector
with three elements \code{major}, \code{minor} and
\code{patch}.
}
\description{
Report the version of Armadillo
}
\details{
The version is defined by Armadillo in the header
\code{arma_version.hpp}.
}
\seealso{
Armadillo header file \code{arma_version.hpp}.
}

