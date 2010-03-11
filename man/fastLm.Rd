\name{fastLm}
\alias{fastLm}
\concept{regression}
\title{Bare-bones linear model fitting function}
\description{
  \code{fastLm} estimates the linear model using the \code{solve}
  function of \code{Armadillo} linear algebra library.
}
\usage{
fastLm(y, X)
}
\arguments{
  \item{y}{a vector containing the explained variable.}
  \item{X}{a data.frame or matrix containing the explanatory variables.}
}
\details{
  Linear models should be estimated using the \code{\link{lm}} function. In
  some cases, \code{\link{lm.fit}} may be appropriate.

  This function provides a reference use case of the \code{Armadillo}
  library.

  It does not provide a formula interface. It does check for missing
  values. It does (yet?) return a proper object with suitable accessor
  functions. That's why we call it bare-bones.
}
\value{
  \code{fastLm} returns a list with two components:
  \item{coefficients}{a vector of coefficients}
  \item{stderr}{a vector of the (estimated) standard errors of the coefficient estimates}

  In other words, no actual \sQuote{fastLm} object is provided unlike
  for \code{\link{lm}} or \code{\link[MASS]{rlm}} but rather just a
  simple list unlike \code{\link{lm.fit}}.
}
\seealso{\code{\link{lm}}, \code{\link{lm.fit}}} 
\references{Armadillo project: \url{http://arma.sourceforge.net/}}
\author{
  Armadillo is written by Conrad Sanderson. RcppArmadillo is written by
  Romain Francois and Dirk Eddelbuettel.
}
\examples{
  data(trees)
  flm <- fastLm( log(trees$Volume),  cbind(rep(1,31), log(trees$Girth)) )
}
\keyword{regression}

