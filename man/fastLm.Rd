\name{fastLm}
\alias{fastLm}
\alias{fastLmCall}
\alias{fastLm.default}
\alias{fastLm.formula}
\alias{print.fastLm}
\alias{summary.fastLm}
\alias{print.summary.fastLm}
\alias{predict.fastLm}
\concept{regression}
\title{Bare-bones linear model fitting function}
\description{
  \code{fastLm} estimates the linear model using the \code{solve}
  function of \code{Armadillo} linear algebra library.
}
\usage{
fastLmCall(y, X)

fastLm(x, ...)
\method{fastLm}{default}(x, y, ...)
\method{fastLm}{formula}(formula, data = list(), ...)
\method{print}{fastLm}(x, ...)
\method{summary}{fastLm}(object, ...)
\method{predict}{fastLm}(object, newdata=NULL, ...)
}
\arguments{
  \item{y}{a vector containing the explained variable.}
  \item{X}{a data.frame or matrix containing the explanatory variables.}
  \item{x}{a numeric design matrix for the model.}
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{newdata}{an optional data frame containing the variables for
    which a model prediction is desired.}
  \item{object}{an object of class \code{"fastLm"}, i.e. a fitted model}
  \item{\ldots}{not used}
}
\details{
  Linear models should be estimated using the \code{\link{lm}} function. In
  some cases, \code{\link{lm.fit}} may be appropriate.

  This function provides a reference use case of the \code{Armadillo}
  library.

  It does not provide a formula interface. It does check for missing
  values. It does (yet?) return a proper object with suitable accessor
  functions. That's why we call it bare-bones.

  Lastly, one must be be careful in timing comparisons of
  \code{\link{lm}} and friends versus this approach based on
  \code{Armadillo}. The reason that \code{Armadillo} can do something
  like \code{\link{lm.fit}} faster than the functions in the stats
  package is because \code{Armadillo} uses the Lapack version of the QR
  decomposition while the stats package uses a \emph{modified} Linpack
  version.  Hence \code{Armadillo} uses level-3 BLAS code whereas the
  stats package uses level-1 BLAS.  However, \code{Armadillo} will choke
  on rank-deficient model matrices whereas the functions from the stats
  package will handle them properly due to the modified Linpack code.
  Statisticians want a pivoting scheme of \dQuote{pivot only on
  (apparent) rank deficiency} and numerical analysts have no idea why
  statisticians want this so it is not part of conventional linear
  algebra software.
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
  Romain Francois, Dirk Eddelbuettel and Douglas Bates.
}
\examples{
  require(datasets)
  data(trees)

  ## bare-bones direct interface
  flm <- fastLmCall( log(trees$Volume),  cbind(1, log(trees$Girth)) )
  print(flm)

  ## standard R interface for formula or data returning object of class fastLm
  flmmod <- fastLm( log(Volume) ~ log(Girth), data=trees)
  summary(flmmod)
}
\keyword{regression}




