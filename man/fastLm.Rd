\name{fastLm}
\alias{fastLm}
\alias{fastLmPure}
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
fastLmPure(y, X)

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

  The \code{fastLmPure} function provides a reference use case of the \code{GSL}
  library via the wrapper functions in the \pkg{RcppGSL} package.
  
  The \code{fastLm} function provides a more standard implementation of
  a linear model fit, offering both a default and a formula interface as
  well as \code{print}, \code{summary} and \code{predict} methods.

  Lastly, one must be be careful in timing comparisons of
  \code{\link{lm}} and friends versus this approach based on
  \code{Armadillo}. The reason that \code{Armadillo} can do something
  like \code{\link{lm.fit}} faster than the functions in the stats
  package is because \code{Armadillo} uses the Lapack version of the QR
  decomposition while the stats package uses a \emph{modified} Linpack
  version.  Hence \code{Armadillo} uses level-3 BLAS code whereas the
  stats package uses level-1 BLAS.  However, \code{Armadillo} will
  either fail or, worse, produce completely incorrect answers
  on rank-deficient model matrices whereas the functions from the stats
  package will handle them properly due to the modified Linpack code.

  An example of the type of situation requiring extra care in checking
  for rank deficiency is a two-way layout with missing cells (see the
  examples section).  These cases require a special pivoting scheme of
  \dQuote{pivot only on (apparent) rank deficiency} which is not part of
  conventional linear algebra software.
}
\value{
  \code{fastLmPure} returns a list with three components:
  \item{coefficients}{a vector of coefficients}
  \item{stderr}{a vector of the (estimated) standard errors of the coefficient estimates}
  \item{df}{a scalar denoting the degrees of freedom in the model}
 
  \code{fastLm} returns a richer object which also includes the
  residuals, fitted values and call argument similar to the \code{\link{lm}} or
  \code{\link[MASS]{rlm}} functions..
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
  flm <- fastLmPure( log(trees$Volume),  cbind(1, log(trees$Girth)) )
  print(flm)

  ## standard R interface for formula or data returning object of class fastLm
  flmmod <- fastLm( log(Volume) ~ log(Girth), data=trees)
  summary(flmmod)

  ## case where fastLm breaks down
  dd <- data.frame(f1 = gl(4, 6, labels = LETTERS[1:4]),
                   f2 = gl(3, 2, labels = letters[1:3]))[-(7:8), ]
  xtabs(~ f2 + f1, dd)     # one missing cell
  mm <- model.matrix(~ f1 * f2, dd)
  kappa(mm)                # large, indicating rank deficiency
  set.seed(1)
  dd$y <- mm \%*\% seq_len(ncol(mm)) + rnorm(nrow(mm), sd = 0.1)
  summary(lm(y ~ f1 * f2, dd))     # detects rank deficiency
  summary(fastLm(y ~ f1 * f2, dd)) # some huge coefficients
}
\keyword{regression}




