\name{RcppArmadillo.package.skeleton}
\alias{RcppArmadillo.package.skeleton}
\title{
Create a skeleton for a new package that intends to use RcppArmadillo
}
\description{
	\code{RcppArmadillo.package.skeleton} automates the creation of 
	a new source package that intends to use features of RcppArmadilo. 
	
	It is based on the \link[utils]{package.skeleton} function
	which it executes first.
}
\usage{
RcppArmadillo.package.skeleton(name = "anRpackage", list = character(), 
	environment = .GlobalEnv, path = ".", force = FALSE, namespace = TRUE, 
	code_files = character(), example_code = TRUE)
}
\arguments{
	\item{name}{See \link[utils]{package.skeleton}}
	\item{list}{See \link[utils]{package.skeleton}}
	\item{environment}{See \link[utils]{package.skeleton}}
	\item{path}{See \link[utils]{package.skeleton}}
	\item{force}{See \link[utils]{package.skeleton}}
	\item{namespace}{See \link[utils]{package.skeleton}}
	\item{code_files}{See \link[utils]{package.skeleton}}
	\item{example_code}{If TRUE, example c++ code using RcppArmadillo is added to the package}
}
\details{
	In addition to \link[utils]{package.skeleton} : 
	
	The \samp{DESCRIPTION} file gains a Depends line requesting that 
	the package depends on Rcpp and RcppArmadillo and 
	a LinkingTo line so that the package finds Rcpp and RcppArmadillo header files.
	
	The \samp{NAMESPACE}, if any, gains a \code{useDynLib} directive.
	
	The \samp{src} directory is created if it does not exists and 
	a \samp{Makevars} file is added setting the environment variable
	\samp{PKG_LIBS} to accomodate the necessary flags
	to link with the Rcpp library. 
	
	If the \code{example_code} argument is set to \code{TRUE}, 
	example files \samp{rcpparma_hello_world.h} and \samp{rcpparma_hello_world.cpp}
	are also created in the \samp{src}. An R file \samp{rcpparma_hello_world.R} is 
	expanded in the \samp{R} directory, the \code{rcpparma_hello_world} function
	defined in this files makes use of the C++ function \samp{rcpparma_hello_world}
	defined in the C++ file. These files are given as an example and should 
	eventually by removed from the generated package.
}
\value{
Nothing, used for its side effects
}
\seealso{
\link[utils]{package.skeleton}
}
\references{
Read the \emph{Writing R Extensions} manual for more details.

Once you have created a \emph{source} package you need to install it:
see the \emph{R Installation and Administration} manual,
\code{\link{INSTALL}} and \code{\link{install.packages}}.
}
\examples{
\dontrun{
RcppArmadillo.package.skeleton( "foobar" )
}
}
\keyword{ programming }

