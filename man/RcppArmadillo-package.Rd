\name{RcppArmadillo-package}
\alias{RcppArmadillo-package}
\alias{RcppArmadillo}
\alias{RcppArmadilloExample}
\docType{package}
\title{
Rcpp/Armadillo bridge
}
\description{
The package eases the integration of Armadillo types with Rcpp
}
\details{
\tabular{ll}{
Package: \tab RcppArmadillo\cr
Type: \tab Package\cr
Version: \tab 0.2.1\cr
Date: \tab 2010-05-19\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\section{Armadillo}{
\code{Armadillo} is a C++ linear algebra library (matrix maths) aiming towards a good
balance between speed and ease of use. Integer, floating point and complex numbers
are supported, as well as a subset of trigonometric and statistics functions. 
Various matrix decompositions are provided through optional integration with 
LAPACK and ATLAS libraries.

A delayed evaluation approach is employed (during compile time) to combine several 
operations into one and reduce (or eliminate) the need for temporaries. 
This is accomplished through recursive templates and template meta-programming.

This library is useful if C++ has been decided as the language of choice 
(due to speed and/or integration capabilities), rather than another language 
like Matlab or Octave. It is distributed under a license that is useful in 
both open-source and commercial contexts.

Armadillo is primarily developed at NICTA (Australia), 
with contributions from around the world.	
}

\section{RcppArmadillo}{
\code{RcppArmadillo} acts as a bridge between \code{Rcpp} and \code{Armadillo}, 
allowing the programmer to write code using armadillo classes that integrate
seemlessly with \code{Rcpp}. 
}

\section{Using RcppArmadillo}{
The simplest way to get started is to create a skeleton of a package
using \code{RcppArmadillo}, this is done by the \code{\link{RcppArmadillo.package.skeleton}}
function. 

The important steps are
\itemize{
	\item Include the \code{RcppArmadillo.h} header file, which also includes
	armadillo.h
	\item Depend and link to Rcpp and RcppArmadillo by adding these lines to 
	the DESCRIPTION file: 

	\preformatted{
	  Depends: Rcpp (>= 0.8.0), RcppArmadillo (>= 0.2.0)
	  LinkingTo: Rcpp, RcppArmadillo
	  }

	\item Link against the Rcpp, blas and lapack libraries, by adding this line
	in the Makevars

\preformatted{PKG_LIBS = $(shell $(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()" ) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) }

and this line to the Makevars.win: 
\preformatted{PKG_LIBS = $(shell $(R_HOME)/bin${R_ARCH_BIN}/Rscript.exe -e "Rcpp:::LdFlags()") $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) }
}
}

\section{Options: ATLAS and Boost}{
Armadillo can optionally use atlas and boost. Please refer to the 
Armadillo documentation \url{http://arma.sourceforge.net/}

\code{RcppArmadillo} does not force these features. it is up to the third party package 
to enable these them for their package. If your package wants to make use of 
these features, you can add a subset of 
these lines before \code{#include <RcppArmadillo.h>}

\preformatted{
#define ARMA_USE_BOOST
#define ARMA_USE_BOOST_DATE
#define ARMA_USE_ATLAS
}

when you do so, you must make sure that the corresponding features are 
available, which typically involves writing a configure file. Refer to 
Writing R Extensions for the procedure.
	
}

\section{Support}{
Please use the Rcpp-devel mailing list on r-forge 
for questions about RcppArmadillo (subscribe first).
\url{https://lists.r-forge.r-project.org/cgi-bin/mailman/listinfo/rcpp-devel}

Questions about armadillo itself should be directed to the armadillo forum
\url{http://sourceforge.net/apps/phpbb/arma/}	
}


\author{
Romain Francois, Dirk Eddelbuettel and Doug Bates
	
Maintainer: Romain, Dirk and Doug <RcppArmadillo-authors@r-enthusiasts.com>
}
\references{
Armadillo project: \url{http://arma.sourceforge.net/}
}
\keyword{ package }
\keyword{ programming }
\keyword{ interface }


