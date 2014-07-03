\name{NEWS}
\title{News for Package 'RcppArmadillo'}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in RcppArmadillo version 0.4.320.0 (2014-07-03)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.320 (Daintree Tea Raider)
    \itemize{
      \item expanded \code{eigs_sym()} and \code{eigs_gen()} to use an
      optional tolerance parameter
      \item expanded \code{eig_sym()} to automatically fall back to standard
      decomposition method if divide-and-conquer fails
      \item automatic installer enables use of C++11 random number generator
      when using gcc 4.8.3+ in C++11 mode
    }
  }
}
      
\section{Changes in RcppArmadillo version 0.4.300.8.0 (2014-05-31)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.300.8 (Medieval Cornea Scraper)
    \itemize{
      \item More robust \code{norm}-related functions
      \item Fixes between interactions between \code{cube} and
      \code{vector} types.
    }
    \item Adds a \code{#define ARMA_DONT_USE_CXX11} to provide an option
    to turn C++11 off for Armadillo (but client packages may still use it)
    \item More robust Windows detection by using \code{_WIN32} as
    well as \code{WIN32} as the latter gets diabled by MinGW with C++11
    \item On Windows, C++11 is turned off as the Armadillo code base uses
    more features of C++11 than g++ 4.6.2 version in Rtools implements
  }
}

\section{Changes in RcppArmadillo version 0.4.300.5.0 (2014-05-19)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.300.5 (Medieval Cornea Scraper)
    \itemize{
      \item Handle possible underflows and overflows in \code{norm()},
      \code{normalise()}, \code{norm_dot()}
      \item Fix for handling of null vectors by \code{norm_dot()}
    }
  }
}

\section{Changes in RcppArmadillo version 0.4.300.2.0 (2014-05-13)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.300.2 (Medieval Cornea Scraper)
    \itemize{
      \item faster \code{find()}
    }
  }
}

\section{Changes in RcppArmadillo version 0.4.300.0 (2014-05-04)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.300 (Medieval Cornea Scraper)
    \itemize{
      \item faster \code{find()}
      \item added \code{find_finite()} and \code{find_nonfinite()} for
      finding indices of finite and non-finite elements 
      \item expressions \code{X=inv(A)*B*C} and \code{X=A.i()*B*C} are
      automatically converted to X=solve(A,B*C) 
    }
    \item Corrected conversion to \code{unsigned int} vectors and
    matrices
    \item Configure script now checks for R version 3.0.3 or newer to
    enable complex divide-and-conquer SVD in case of R-supplied LAPACK
  }
}
  
\section{Changes in RcppArmadillo version 0.4.200.0 (2014-04-07)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.200 (Flintlock Swoop)
    \itemize{
      \item faster transpose of sparse matrices
      \item more efficient handling of aliasing during matrix multiplication
      \item faster inverse of matrices marked as diagonal
    }
  }
}

\section{Changes in RcppArmadillo version 0.4.100.2 (2014-03-07)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.100.2 
    \itemize{
      \item fix for handling null vectors by \code{normalise()}
      \item fix for memory handling by sparse matrices
    }
    \item Correct use of \code{[[ depends()]]} in skeleton example file
    \item Prepare \code{src/Makevars} for C++11 support from R 3.1.0 by
    defining \code{USE_CXX11} which is currently commented out
    \item In the Armadillo configuration, turn on C++11 support if
    \code{USE_CXX11} is defined
  }
}

\section{Changes in RcppArmadillo version 0.4.100.0 (2014-02-28)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.100.0 (Dirt Cruiser)
    \itemize{
      \item added \code{normalise()} for normalising vectors to unit p-norm
      \item extended the \code{field} class to handle 3D layout
      \item extended \code{eigs_sym()} and \code{eigs_gen()} to obtain
      eigenvalues of various forms (eg. largest or smallest magnitude) 
      \item automatic SIMD vectorisation of elementary expressions
      (eg. matrix addition) when using Clang 3.4+ with -O3 optimisation 
      \item faster handling of sparse submatrix views
      \item workaround for a bug in LAPACK 3.4
    }
  }
}

\section{Changes in RcppArmadillo version 0.4.000.4 (2014-02-19)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.000.4
    \itemize{
      \item fix for \code{randi()} generating out-of-interval values
      \item fix for saving field objects
      \item workaround for a bug in the Intel compiler
    }
    \item Updated for \cpkg{Rcpp} (>= 0.11.0) by removing linking step
    from build process, added appropriate dependency on \cpkg{Rcpp}
    \item Updated \code{RcppArmadillo.package.skeleton} function
    accordingly to use proper \code{NAMESPACE} import
    \item Rewritten \code{rcpparma_hello_world} (which is used by the
    \code{RcppArmadillo.package.skeleton} function) to use Rcpp
    Attributes, and added more examples
    \item Added two functions to set Armadillo's RNGs (ie the system
    RNG) from a given value and to randomize it, as suggested by Gábor
    Csárdi. Note that these do not work within RStudio (which itself
    also uses the same system RNG). 
  }
}

\section{Changes in RcppArmadillo version 0.4.000.2 (2014-01-21)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.000.2
    \itemize{
      \item fix for \code{randi()} generating out-of-interval values
      \item workaround for a bug in the Intel compiler
    }
  }
}

\section{Changes in RcppArmadillo version 0.4.000 (2014-01-05)}{
  \itemize{
    \item Upgraded to Armadillo release Version 4.000 (Feral Steamroller)
    \itemize{
      \item added eigen decompositions of sparse matrices:
      \code{eigs_sym()} and \code{eigs_gen()} [ but this requires
      linking against ARPACK which \cpkg{RcppArmadillo} as a
      pure-template package does not do, and R is not linked against
      ARPACK either. ]
      \item added eigen decomposition for pair of matrices: \code{eig_pair()}
      \item added simpler forms of \code{eig_gen()}
      \item added condition number of matrices: \code{cond()}
      \item expanded \code{find()} to handle cubes
      \item expanded subcube views to access elements specified in a vector
      \item template argument for \code{running_stat_vec} expanded to
      accept vector types 
      \item more robust fast inverse of 4x4 matrices
      \item faster divide-and-conquer decompositions are now used by
      default for \code{eig_sym()}, \code{pinv()}, \code{princomp()},
      \code{rank()}, \code{svd()}, \code{svd_econ()} 
      \item the form \code{inv(sympd(X))} no longer assumes that X is
      positive definite; use \code{inv_sympd()} instead 
      \item added MEX connector for interfacing Octave/Matlab with
      Armadillo matrices (contributed by George Yammine) 
    }
  }
}

\section{Changes in RcppArmadillo version 0.3.930.1 (2013-12-09)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.930.1
    \itemize{
      \item Armadillo falls back to standard complex svd if the more
      performant divide-and-conquer variant is unavailable
    }
    \item Added detection for Lapack library and distinguish between R's
    own version (withhout zgesdd) and system Lapack; a preprocessor define is
    set accordingly
  }
}

\section{Changes in RcppArmadillo version 0.3.930.0 (2013-12-06)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.930 ("Dragon's Back")
    \itemize{
      \item added divide-and-conquer variant of \code{svd_econ()}, for
      faster SVD 
      \item added divide-and-conquer variant of \code{pinv()}, for
      faster pseudo-inverse 
      \item added element-wise variants of \code{min()} and \code{max()} 
      \item added \code{size()} based specifications of submatrix view
      sizes 
      \item added \code{randi()} for generating matrices with random
      integer values 
      \item added more intuitive specification of sort direction in
      \code{sort()} and \code{sort_index()} 
      \item added more intuitive specification of method in
      \code{det()}, \code{.i()}, \code{inv()} and \code{solve()} 
      \item added more precise timer for the \code{wall_clock} class
      when using C++11 
    }
    \item New unit tests for complex matrices and vectors
  }
}

\section{Changes in RcppArmadillo version 0.3.920.3 (2013-11-20)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.920.3
    \itemize{
      \item fix for handling of tiny matrices by \code{.swap()}
    }
  }
}

\section{Changes in RcppArmadillo version 0.3.920.1 (2013-09-27)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.920.1 (Agencia
    Nacional Stasi)
    \itemize{
      \item faster \code{.zeros()}
      \item faster \code{round()}, \code{exp2()} and \code{log2()} when using C++11
      \item added signum function: \code{sign()}
      \item added move constructors when using C++11
      \item added 2D fast Fourier transform: \code{fft2()}
      \item added \code{.tube()} for easier extraction of vectors and subcubes from
      cubes
      \item added specification of a fill type during construction of Mat, Col,
      Row and Cube classes, eg. \code{mat X(4, 5, fill::zeros)}
    }
    \item Initial implementation of \code{wrap<subview>}
    \item Improved implementation of \code{as<>()} and \code{wrap()}
    for sparse matrices
    \item Converted main vignette from \code{LaTeX} style \code{minted}
    to \code{lstlisting} which permits builds on CRAN; removed
    set \code{BuildVignettes: FALSE}.
  }
}

\section{Changes in RcppArmadillo version 0.3.910.0 (2013-08-12)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.910.0 (Pyrenees)
    \itemize{
      \item faster multiplication of a matrix with a transpose of
      itself, ie. \code{X*X.t()} and \code{X.t()*X} 
      \item added \code{vectorise()} for reshaping matrices into vectors
      \item added \code{all()} and \code{any()} for indicating presence
      of elements satisfying a relational condition 
    }
    \item Added conversion support for sparse matrices (of type double)
    created by the \cpkg{Matrix} package as class \code{dgCMatrix}
    \item Moved vignette sources from \code{inst/doc} to \code{vignettes};
    set \code{BuildVignettes: FALSE} as the \code{minted} mode for
    \code{LaTeX} upsets the CRAN builders.
  }
}

\section{Changes in RcppArmadillo version 0.3.900.7 (2013-08-02)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.900.7 (Bavarian
    Sunflower)
    \itemize{
      \item minor fix for inplace \code{reshape()}
      \item minor corrections for compilation issues under GCC 4.8+ and MSVC
    }
    \item Corrected setting of \code{vec_stats} in intialization of row,
    vector and matrix objects
    \item The \pkg{inline} package is no longer used in the examples and
    unit tests which have all been converted to using Rcpp attributes 
  }
}

\section{Changes in RcppArmadillo version 0.3.900 (2013-06-04)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.900.0 (Bavarian
    Sunflower)
    \itemize{
      \item added automatic SSE2 vectorisation of elementary expressions
      (eg. matrix addition) when using GCC 4.7+ with -O3 optimisation
      \item added support for saving & loading of cubes in HDF5 format,
      contributed by Szabolcs Horvat
      \item faster \code{median()}, contributed by Ruslan Shestopalyuk
      \item faster handling of compound expressions with transposes of
      submatrix rows
      \item faster handling of compound expressions with transposes of
      complex vectors
    }
    \item Kalman filter example switched from inline to \code{sourceCpp},
    which simplifies / streamlines the C++ side a little
  }
}

\section{Changes in RcppArmadillo version 0.3.820 (2013-05-12)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.820 (Mt Cootha)
    \itemize{
      \item faster \code{as_scalar()} for compound expressions
      \item faster transpose of small vectors
      \item faster matrix-vector product for small vectors
      \item faster multiplication of small fixed size matrices
    }
  }
}

\section{Changes in RcppArmadillo version 0.3.810.2 (2013-04-30)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.810.2
    \itemize{
      \item minor fix for initialisation of sparse matrices
    }
  }
}

\section{Changes in RcppArmadillo version 0.3.810.0 (2013-04-19)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.810.0 (Newell Highway)
    \itemize{
      \item added fast Fourier transform: \code{fft()}
      \item added handling of \code{.imbue()} and \code{.transform()} by
      submatrices and subcubes 
      \item added batch insertion constructors for sparse matrices
      \item minor fix for multiplication of complex sparse matrices
    }
    \item Updated sample() function and test again contributed by
    Christian Gunning 
  }
}

\section{Changes in RcppArmadillo version 0.3.800.1 (2013-03-12)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.800.1 (Miami Beach)
    \itemize{
      \item workaround for a bug in ATLAS 3.8 on 64 bit systems
      \item faster matrix-vector multiply for small matrices
    }
    \item Added new sample() function and tests contributed by Christian Gunning
    \item Refactored unit testing code for faster unit test performance
  }
}

\section{Changes in RcppArmadillo version 0.3.800.0 (2013-03-01)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.800.0 (Miami Beach)
    \itemize{
      \item Armadillo is now licensed using the Mozilla Public License 2.0
      \item added \code{.imbue()} for filling a matrix/cube with values provided by a functor or lambda expression
      \item added \code{.swap()} for swapping contents with another matrix
      \item added \code{.transform()} for transforming a matrix/cube using a functor or lambda expression
      \item added \code{round()} for rounding matrix elements towards nearest integer
      \item faster \code{find()}
      \item fixes for handling non-square matrices by \code{qr()} and \code{qr_econ()}
      \item minor fixes for handling empty matrices
      \item reduction of pedantic compiler warnings
    }
    \item Updated vignette to paper now in press at CSDA
    \item Added CITATION file with reference to CSDA paper
  }
}

\section{Changes in RcppArmadillo version 0.3.6.3 (2013-02-20)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.6.3
    \itemize{
      \item faster \code{find()}
      \item minor fix for non-contiguous submatrix views to handle empty vectors of indices
      \item reduction of pedantic compiler warnings
    }
  }
}

\section{Changes in RcppArmadillo version 0.3.6.2 (2013-01-29)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.6.2
    \itemize{
      \item faster determinant for matrices marked as diagonal or triangular
      \item more fine-grained handling of 64 bit integers
    }
    \item Added a new example of a Kalman filter implementation in R, and C++
    using Armadillo via RcppArmadillo, complete with timing comparison
  }
}

\section{Changes in RcppArmadillo version 0.3.6.1 (2012-12-17)}{
  \itemize{
    \item Upgraded to Armadillo release Version 3.6.1 (Piazza del Duomo)
    \itemize{
      \item faster \code{trace()}
      \item fix for handling sparse matrices by \code{dot()}
      \item fixes for interactions between sparse and dense matrices
    }
    \item Now throws compiler error if \code{Rcpp.h} is included before
    \code{RcppArmadillo.h} (as the former is included automatically by the
    latter anyway, but template logic prefers this ordering).
  }
}

\section{Changes in RcppArmadillo version 0.3.4.3 (2012-10-04)}{
  \itemize{
    \item Upgraded to Armadillo release 3.4.3
    \itemize{
      \item fix for aliasing issue in \code{diagmat()}
      \item fix for \code{speye()} signature
    }
  }
}

\section{Changes in RcppArmadillo version 0.3.4.2 (2012-09-25)}{
  \itemize{
    \item Upgraded to Armadillo release 3.4.2
    \itemize{
      \item minor fixes for handling sparse submatrix views
      \item minor speedups for sparse matrices
    }
  }
}

\section{Changes in RcppArmadillo version 0.3.4.1 (2012-09-18)}{
  \itemize{
    \item Upgraded to Armadillo release 3.4.1
    \itemize{
      \item workaround for a bug in the Mac OS X accelerate framework
      \item fixes for handling empty sparse matrices
      \item added documentation for saving & loading matrices in HDF5 format
      \item faster dot() and cdot() for complex numbers
    }
  }
}

\section{Changes in RcppArmadillo version 0.3.4.0 (2012-09-06)}{
  \itemize{
    \item Upgraded to Armadillo release 3.4.0 (Ku De Ta)
    \itemize{
      \item added economical QR decomposition: qr_econ()
      \item added .each_col() & .each_row() for vector operations repeated on each column or row
      \item added preliminary support for sparse matrices, contributed by Ryan Curtin et al. (Georgia Institute of Technology)
      \item faster singular value decomposition via divide-and-conquer algorithm
      \item faster .randn()
    }
    \item NEWS file converted to Rd format
  }
}
\section{Changes in RcppArmadillo version 0.3.3.91 (2012-08-30)}{
  \itemize{
    \item Upgraded to Armadillo release 3.3.91
    \itemize{
      \item faster singular value decomposition via "divide and conquer" algorithm
      \item added economical QR decomposition: qr_econ()
      \item added .each_col() & .each_row() for vector operations repeated on each column or row
      \item added preliminary support for sparse matrices, contributed by Ryan Curtin, James Cline and Matthew Amidon (Georgia Institute of Technology)
    }
    \item Corrected summary method to deal with the no intercept case when using a formula; also display residual summary() statistics
    \item Expanded unit tests for fastLm
  }
}
\section{Changes in RcppArmadillo version 0.3.2.4 (2012-07-11)}{
  \itemize{
    \item Upgraded to Armadillo release 3.2.4
    \itemize{
      \item workaround for a regression (bug) in GCC 4.7.0 and 4.7.1
    }
  }
}
\section{Changes in RcppArmadillo version 0.3.2.3 (2012-07-01)}{
  \itemize{
    \item Upgraded to Armadillo release 3.2.3 
    \itemize{
      \item minor correction for declaration of fixed size vectors and matrices
      \item Reverted three header files \{Mat,Row,Col\}_bones.hpp back to previous release due to compilation failures under g++-4.7
      \item Added new vignette 'RcppArmadillo-intro' based on a just-submitted introductory paper (by Eddelbuettel and Sanderson) about RcppArmadillo 
      \item Change from release 3.2.2 which we skipped as it did not really affect builds under R:
      \itemize{
	\item minor fix for compiling without debugging enabled (aka release mode)
	\item better detection of ATLAS during installation on Fedora and Red Hat systems
      }
      \item Small enhancement to fastLm 
    }
  }
}
\section{Changes in RcppArmadillo version 0.3.2.0 (2012-05-21)}{
  \itemize{
    \item Upgraded to Armadillo release 3.2.0 "Creamfields"
    \itemize{
      \item faster eigen decomposition via "divide and conquer" algorithm
      \item faster transpose of vectors and compound expressions
      \item faster handling of diagonal views
      \item faster handling of tiny fixed size vectors (≤ 4 elements)
      \item added unique(), for finding unique elements of a matrix
    }
  }
}
\section{Changes in RcppArmadillo version 0.3.1.94 (2012-05-15)}{
  \itemize{
    \item Upgraded to Armadillo release 3.1.94 "v3.2 beta 2"
    \itemize{
      \item added unique(), for finding unique elements of a matrix
      \item faster eigen decomposition via "divide and conquer" algorithm
      \item faster transpose of vectors and compound expressions
      \item faster handling of tiny fixed size vectors (≤ 4 elements)
    }
  }
}
\section{Changes in RcppArmadillo version 0.3.1.92 (2012-05-10)}{
  \itemize{
    \item Upgraded to Armadillo release 3.1.92 "v3.2 beta 2"
    \itemize{
      \item added unique(), for finding unique elements of a matrix
      \item faster eigen decomposition via optional use of "divide and conquer" by eig_sym() 
      \item faster transpose of vectors and compound expressions
    }
  }
}
\section{Changes in RcppArmadillo version 0.3.0.3 (2012-05-03)}{
  \itemize{
    \item Upgraded to Armadillo release 3.0.3
    \itemize{
      \item fixes for inplace transpose of complex number matrices
      \item fixes for complex number version of svd_econ()
      \item fixes for potential aliasing issues with submatrix views
    }
    \item New example script fastLm 
  }
}
\section{Changes in RcppArmadillo version 0.3.0.2 (2012-04-19)}{
  \itemize{
    \item Upgraded to Armadillo release 3.0.2
    \itemize{
      \item fixes for handling diagonal matrices
    }
    \item Undefine NDEBUG if it has been set (as R does) as this prevents a number of useful debugging checks. Users can still define it or define ARMA_NO_DEBUG if they want a 'non-development' build
  }
}
\section{Changes in RcppArmadillo version 0.3.0.1 (2012-04-12)}{
  \itemize{
    \item Upgraded to Armadillo release 3.0.1
    \itemize{
      \item fixes for compilation errors
      \item fixes for potential aliasing issues
    }
  }
}
\section{Changes in RcppArmadillo version 0.3.0 (2012-04-10)}{
  \itemize{
    \item Upgraded to Armadillo release 3.0.0 "Antarctic Chilli Ranch"
    \itemize{
      \item added non-contiguous submatrix views
      \item added shorthand for inverse: .i()
      \item added hist() and histc()
      \item faster repmat()
      \item faster handling of submatrix views with one row or column 
      \item faster generation of random numbers
      \item faster element access in fixed size matrices
      \item better detection of vector expressions by sum(), cumsum(), prod(), min(), max(), mean(), median(), stddev(), var() 
      \item expressions X=A.i()*B and X=inv(A)*B are automatically converted to X=solve(A,B) 
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.40 (2012-04-04)}{
  \itemize{
    \item Upgraded to Armadillo release 2.99.4 "Antarctic Chilli Ranch (Beta 4)"
    \itemize{
      \item fixes for handling expressions with fixed size matrices
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.39 (2012-04-02)}{
  \itemize{
    \item Upgraded to Armadillo release 2.99.3 "Antarctic Chilli Ranch (Beta 3)"
    \itemize{
      \item faster repmat()
      \item workarounds for braindead compilers (eg. Visual Studio)
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.38 (2012-03-28)}{
  \itemize{
    \item Upgraded to Armadillo release 2.99.2 "Antarctic Chilli Ranch (Beta 2)"
    \itemize{
      \item added .i()
      \item much faster handling of .col() and .row()
      \item expressions X=A.i()*B and X=inv(A)*B are automatically converted to X=solve(A,B) 
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.37 (2012-03-19)}{
  \itemize{
    \item Upgraded to Armadillo release 2.99.1 "Antarctic Chilli Ranch (Beta 1)"
    \itemize{
      \item added non-contiguous submatrix views
      \item added hist() and histc()
      \item faster handling of submatrix views
      \item faster generation of random numbers
      \item faster element access in fixed size matrices
      \item better detection of vector expressions by sum(), cumsum(), prod(), min(), max(), mean(), median(), stddev(), var() 
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.36 (2012-03-05)}{
  \itemize{
    \item Upgraded to Armadillo release 2.4.4
    \itemize{
      \item fixes for qr() and syl()
      \item more portable wall_clock class
      \item faster relational operators on submatrices
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.35 (2012-02-17)}{
  \itemize{
    \item Upgraded to Armadillo release 2.4.3
    \itemize{
      \item Support for ARMA_DEFAULT_OSTREAM using Rcpp::Rcout added
      \item Minor bug fix release improving corner cases affecting builds:
      \itemize{
	\item Missing semicolon added in Mat_meat (when in C++0x mode), with thanks to Teo Guo Ci 
	\item Armadillo version vars now instantiated in RcppArmadillo.cpp which helps older g++ versions, with thanks to Gershon Bialer
	\item Thanks also to Martin Renner for testing these changes
	\item Unit tests output fallback directory changed per Brian Ripley's request to not ever use /tmp
	\item Minor update to version numbers in RcppArmadillo-package.Rd
      }
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.34 (2011-12-12)}{
  \itemize{
    \item Upgraded to Armadillo release 2.4.2
    \itemize{
      \item clarified documentation for .reshape()
      \item fix for handling of empty matrices by .resize()
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.33 (2011-12-07)}{
  \itemize{
    \item Upgraded to Armadillo release 2.4.1
    \itemize{
      \item added .resize()
      \item fix for vector initialisation
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.32 (2011-12-04)}{
  \itemize{
    \item Upgraded to Armadillo test release 2.4.0 "Loco Lounge Lizard"
    \item Minimal changes relative to 0.2.31 based on 2.3.92, next section is relative to the previous stable release series 2.2.* of Armadillo
    \itemize{
      \item added shorter forms of transposes: .t() and .st()
      \item added optional use of 64 bit indices, allowing matrices to have more than 4 billion elements 
      \item added experimental support for C++11 initialiser lists
      \item faster pinv()
      \item faster inplace transpose
      \item faster handling of expressions with diagonal views
      \item fixes for handling expressions with aliasing and submatrices
      \item fixes for linking on Ubuntu and Debian systems
      \item fixes for inconsistencies in interactions between matrices and cubes 
      \item refactored code to eliminate warnings when using the Clang C++ compiler 
      \item .print_trans() and .raw_print_trans() are deprecated
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.31 (2011-11-28)}{
  \itemize{
    \item Upgraded to Armadillo test release 2.3.92 "Loco Lounge Lizard (Beta 2)"
    \itemize{
      \item fixes for linking on Ubuntu and Debian systems
      \item fixes for inconsistencies in interactions between matrices and cubes
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.30 (2011-11-19)}{
  \itemize{
    \item Upgraded to Armadillo test release 2.3.91 "Loco Lounge Lizard (Beta 1)"
    \itemize{
      \item added shorter forms of transposes: .t() and .st()
      \item added optional use of 64 bit indices, allowing matrices to have more than 4 billion elements
      \item added experimental support for C++11 initialiser lists
      \item faster pinv()
      \item faster inplace transpose
      \item bugfixes for handling expressions with aliasing and submatrices
      \item refactored code to eliminate warnings when using the Clang C++ compiler
      \item .print_trans() and .raw_print_trans() are deprecated 
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.29 (2011-09-01)}{
  \itemize{
    \item Upgraded to Armadillo release 2.2.3
    \itemize{
      \item Release fixes a speed issue in the as_scalar() function.
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.28 (2011-08-02)}{
  \itemize{
    \item Upgraded to Armadillo release 2.2.1 "Blue Skies Debauchery"
    \itemize{
      \item faster multiplication of small matrices
      \item faster trans()
      \item faster handling of submatrices by norm()
      \item added economical singular value decomposition: svd_thin()
      \item added circ_toeplitz()
      \item added .is_colvec() & .is_rowvec()
      \item fixes for handling of complex numbers by cov(), cor(), running_stat_vec
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.27 (2011-07-22)}{
  \itemize{
    \item Upgraded to Armadillo release 2.1.91 "v2.2 beta 1"
    \itemize{
      \item faster multiplication of small matrices
      \item faster trans()
      \item faster handling of submatrices by norm()
      \item added economical singular value decomposition: svd_thin()
      \item added circ_toeplitz()
      \item added .is_colvec() & .is_rowvec()
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.26 (2011-07-17)}{
  \itemize{
    \item Upgraded to Armadillo release 2.0.2 
    \itemize{
      \item fix for handling of conjugate transpose by as_scalar()
      \item fix for handling of aliasing by diagmat()
      \item fix for handling of empty matrices by symmatu()/symmatl()
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.25 (2011-06-30)}{
  \itemize{
    \item Upgraded to Armadillo 2.0.1 which fixes two minor compilation issues
  }
}
\section{Changes in RcppArmadillo version 0.2.24 (2011-06-29)}{
  \itemize{
    \item Upgraded to Armadillo release 2.0.0 "Carnivorous Sugar Glider"
    \itemize{
      \item faster multiplication of tiny matrices (≤ 4x4)
      \item faster compound expressions containing submatrices
      \item faster inverse of symmetric positive definite matrices
      \item faster element access for fixed size matrices
      \item added handling of arbitrarily sized empty matrices (eg. 5x0)
      \item added loading & saving of matrices as CSV text files
      \item added .count() member function to running_stat and running_stat_vec 
      \item added syl(), strans(), symmatu()/symmatl()
      \item added submatrices of submatrices
      \item det(), inv() and solve() can be forced to use more precise
      \item algorithms for tiny matrices (≤ 4x4)
      \item htrans() has been deprecated; use trans() instead
      \item API change: trans() now takes the complex conjugate when transposing a complex matrix 
      \item API change: .is_vec() now outputs true for empty vectors (eg. 0x1) 
      \item API change: forms of chol(), eig_sym(), eig_gen(), inv(), lu(), pinv(), princomp(), qr(), solve(), svd(), syl() that do not return a bool indicating success now throw std::runtime_error exceptions when failures are detected 
      \item API change: princomp_cov() has been removed; princomp() in conjunction with cov() can be used instead
      \item API change: set_log_stream() & get_log_stream() have been replaced by set_stream_err1() & get_stream_err1()
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.23 (2011-06-23)}{
  \itemize{
    \item Upgraded to Armadillo release 1.99.5 "v2.0 beta 5"
    \itemize{
      \item Forms of chol(), eig_sym(), eig_gen(), inv(), lu(), pinv(), princomp(), qr(), solve(), svd(), syl() that do not return a bool indicating success now throw std::runtime_error exceptions when failures are detected 
      \item princomp_cov() has been removed; princomp() in conjunction with cov() can be used instead
      \item set_log_stream() & get_log_stream() have been replaced by set_stream_err1() & get_stream_err1()
      \item det(), inv() and solve() can be forced to use more precise algorithms for tiny matrices (≤ 4x4)
      \item Added loading & saving of matrices as CSV text files 
    }
    \item fastLmPure() now uses same argument order as R's lm.fit()
    \item Export and document S3 methods in NAMESPACE and manual page as such
  }
}
\section{Changes in RcppArmadillo version 0.2.22 (2011-06-06)}{
  \itemize{
    \item Upgraded to Armadillo release 1.99.4 "v2.0 beta 4" 
    \itemize{
      \item fixes for handling of tiny matrices
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.21 (2011-05-27)}{
  \itemize{
    \item Upgraded to Armadillo release 1.99.3 "v2.0 beta 3" 
    \itemize{
      \item stricter size checking for row and column vectors
      \item added .count() member function to running_stat and running_stat_vec
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.20 (2011-05-25)}{
  \itemize{
    \item Upgraded to Armadillo release 1.99.2 "v2.0 beta 2" (and 1.99.1 before)
    \itemize{
      \item faster inverse of symmetric matrices
      \item faster element access for fixed size matrices
      \item faster multiplication of tiny matrices (eg. 4x4)
      \item faster compund expressions containing submatrices
      \item added handling of arbitrarily sized empty matrices (eg. 5x0)
      \item added syl()
      \item added strans()
      \item added symmatu()/symmatl()
      \item added submatrices of submatrices
      \item htrans() has been deprecated; use trans() instead
      \item trans() now takes the complex conjugate when transposing a complex matrix
      \item .is_vec() now outputs true for empty matrices
      \item most functions with matrix inputs no longer throw exceptions when given empty matrices (eg. 5x0)
    }
    \item Added a new subdirectory examples/ seeded with a nice Vector Autoregression simulation simulation example by Lance Bachmeier
    \item Rewrote armadillo_version as to no longer require an instance of arma::arma_version, with tanks to Conrad for the suggestion
  }
}
\section{Changes in RcppArmadillo version 0.2.19 (2011-04-18)}{
  \itemize{
    \item Upgraded to Armadillo version 1.2.0 "Unscrupulous Carbon Emitter"
    \itemize{
      \item Added ability to use Blas & Lapack libraries with capitalised function names 
      \item Reduction of pedantic compiler warnings
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.18 (2011-04-03)}{
  \itemize{
    \item Upgraded to Armadillo version 1.1.92  "Jurassic Barbecue"
    \itemize{
      \item Bugfix in cor()
      \item Automatic installation now requires CMake >= 2.6
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.17 (2011-03-22)}{
  \itemize{
    \item Upgraded to Armadillo version 1.1.90  "Inside Job"
    \itemize{
      \item Added .min() & .max(), which can provide the extremum's location
      \item More robust mean(), var(), stddev()
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.16 (2011-03-10)}{
  \itemize{
    \item Upgraded to Armadillo version 1.1.8   "Kangaroo Steak"
    \itemize{
      \item Added floor() and ceil()
      \item Added “not a number”: math::nan()
      \item Added infinity: math::inf()
      \item Added standalone is_finite()
      \item Faster min(), max(), mean()
      \item Bugfix for a corner case with NaNs in min() and max()
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.15 (2011-03-04)}{
  \itemize{
    \item Upgraded to Armadillo version 1.1.6   “Baby Carpet Shark”
    \itemize{
      \item fixed size matrices and vectors can use auxiliary (external) memory
      \item .in_range() can use span() arguments
      \item subfields can use span() arguments
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.14 (2011-03-02)}{
  \itemize{
    \item Support Run-Time Type Information (RTTI) on matrices by setting the state variable vec_state in Row and Col instantiation, with thanks to Conrad Sanderson for the hint
    \item fastLm code simplified further by instantiating the Armadillo matrix and vector directly from the SEXP coming from R
    \item inst/doc/Makefile now respects $R_HOME environment variable
  }
}
\section{Changes in RcppArmadillo version 0.2.13 (2011-02-18)}{
  \itemize{
    \item Upgraded to Armadillo version 1.1.4   “Manta Lodge”
    \itemize{
      \item Faster sort()
      \item Updated installation to detect recent versions of Intel's MKL
      \item Added interpretation of arbitrary "flat" subcubes as matrices
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.12 (2011-02-15)}{
  \itemize{
    \item Upgraded to Armadillo version 1.1.2   “Flood Kayak”
    \itemize{
      \item Faster prod()
      \item Faster solve() for compound expressions
      \item Fix for compilation using GCC's C++0x mode
      \item Fix for matrix handling by subcubes
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.11 (2011-01-06)}{
  \itemize{
    \item Upgraded to Armadillo version 1.1.0   “Climate Vandal”
    \itemize{
      \item Extended submatrix views, including access to elements whose indices are specified in a separate vector
      \item Added handling of raw binary files by save/load functions
      \item Added cumsum()
      \item Added interpretation of matrices as triangular via trimatu()/trimatl() 
      \item Faster solve(), inv() via explicit handling of triangular matrices 
      \item The stream for logging of errors and warnings can now be changed 
    }
    \item New unexported R function SHLIB, a small wrapper around R CMD SHLIB, which can be used as   Rscript -e "RcppArmadillo:::SHLIB('foo.cpp')"
  }
}
\section{Changes in RcppArmadillo version 0.2.10 (2010-11-25)}{
  \itemize{
    \item Upgraded to Armadillo 1.0.0 "Antipodean Antileech"
    \itemize{
      \item After 2 1/2 years of collaborative development, we are proud to release the 1.0 milestone version. 
      \item Many thanks are extended to all contributors and bug reporters.
    }
    \item R/RcppArmadillo.package.skeleton.R: Updated to no longer rely on GNU make for builds of packages using RcppArmadillo
    \item summary() for fastLm() objects now returns r.squared and adj.r.squared
  }
}
\section{Changes in RcppArmadillo version 0.2.9 (2010-11-11)}{
  \itemize{
    \item Upgraded to Armadillo 0.9.92 "Wall Street Gangster":
    \itemize{
      \item Fixes for compilation issues under the Intel C++ compiler
      \item Added matrix norms
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.8 (2010-10-16)}{
  \itemize{
    \item Upgraded to Armadillo 0.9.90 "Water Dragon":
    \itemize{
      \item Added unsafe_col()
      \item Speedups and bugfixes in lu()
      \item Minimisation of pedantic compiler warnings
    }
    \item Switched NEWS and ChangeLog between inst/ and the top-level directory so that NEWS (this file) gets installed with the package
  }
}
\section{Changes in RcppArmadillo version 0.2.7 (2010-09-25)}{
  \itemize{
    \item Upgraded to Armadillo 0.9.80 "Chihuahua Muncher":
    \itemize{
      \item Added join_slices(), insert_slices(), shed_slices()
      \item Added in-place operations on diagonals
      \item Various speedups due to internal architecture improvements
    }
  }
}
\section{Changes in RcppArmadillo version 0.2.6 (2010-09-12)}{
  \itemize{
    \item Upgraded to Armadillo 0.9.70 "Subtropical Winter Safari"
    \item arma::Mat, arma::Row and arma::Col get constructor that take vector or matrix sugar expressions. See the unit test "test.armadillo.sugar.ctor" and "test.armadillo.sugar.matrix.ctor" for examples.
  }
}
\section{Changes in RcppArmadillo version 0.2.5 (2010-08-05)}{
  \itemize{
    \item Upgraded to Armadillo 0.9.60 "Killer Bush Turkey"
  }
}
\section{Changes in RcppArmadillo version 0.2.4 (2010-07-27)}{
  \itemize{
    \item Upgraded to Armadillo 0.9.52 'Monkey Wrench'
    \item src/fastLm.cpp: Switch from inv() to pinv() as inv() now tests for singular matrices and warns and returns an empty matrix which stops the example fastLm() implementation on the manual page -- and while this is generally reasonably it makes sense here to continue which the Moore-Penrose pseudo-inverse allows us to do this
  }
}
\section{Changes in RcppArmadillo version 0.2.3 (2010-06-14)}{
  \itemize{
    \item Better configuration to detect suncc (which does not have std::isfinite)
  }
}
\section{Changes in RcppArmadillo version 0.2.2 (2010-06-09)}{
  \itemize{
    \item Added RcppArmadillo:::CxxFlags for cases where RcppArmadillo is not used via a package
    \item Upgraded to Armadillo 0.9.10 'Chilli Espresso'
    \item Wrap support for mtOp, i.e. operations involving mixed types such as a complex and an arma::mat, which have been introduced in armadillo 0.9.10
    \item Wrap support for mtGlue, i.e. operations involving matrices of mixed types such as an arma::mat and an arma::imat, which have been introduced in armadillo 0.9.10
    \item Included an inline plugin to support the plugin system introduced in inline 0.3.5. The unit tests have moved from the src directory to the unit test directory (similar to Rcpp) using cxxfunction with the RcppArmadillo plugin. 
  }
}
\section{Changes in RcppArmadillo version 0.2.1 (2010-05-19)}{
  \itemize{
    \item Bug-fix release permitting compilation on Windows
  }
}
\section{Changes in RcppArmadillo version 0.2.0 (2010-05-18)}{
  \itemize{
    \item fastLm() is now generic and has a formula interface as well as methods for print, summary, predict to behave like a standard model fitting function
    \item Armadillo sources (using release 0.9.8) are now included in the package using a standardized build suitable for our purposes (not assuming Boost or Atlas) -- see ?RcppArmadillo for details
    \item New R function RcppArmadillo.package.skeleton, similar to Rcpp::Rcpp.package.skeleton, but targetting use of RcppArmadillo
  }
}
\section{Changes in RcppArmadillo version 0.1.0 (2010-03-11)}{
  \itemize{
    \item the fastLm() implementation of a bare-bones lm() fit (using Armadillo's solve() function) provides an example of how efficient code can be written compactly using the combination of Rcpp, RcppAramadillo and Armadillo
    \item support for Rcpp implicit wrap of these types : Mat<T>, Col<T>, Row<T>, Cube<T> where T is one of : int, unsigned int, double, float
    \item support for Rcpp implicit as of these types : Mat<T>, Col<T>, Row<T> where R is one of : int, unsigned int, double, float
  }
}
