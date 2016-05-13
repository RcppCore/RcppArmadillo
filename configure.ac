## -*- mode: autoconf; autoconf-indentation: 4; -*-
##
##  RcppArmadillo configure.ac
##
##  'Rcpp' Integration for the 'Armadillo' Templated Linear Algebra Library
##
##  Copyright (C) 2016  Dirk Eddelbuettel 
##
##  Licensed under GPL-2 or later

## require at least autoconf 2.61
AC_PREREQ(2.61)

## Process this file with autoconf to produce a configure script.
AC_INIT(RcppArmadillo, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

## Check the C++ compiler
AC_PROG_CXX
## If it is g++, we have GXX set so let's examine it
if test "${GXX}" = yes; then
    AC_MSG_CHECKING([whether g++ version is sufficient])
    gxx_version=$(${CXX} -v 2>&1 | awk '/^.*g.. version/ {print $3}')
    case ${gxx_version} in
        1.*|2.*|3.*|4.0.*|4.1.*|4.2.*|4.3.*|4.4.*|4.5.*)
             AC_MSG_RESULT([no])
	     AC_MSG_WARN([Only g++ version 4.6 or greater can be used with RcppArmadillo.])
	     AC_MSG_ERROR([Please use a different compiler.])   
        ;;
	4.6.*|4.7.*|4.8.*|4.9.*|5.*|6.*)
	     gxx_newer_than_45="-fpermissive"
             AC_MSG_RESULT([yes (${gxx_version})])
	;;
    esac
fi


## Next, set R_HOME, respecting an environment variable if set 
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi

## Check for suitable LAPACK_LIBS
AC_MSG_CHECKING([LAPACK_LIBS])

## external LAPACK has the required function
lapack=$(${R_HOME}/bin/R CMD config LAPACK_LIBS)
hasRlapack=$(echo ${lapack} | grep lRlapack)
## internal Rlapack now the required functions if "new enough": R 3.3.0 or later
newR=$("${R_HOME}/bin/Rscript" --vanilla -e 'if (getRversion() >= "3.3.0") cat("yes") else cat("no")')

## in what follows below we substitute both side of the define/undef
## while this may seem a little unusual we do it to fully reproduce the
## previous bash-based implementation

if test x"${hasRlapack}" == x""; then
    ## We are using a full Lapack and can use zgesdd -- so #undef remains
    AC_MSG_RESULT([system LAPACK found])
    have_good_lapack="#define"
    arma_crippled_lapack="ARMA_CRIPPLED_LAPACK 1"
else
    if test x"${newR}" == x"yes"; then
        ## The R versions are recent enough and has an augmented internal Rlapack
        AC_MSG_RESULT([fallback LAPACK from R 3.3.0 or later used])
        have_good_lapack="#define"
        arma_crippled_lapack="ARMA_CRIPPLED_LAPACK 1"
    else     
        ## We are using an older R's subset of Lapack and CANNOT use zgesdd etc, so we mark it
        AC_MSG_RESULT([R-supplied partial LAPACK found])
        AC_MSG_WARN([Some complex-valued LAPACK functions may not be available])
        have_good_lapack="#undef"
        arma_crippled_lapack="ARMA_CRIPPLED_LAPACK"
    fi
fi

## now use all these
AC_SUBST([HAVE_GOOD_LAPACK],["${have_good_lapack}"])
AC_SUBST([ARMA_CRIPPLED_LAPACK], ["${arma_crippled_lapack}"])
AC_CONFIG_FILES([inst/include/RcppArmadilloLapack.h])
AC_OUTPUT
