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
AC_INIT(RcppArmadillo, 0.6.700.6.0)

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

## Check the C++ compiler
AC_PROG_CXX
## If it is g++, we have GXX set so let's examine it
if test "${GXX}" = yes; then
    AC_MSG_CHECKING([whether g++ version is sufficient])
    gxx_version=`${CXX} -v 2>&1 | awk '/^.*g.. version/ {print $3}'`
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

dnl ## external LAPACK has the required function
dnl lapack=$(${R_HOME}/bin/R CMD config LAPACK_LIBS)
dnl hasRlapack=$(echo ${lapack} | grep lRlapack)

dnl ## internal Rlapack now the required functions if "new enough": R 3.3.0 or later
dnl newR=$(${R_HOME}/bin/R --slave -q -e 'cat(ifelse(getRversion() >= "3.3.0","yes","no"))')

dnl if [ "${hasRlapack}" == "" ]; then
dnl     ## We are using a full Lapack and can use zgesdd -- so #undef remains
dnl     echo "system LAPACK found"
dnl     cp inst/include/RcppArmadilloLapack.h.in inst/include/RcppArmadilloLapack.h 
dnl elif [ "$newR" == "yes" ]; then
dnl     ## The R versions are recent enough and has an augmented internal Rlapack
dnl     echo "fallback LAPACK from R 3.3.0 or later used"
dnl     cp inst/include/RcppArmadilloLapack.h.in inst/include/RcppArmadilloLapack.h 
dnl else
dnl     ## We are using a R's subset of Lapack and CANNOT use zgesdd etc, so we mark it
dnl     echo "R-supplied partial LAPACK found"
dnl     echo "* some operations may not be available"
dnl     sed -e 's/\/\/ \#undef ARMA_CRIPPLED_LAPACK/\#define ARMA_CRIPPLED_LAPACK 1/' \
dnl         inst/include/RcppArmadilloLapack.h.in > inst/include/RcppArmadilloLapack.h 
dnl fi

## BAD
have_good_lapack="#undef"
arma_crippled_lapack="ARMA_CRIPPLED_LAPACK"

## GOOD
have_good_lapack="#define"
arma_crippled_lapack="ARMA_CRIPPLED_LAPACK 1"

## now use all these
AC_SUBST([HAVE_GOOD_LAPACK],["${have_good_lapack}"])
AC_SUBST([ARMA_CRIPPLED_LAPACK], ["${arma_crippled_lapack}"])
AC_CONFIG_FILES([inst/include/RcppArmadilloLapack.h])
AC_OUTPUT
