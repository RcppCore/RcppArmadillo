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
             AC_MSG_RESULT([yes: ${gxx_version}])
	;;
    esac
fi
