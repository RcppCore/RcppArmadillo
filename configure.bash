#!/usr/bin/env bash
#
# Copyright 2013 - 2016  Dirk Eddelbuettel
#
# Licensed under GPL 2 or later

# This file uses bash.
#
# If this file is not suitable for your system for lack of bash or
# another suitable /bin/sh implementation, I recommend deletion of the file
# along with *manual* adjustment to RcppArmadilloLapack.h as per the
# test below. 
#
# In most case you can just copy RcppArmadilloLapack.h.in over to
# RcppArmadilloLapack.h -- and in aggregate it is not worth rewriting this.
# At some point in the future we will be able to just rely on recent enough R
# with good enough embeded Lapack, or use of external Lapack.

if [ "${R_HOME}" == "" ]; then
   R_HOME=$(R RHOME)
fi

echo -n "* checking LAPACK_LIBS: "

## external LAPACK has the required function
lapack=$(${R_HOME}/bin/R CMD config LAPACK_LIBS)
hasRlapack=$(echo ${lapack} | grep lRlapack)

## internal Rlapack now the required functions if "new enough": R 3.3.0 or later
newR=$(${R_HOME}/bin/R --slave -q -e 'cat(ifelse(getRversion() >= "3.3.0","yes","no"))')

if [ "${hasRlapack}" == "" ]; then
    ## We are using a full Lapack and can use zgesdd -- so #undef remains
    echo "system LAPACK found"
    cp inst/include/RcppArmadilloLapack.h.in inst/include/RcppArmadilloLapack.h 
elif [ "$newR" == "yes" ]; then
    ## The R versions are recent enough and has an augmented internal Rlapack
    echo "fallback LAPACK from R 3.3.0 or later used"
    cp inst/include/RcppArmadilloLapack.h.in inst/include/RcppArmadilloLapack.h 
else
    ## We are using a R's subset of Lapack and CANNOT use zgesdd etc, so we mark it
    echo "R-supplied partial LAPACK found"
    echo "* some operations may not be available"
    sed -e 's/\/\/ \#undef ARMA_CRIPPLED_LAPACK/\#define ARMA_CRIPPLED_LAPACK 1/' \
        inst/include/RcppArmadilloLapack.h.in > inst/include/RcppArmadilloLapack.h 
fi

function checkoldcompiler() {
    ## store compiler argument (obtained via, say, R CMD config CXX)
    cxx=${1}
    ## try the compiler, and see if we can extract a version following g..
    ## not that this only works for 
    match=$(${cxx} -v 2>&1 | awk '/^.*g.. version/ {print $3}')
    ## now test the match (if we have one)
    if [ x${match} != x ]; then
        case ${match} in
            1.*|2.*|3.*|4.0.*|4.1.*|4.2.*|4.3.*|4.4.*|4.5.*)
	        echo "* ERROR: Your g++ version ${match} is too old to build RcppArmadillo. "
                echo "* ERROR: You must upgrade to a version supporting a minimum of C++11."
                echo "* ERROR: At least g++ 4.6.* or later is recommended."
                exit -1
                ;;
	    4.6.*|4.7.*|4.8.*|4.9.*|5.*|6.*)
	        #echo "GOOD: g++ version ${match} is sufficient."
	        ;;
        esac
    fi
}

checkoldcompiler "$(R CMD config CXX)"
checkoldcompiler "$(R CMD config CXX1X)"
exit 0
