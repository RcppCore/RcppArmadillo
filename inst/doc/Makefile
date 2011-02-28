
all: 	clean RcppArmadillo-unitTests.pdf

clean:
	rm -f RcppArmadillo-unitTests.pdf

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif
RPROG=	$(R_HOME)/bin/R
RSCRIPT=$(R_HOME)/bin/Rscript

RcppArmadillo-unitTests.pdf:
	rm -fr unitTests-results/*
	$(RSCRIPT) unitTests/RcppArmadillo-unitTests.R
	$(RPROG) CMD Sweave RcppArmadillo-unitTests.Rnw
	$(RSCRIPT) -e "tools::texi2dvi( 'RcppArmadillo-unitTests.tex', pdf = TRUE, clean = TRUE )"
	rm -fr RcppArmadillo-unitTests.tex

