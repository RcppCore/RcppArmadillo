
all: 	RcppArmadillo-unitTests.pdf RcppArmadillo-intro.pdf

pdfclean:
	rm -f RcppArmadillo-unitTests.pdf RcppArmadillo-intro.pdf

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

RcppArmadillo-intro.pdf:
	${RPROG} CMD Sweave RcppArmadillo-intro.Rnw
	pdflatex -shell-escape RcppArmadillo-intro
	bibtex RcppArmadillo-intro
	pdflatex -shell-escape RcppArmadillo-intro
	pdflatex -shell-escape RcppArmadillo-intro
	rm RcppArmadillo-intro.aux RcppArmadillo-intro.log RcppArmadillo-intro.bbl RcppArmadillo-intro.blg RcppArmadillo-intro.tex RcppArmadillo-intro.spl
