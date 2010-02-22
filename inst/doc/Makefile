
all: clean RcppArmadillo-unitTests.pdf

clean:
	touch RcppArmadillo-unitTests.pdf
	rm RcppArmadillo-unitTests.pdf
	
RcppArmadillo-unitTests.pdf:
	rm -fr unitTests-results/*
	Rscript unitTests/RcppArmadillo-unitTests.R
	R CMD Sweave RcppArmadillo-unitTests.Rnw
	Rscript -e "tools::texi2dvi( 'RcppArmadillo-unitTests.tex', pdf = TRUE, clean = TRUE )"
	rm -fr RcppArmadillo-unitTests.tex

