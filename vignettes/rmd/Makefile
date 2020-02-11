
## This Makefile is not included in the package sources so we can use GNUmake idioms. Yay.

rmdsources := 		$(wildcard *.Rmd)
rmdvignettes := 	$(rmdsources:.Rmd=.pdf)

%.pdf:			%.Rmd
			Rscript -e 'rmarkdown::render("$<")'
			Rscript -e 'tools::compactPDF("$@", gs_quality="ebook")'
			cp -vax $@ ../pdf

all:			${rmdvignettes} RcppArmadillo-intro.pdf

RcppArmadillo-intro.tex: RcppArmadillo-intro.Rnw
			Rscript -e 'Sweave("$<")'

RcppArmadillo-intro.pdf: RcppArmadillo-intro.tex
			Rscript -e 'tools::texi2pdf("$<", texi2dvi="pdflatex")'
			Rscript -e 'tools::texi2pdf("$<", texi2dvi="pdflatex")'
			Rscript -e 'tools::compactPDF("$@", gs_quality="ebook")'
			cp -vax $@ ../pdf

clean:
			@rm -rf *.aux *.log *.out *.tex *.spl Rcpp*.pdf \
				RcppArmadillo*_cache pinp.cls jss.bst
