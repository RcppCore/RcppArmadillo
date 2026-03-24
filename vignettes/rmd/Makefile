
## This Makefile is not included in the package sources so we can use GNUmake idioms. Yay.

rmdsources := 		$(wildcard *.Rmd)
rmdvignettes := 	$(rmdsources:.Rmd=.pdf)

%.pdf:			%.Rmd
			Rscript -e 'rmarkdown::render("$<")'
			Rscript -e 'tools::compactPDF("$@", gs_quality="ebook")'
			cp -vax $@ ..

all:			${rmdvignettes} RcppArmadillo-intro.pdf

RcppArmadillo-intro.pdf: RcppArmadillo-intro.Rnw
			Rscript -e 'Sweave("$<")'
			pdflatex RcppArmadillo-intro.tex
			bibtex RcppArmadillo-intro.aux
			pdflatex RcppArmadillo-intro.tex
			pdflatex RcppArmadillo-intro.tex
			Rscript -e 'tools::compactPDF("$@", gs_quality="ebook")'
			cp -vax $@ ..

clean:
			@rm -rf *.aux *.log *.out *.tex *.spl *.bbl *.blg Rcpp*.pdf \
				RcppArmadillo*_cache pinp.cls jss.bst auto/
