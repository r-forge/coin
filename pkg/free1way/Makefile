
free1way.Rnw:	free1way.w
	nuweb -r free1way.w
	mv free1way.tex free1way.Rnw
	### this version doesn't contain proper page numbers
	### but "make" will allow a quick R CMD check run
	cp free1way.Rnw vignettes
	mv free1way.R R
	mv linkfun.R R

free1way.R: free1way.Rnw
	Stangle free1way.Rnw

free1way.tex:	free1way.Rnw
	Rdev CMD Sweave free1way.Rnw
	
free1way.pdf:	free1way.tex
	texi2dvi --pdf free1way.tex
	### NEED to run nuweb TWICE to
	### sort-out page numbers
	nuweb -r free1way.w
	mv free1way.tex free1way.Rnw
	### copy again with correct page numbers
	cp free1way.Rnw vignettes
	Rdev CMD Sweave free1way.Rnw
	pdflatex free1way.tex
	pdflatex free1way.tex

all:	free1way.pdf

clean:	
	./cleanup
	rm -rf *~
	rm -rf *pdf
	rm -rf free1way.pdf
	rm -rf free1way.R
	texclean free1way
	rm -rf free1way.tex
	rm -rf free1way.Rnw
	rm -rf *.c
	rm -rf *.h
	rm -rf *.o
	rm -rf *.so	
	rm -rf *.def
	rm -rf *.R
	rm -rf Makevars
	rm -rf nw*
	rm -rf *bbl
	rm -rf *aux
	rm -rf *log
	rm -rf *blg
	rm -rf *out
	rm -rf *toc
	rm -rf *brf
	rm -rf *wk

distclean: clean
	rm -rf src/Schur.c
	rm -rf R/*
	rm -rf vignettes/*Rnw
	rm -rf tests/free1way.R*
	