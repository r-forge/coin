
free1way.Rnw:	free1way.w
	nuweb -r free1way.w
	mv free1way.tex free1way.Rnw
	cp free1way.Rnw vignettes
	# mv Schur.c src
	#mv Makevars src
	#mv *def src
	mv free1way.R R
	mv linkfun.R R
	#mv *.Rd man

free1way.R: free1way.Rnw
	Stangle free1way.Rnw
#	mv free1way.R tests

free1way.tex:	free1way.Rnw
	R CMD Sweave free1way.Rnw
	
free1way.pdf:	free1way.tex
	texi2dvi --pdf free1way.tex
	### NEED to run nuweb TWICE to
	### sort-out page numbers
	nuweb -r free1way.w
	mv free1way.tex free1way.Rnw
	R CMD Sweave free1way.Rnw
	pdflatex free1way.tex
	pdflatex free1way.tex
	# cp free1way.Rnw vignettes

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
	