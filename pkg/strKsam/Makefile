
strKsam.Rnw:	strKsam.w
	nuweb -r strKsam.w
	mv strKsam.tex strKsam.Rnw
	#mv *.h src
	mv Schur_src.c src/Schur.c
	#mv Makevars src
	#mv *def src
	mv free1way.R R
	mv linkfun.R R
	#mv *.Rd man

strKsam.R: strKsam.Rnw
	Stangle strKsam.Rnw
#	mv strKsam.R tests

strKsam.tex:	strKsam.Rnw
	R CMD Sweave strKsam.Rnw
	
strKsam.pdf:	strKsam.tex
	texi2dvi --pdf strKsam.tex
	### NEED to run nuweb TWICE to
	### sort-out page numbers
	nuweb -r strKsam.w
	mv strKsam.tex strKsam.Rnw
	R CMD Sweave strKsam.Rnw
	pdflatex strKsam.tex
	pdflatex strKsam.tex
	# cp strKsam.Rnw vignettes

all:	strKsam.pdf

clean:	
	./cleanup
	rm -rf *~
	rm -rf strKsam.pdf
	rm -rf strKsam.R
	texclean strKsam
	rm -rf strKsam.tex
	rm -rf strKsam.Rnw
	rm -rf *.c
	rm -rf *.h
	rm -rf *.o
	rm -rf *.so	
	rm -rf *.def
	rm -rf *.Rd
	rm -rf *.R
	rm -rf Makevars
	rm -rf nw*
	rm -rf *bbl
	rm -rf *brf
	rm -rf *wk

distclean: clean
	rm -rf src/*
	rm -rf R/*
	rm -rf man/*
	rm -rf vignettes/*Rnw
	rm -rf tests/strKsam.R*
	