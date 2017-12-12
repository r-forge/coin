
libcoin.Rnw:	libcoin.w
	nuweb -r libcoin.w
	mv libcoin.tex libcoin.Rnw
	cp libcoin.Rnw vignettes
	mv *.h src
	mv *.c src
	mv Makevars src
	mv *def src
	mv *.R R
	mv *.Rd man

libcoin.R: libcoin.Rnw
	Stangle libcoin.Rnw
	mv libcoin.R tests

libcoin.tex:	libcoin.Rnw
	R CMD Sweave libcoin.Rnw
	
libcoin.pdf:	libcoin.tex
	texi2dvi --pdf libcoin.tex
	### NEED to run nuweb TWICE to
	### sort-out page numbers
	nuweb -r libcoin.w
	mv libcoin.tex libcoin.Rnw
	R CMD Sweave libcoin.Rnw
	pdflatex libcoin.tex
	rm -rf libcoin.Rnw

all:	libcoin.R

clean:	
	rm -rf *~
	rm -rf libcoin.pdf
	rm -rf libcoin.R
	texclean libcoin
	rm -rf libcoin.tex
	rm -rf libcoin.Rnw
	rm -rf *.c
	rm -rf *.h
	rm -rf *.o
	rm -rf *.so	
	rm -rf *.def
	rm -rf *.Rd
	rm -rf *.R
	rm -rf Makevars
	rm -rf nw*
	rm -rf src/*
	rm -rf R/*
	rm -rf man/*
	rm -rf vignettes/*Rnw
	rm -rf tests/libcoin.R*
	rm -rf *bbl
	rm -rf *brf
	rm -rf *wk
		
	