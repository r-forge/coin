
libcoin.pdf:	libcoin.w
	nuweb libcoin.w
	mv libcoin.tex libcoin.Rnw
	Sweave libcoin.Rnw
	texi2dvi --pdf libcoin.tex
	nuweb libcoin.w
	mv libcoin.tex libcoin.Rnw
	Sweave libcoin.Rnw
	texi2dvi --pdf libcoin.tex
	pdflatex libcoin.tex
	R CMD SHLIB Sums.c

libcoin.R:	libcoin.pdf
	Stangle libcoin.Rnw

all:	libcoin.pdf

clean:	
	rm -rf libcoin.pdf
	rm -rf *c	