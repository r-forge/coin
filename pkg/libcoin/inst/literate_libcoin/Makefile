
libcoin.pdf:	libcoin.w
	nuweb libcoin.w
	texi2dvi --pdf libcoin.tex
	nuweb libcoin.w
	texi2dvi --pdf libcoin.tex
	pdflatex libcoin.tex
	R CMD SHLIB Sums.c

all:	libcoin.pdf

clean:	
	rm -rf libcoin.pdf
	rm -rf *c	