
libcoin.pdf:	libcoin.w
	nuweb -r libcoin.w
	R CMD SHLIB Sums.c
	mv libcoin.tex libcoin.Rnw
	Sweave libcoin.Rnw
	texi2dvi --pdf libcoin.tex
	nuweb -r libcoin.w
	mv libcoin.tex libcoin.Rnw
	Sweave libcoin.Rnw
	texi2dvi --pdf libcoin.tex
	pdflatex libcoin.tex

libcoin.R:	libcoin.pdf
	Stangle libcoin.Rnw

all:	libcoin.pdf libcoin.R

clean:	
	rm -rf libcoin.pdf
	rm -rf libcoin.R
	texclean libcoin
	rm -rf *c
	rm -rf *h	