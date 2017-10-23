
libcoin.Rnw:	libcoin.w
	nuweb -r libcoin.w
	mv libcoin.tex libcoin.Rnw
	mv *.h src
	mv *.c src
	mv Makevars src
	mv *def src
	mv *.R R

libcoin.R: libcoin.Rnw
	Stangle libcoin.Rnw
	mv libcoin.R tests

libcoin.tex:	libcoin.Rnw
	R CMD Sweave libcoin.Rnw
	
libcoin.pdf:	libcoin.tex
	texi2dvi --pdf libcoin.tex
	nuweb -r libcoin.w
	mv libcoin.tex libcoin.Rnw
	R CMD Sweave libcoin.Rnw
	pdflatex libcoin.tex

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
	rm -rf src/*
	rm -rf R/*