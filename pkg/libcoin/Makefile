ifeq ($(OS), Windows_NT)
	TEXI2DVI = texify --pdf --quiet
else
	TEXI2DVI = texi2dvi --pdf --quiet
endif

libcoin.Rnw: libcoin.w
	nuweb -r libcoin.w
	mv libcoin.tex libcoin.Rnw
	mv NEWS.Rd inst/
	mv *.Rd man/
	mv *.R R/
	mv *.c *.def *.h Makevars src/

libcoin.R: libcoin.Rnw
	R CMD Stangle libcoin.Rnw
#	mv libcoin.R tests

libcoin.tex: libcoin.Rnw
	R CMD Sweave libcoin.Rnw

libcoin.pdf: libcoin.tex
	$(TEXI2DVI) libcoin.tex
	### NEED to run nuweb AGAIN to sort out page numbers
	### WITHOUT regeneration of the output files
	nuweb -or libcoin.w
	mv libcoin.tex libcoin.Rnw
	R CMD Sweave libcoin.Rnw
	$(TEXI2DVI) libcoin.tex
	cp libcoin.Rnw vignettes

all: libcoin.pdf

clean:
	./cleanup
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
	rm -rf *bbl
	rm -rf *brf
	rm -rf *wk

distclean: clean
	rm -rf src/*
	rm -rf R/*
	rm -rf man/*
	rm -rf vignettes/*Rnw
	rm -rf tests/libcoin.R*
