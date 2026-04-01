ifeq ($(OS), Windows_NT)
	TEXI2DVI = texify --pdf --quiet
else
	TEXI2DVI = texi2dvi --pdf --quiet
endif

libcoin.Rnw: libcoin.w
	mkdir -p inst/C_API_example/{man,R,src} inst/{include,nuweb} man R src
	nuweb -r libcoin.w
	cp src/*.h inst/include/
	cp libcoin.w inst/nuweb/
	mv libcoin.tex libcoin.Rnw

libcoin.R: libcoin.Rnw
	R CMD Stangle libcoin.Rnw
#	mv libcoin.R tests

libcoin.tex: libcoin.Rnw
	R CMD Sweave libcoin.Rnw

libcoin.pdf: libcoin.tex
	pdflatex --draftmode --quiet libcoin.tex
	### NEED to run nuweb AGAIN to sort out page numbers
	### WITHOUT regeneration of the output files
	nuweb -or libcoin.w
	mv libcoin.tex libcoin.Rnw
	R CMD Sweave libcoin.Rnw
	mkdir -p vignettes
	cp libcoin.Rnw vignettes/
	$(TEXI2DVI) libcoin.tex

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
