
maxstat.pdf:	maxstat.Rnw
	Rnw2pdf maxstat.Rnw
	
maxstat.tex:	maxstat.Rnw
	R CMD Sweave maxstat.Rnw
	
maxstat.zip:	maxstat.tex
		zip maxstat.zip abstract.tex  acknowledgement.tex  definitions.tex  headerBiometrics.tex  maxstat.tex maxstat-approx.pdf maxstat-maxstat.pdf biometrics.sty biometrics.bst maxstat.bib

all:	maxstat.pdf

zip:	maxstat.zip
	