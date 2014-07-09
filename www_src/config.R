
pkg <- "coin"
dest <- "html"
publish <- "../www"

download.file("http://user.math.uzh.ch/hothorn/TH.bib", dest = "TH.bib")
system("cat coin.bib >> TH.bib")
bib <- "TH.bib"