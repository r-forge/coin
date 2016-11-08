
library("libcoin")

### by Henric Winell
X <- runif(10)
Y <- runif(10)
o <- LinStatExpCov(X, Y)
ov <- LinStatExpCov(X, Y, varonly = TRUE)
stopifnot(all.equal(doTest(o, teststat = "scalar"),
                    doTest(ov, teststat = "scalar")))
