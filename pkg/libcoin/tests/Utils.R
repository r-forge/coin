
library("libcoin")
set.seed(29)

n <- 50
X <- crossprod(matrix(runif(n^2), nr = n))
x <- X[!upper.tri(X)]
system.time(g1 <- partykit:::.MPinv(X))
system.time(g2 <- .Call("R_MPinv_sym", x, sqrt(.Machine$double.eps)))
g1$Xplus <- g1$Xplus[!upper.tri(g1$Xplus)]
names(g2) <- c("Xplus", "rank")
all.equal(g1, g2)
