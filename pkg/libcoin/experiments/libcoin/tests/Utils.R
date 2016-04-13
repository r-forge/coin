
library("libcoin")

set.seed(29)
p <- 3
x <- crossprod(matrix(runif(p * p), p))

svd(x)

.Call("R_svd", x)

.Call("R_MPinv", x, sqrt(.Machine$double.eps))

partykit:::.MPinv(x)
