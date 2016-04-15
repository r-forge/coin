
library("libcoin")

set.seed(29)
p <- 3
x <- crossprod(matrix(runif(p * p), p))

svd(x)

.Call("R_svd", x)

.Call("R_MPinv", x, sqrt(.Machine$double.eps))

partykit:::.MPinv(x)

set.seed(29)
n <- 100
k <- 20
subset <- sample(1:n, k, replace = FALSE)
block <- gl(2, 2, length = n)

subset[block[subset] == 1] - 1
subset[block[subset] == 2] - 1

.Call("R_PermuteBlock_subset", subset - 1L, block)
gl(3, 5)
.Call("R_PermuteBlock", gl(3, 5))
