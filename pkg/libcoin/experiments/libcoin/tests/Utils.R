
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


set.seed(29)
n <- 50
X <- crossprod(matrix(runif(n^2), nr = n))
x <- X[!upper.tri(X)]
system.time(g1 <- partykit:::.MPinv(X))
system.time(g2 <- .Call("R_MPinv_sym", x, sqrt(.Machine$double.eps)))
names(g2) <- c("Xplus", "rank")
all.equal(g1, g2)

