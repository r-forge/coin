
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


n <- 1000L
x <- .Call("R_Permute", n)
stopifnot(all.equal(x[[1]], 1:n - 1L))
stopifnot(all.equal(table(x[[1]]), table(x[[2]])))

i <- sample(1:n, floor(n/2)) - 1L
x <- .Call("R_Permute_subset", i)
stopifnot(all.equal(x[[1]], i))
stopifnot(all.equal(table(x[[1]]), table(x[[2]])))

w <- sample(0:3, n, replace = TRUE)
x <- .Call("R_Permute_weights", w)
stopifnot(all.equal(x[[1]], rep(1:n, w) - 1L))
stopifnot(all.equal(table(x[[1]]), table(x[[2]])))

w <- sample(0:3, n, replace = TRUE)
x <- .Call("R_Permute_weights_subset", w, i)
stopifnot(all.equal(x[[1]], rep(i + 1, w[i + 1]) - 1L))
stopifnot(all.equal(table(x[[1]]), table(x[[2]])))

b <- sample(gl(5, n / 5))

x <- .Call("R_PermuteBlock", b)
stopifnot(all.equal(b[x[[1]] + 1], b[x[[2]] + 1]))

i <- sample(1:n, floor(n/2)) - 1L
x <- .Call("R_PermuteBlock_subset", i, b)
stopifnot(all.equal(b[x[[1]] + 1], b[x[[2]] + 1]))

w <- sample(0:3, n, replace = TRUE)
x <- .Call("R_PermuteBlock_weights", w, b)
stopifnot(all.equal(b[x[[1]] + 1], b[x[[2]] + 1]))

w <- sample(0:3, n, replace = TRUE)
x <- .Call("R_PermuteBlock_weights_subset", w, i, b)
stopifnot(all.equal(b[x[[1]] + 1], b[x[[2]] + 1]))


