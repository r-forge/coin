
library("lehmann")

set.seed(29)
N <- 20
P <- 10
a <- rnorm(N)
b <- rnorm(N - 1)
#b[floor(1:3 / 4 * N)] <- 0
b[floor(N / 2)] <- 0
X <- matrix(rnorm(N * P), nrow = N)
Z <- matrix(rnorm(P * P), nrow = P)

A <- diag(a / 2)
A[cbind(1:(N-1),2:N)] <- b
A <- t(A) + A

Z1a <- Z - t(X) %*% solve(A) %*% X
Z1b <- lehmann:::Schur_symtri(a = a, b = b, X = X, Z = Z)

stopifnot(isTRUE(all.equal(Z1a, Z1b)))
