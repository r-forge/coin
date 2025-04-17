
library("lehmann")

set.seed(29)
N <- 10
P <- 3
a <- runif(N)
b <- runif(N - 1)
X <- matrix(runif(N * P), nrow = N)
Z <- matrix(runif(P * P), nrow = P)

A <- diag(a / 2)
A[cbind(1:(N-1),2:N)] <- b
A <- t(A) + A
A <- cbind(A, X)
A <- rbind(A, cbind(t(X), Z))

Z1a <- solve(A)[-(1:N),-(1:N),drop = FALSE]
Z1b <- Schur_symtri(a = a, b = b, X = X, Z = Z)

stopifnot(isTRUE(all.equal(Z1a, Z1b)))
