
library("libcoin")
library("coin")
set.seed(29)

n <- 100
p <- 4
q <- 2
X <- matrix(runif(p * n), nc = p)
Y <- matrix(runif(q * n), nc = q)
w <- as.integer(floor(runif(n, max = 4)))
s <- sample(1:n, floor(n/2), replace = TRUE)
b <- sample(gl(2, 2, length = n))

cmp <- function(t1, t2) {
    if (is.null(t1$Covariance)) {
        var1 <- t1$Variance
        var2 <- diag(covariance(t2))
    } else {
        var1 <- t1$Covariance
        var2 <- covariance(t2)
    }
    c(max(abs(t1$LinStat - statistic(t2, "linear"))),
      max(abs(t1$Expectation - expectation(t2))),
      max(abs(var1 - var2)))
}

t1 <-LinStatExpCov(X, Y)
t1v <-LinStatExpCov(X, Y, varonly = TRUE)
t2 <- independence_test(Y ~ X)
cmp(t1, t2)
cmp(t1v, t2)

t1 <-LinStatExpCov(X, Y, weights = w)
t1v <-LinStatExpCov(X, Y, weights = w, varonly = TRUE)
t2 <- independence_test(Y ~ X, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, subset = s)
t1v <- LinStatExpCov(X, Y, subset = s, varonly = TRUE)
t2 <- independence_test(Y ~ X, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, weights = w, subset = s)
t1v <- LinStatExpCov(X, Y, weights = w, subset = s, varonly = TRUE)
t2 <- independence_test(Y ~ X, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, block = b)
t1v <- LinStatExpCov(X, Y, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X  | b)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, weights = w, block = b)
t1v <- LinStatExpCov(X, Y, weights = w, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov(X, Y, weights = w, subset = s, block = b)
t1v <- LinStatExpCov(X, Y, weights = w, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y ~ X | b, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

n <- 100
n1 <- 5
n2 <- 4
p <- 4
q <- 2
X <- rbind(0, matrix(runif(p * n1), nc = p))
Y <- rbind(0, matrix(runif(q * n2), nc = q))
ix <- sample(1:n1, n, replace = TRUE)
iy <- sample(1:n2, n, replace = TRUE)
w <- as.integer(floor(runif(n, max = 4)))
s <- sample(1:n, floor(n/2), replace = TRUE)
b <- sample(gl(2, 2, length = n))

t1 <- LinStatExpCov2d(X, Y, ix, iy)
t1v <- LinStatExpCov2d(X, Y, ix, iy, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,])
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov2d(X, Y, ix, iy, weights = w)
t1v <- LinStatExpCov2d(X, Y, ix, iy, weights = w, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,], weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov2d(X, Y, ix, iy, subset = s)
t1v <- LinStatExpCov2d(X, Y, ix, iy, subset = s, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,], subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov2d(X, Y, ix, iy, weights = w, subset = s)
t1v <- LinStatExpCov2d(X, Y, ix, iy, weights = w, subset = s, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,], weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov2d(X, Y, ix, iy, block = b)
t1v <- LinStatExpCov2d(X, Y, ix, iy, block = b, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,]  | b)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov2d(X, Y, ix, iy, weights = w, block = b)
t1v <- LinStatExpCov2d(X, Y, ix, iy, weights = w, block = b, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,] | b, weights = ~ w)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov2d(X, Y, ix, iy, subset = s, block = b)
t1v <- LinStatExpCov2d(X, Y, ix, iy, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,] | b, subset = s)
cmp(t1, t2)
cmp(t1v, t2)

t1 <- LinStatExpCov2d(X, Y, ix, iy, weights = w, subset = s, block = b)
t1v <- LinStatExpCov2d(X, Y, ix, iy, weights = w, subset = s, block = b, varonly = TRUE)
t2 <- independence_test(Y[iy + 1,] ~ X[ix + 1,]| b, weights = ~w, subset = s)
cmp(t1, t2)
cmp(t1v, t2)
