### R code from vignette source 'libcoin.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: 1dex-1
###################################################
library("libcoin")
set.seed(290875)
x <- gl(5, 20)
y <- round(runif(length(x)), 1)
ls1 <- LinStatExpCov(X = model.matrix(~ x - 1), Y = matrix(y, ncol = 1))
ls1$LinearStatistic
tapply(y, x, sum)


###################################################
### code chunk number 2: 1dex-2
###################################################
ls2 <- LinStatExpCov(X = x, Y = matrix(y, ncol = 1))
all.equal(ls1, ls2)


###################################################
### code chunk number 3: 2dex-1
###################################################
X <- rbind(0, diag(nlevels(x)))
ix <- unclass(x)
ylev <- sort(unique(y))
Y <- rbind(0, matrix(ylev, ncol = 1))
iy <- .bincode(y, breaks = c(-Inf, ylev, Inf))
ls3 <- LinStatExpCov(X = X, ix = ix, Y = Y, iy = iy)
all.equal(ls1, ls3)


###################################################
### code chunk number 4: 2dex-2
###################################################
ls4 <- LinStatExpCov(ix = ix, Y = Y, iy = iy)
all.equal(ls3, ls4)


###################################################
### code chunk number 5: 2dex-3
###################################################
ls3$Table
xtabs(~ ix + iy)


###################################################
### code chunk number 6: vcov-1
###################################################
ls1$Covariance
vcov(ls1)


###################################################
### code chunk number 7: doTest-1
###################################################
### c_max test statistic
### no p-value
doTest(ls1, teststat = "maximum", pvalue = FALSE)
### p-value
doTest(ls1, teststat = "maximum")
### log(p)-value
doTest(ls1, teststat = "maximum", log = TRUE)
### (1-p)-value
doTest(ls1, teststat = "maximum", lower = TRUE)
### log(1 - p)-value
doTest(ls1, teststat = "maximum", lower = TRUE, log = TRUE)
### quadratic
doTest(ls1, teststat = "quadratic")


###################################################
### code chunk number 8: ctabsex-1
###################################################
t1 <- ctabs(ix = ix, iy = iy)
t2 <- xtabs(~ ix + iy)
max(abs(t1[-1, -1] - t2))


###################################################
### code chunk number 9: ex-setup
###################################################
N <- 20L
P <- 3L
Lx <- 10L
Ly <- 5L
Q <- 4L
B <- 2L
iX2d <- rbind(0, matrix(runif(Lx * P), nrow = Lx))
ix <- sample(1:Lx, size = N, replace = TRUE)
levels(ix) <- 1:Lx
ixf <- factor(ix, levels = 1:Lx, labels = 1:Lx)
x <- iX2d[ix + 1,]
Xfactor <- diag(Lx)[ix,]
iY2d <- rbind(0, matrix(runif(Ly * Q), nrow = Ly))
iy <- sample(1:Ly, size = N, replace = TRUE)
levels(iy) <- 1:Ly
iyf <- factor(iy, levels = 1:Ly, labels = 1:Ly)
y <- iY2d[iy + 1,]
weights <- sample(0:5, size = N, replace = TRUE)
block <- sample(gl(B, ceiling(N / B))[1:N])
subset <- sort(sample(1:N, floor(N * 1.5), replace = TRUE))
subsety <- sample(1:N, floor(N * 1.5), replace = TRUE)
r1 <- rep(1:ncol(x), ncol(y))
r1Xfactor <- rep(1:ncol(Xfactor), ncol(y))
r2 <- rep(1:ncol(y), each = ncol(x))
r2Xfactor <- rep(1:ncol(y), each = ncol(Xfactor))


###################################################
### code chunk number 10: Rlibcoin
###################################################
LECV <- function(X, Y, weights = integer(0), subset = integer(0), block = integer(0)) {

    if (length(weights) == 0) weights <- rep(1, nrow(X))
    if (length(subset) == 0) subset <- 1:nrow(X)
    idx <- rep(subset, weights[subset])
    X <- X[idx,,drop = FALSE]
    Y <- Y[idx,,drop = FALSE]
    sumweights <- length(idx)

    if (length(block) == 0) {
        ExpX <- colSums(X)
        ExpY <- colSums(Y) / sumweights
        yc <- t(t(Y) - ExpY)
        CovY <- crossprod(yc) / sumweights
        CovX <- crossprod(X)
        Exp <- kronecker(ExpY, ExpX) 
        Cov <- sumweights / (sumweights - 1) * kronecker(CovY, CovX) - 
               1 / (sumweights - 1) * kronecker(CovY, tcrossprod(ExpX))
 
        ret <- list(LinearStatistic = as.vector(crossprod(X, Y)),
                    Expectation = as.vector(Exp), 
                    Covariance = Cov,
                    Variance = diag(Cov))
   } else {
        block <- block[idx]
        ret <- list(LinearStatistic = 0, Expectation = 0, Covariance = 0, Variance = 0)
        for (b in levels(block)) {
            tmp <- LECV(X = X, Y = Y, subset = which(block == b))
            for (l in names(ret)) ret[[l]] <- ret[[l]] + tmp[[l]]
        }
   }
   return(ret)
}


###################################################
### code chunk number 11: cmpr
###################################################
cmpr <- function(ret1, ret2) {
    if (inherits(ret1, "LinStatExpCov")) {
        if (!ret1$varonly)
            ret1$Covariance <- vcov(ret1)
    }
    ret1 <- ret1[!sapply(ret1, is.null)]
    ret2 <- ret2[!sapply(ret2, is.null)]
    nm1 <- names(ret1)
    nm2 <- names(ret2)
    nm <- c(nm1, nm2)
    nm <- names(table(nm))[table(nm) == 2]
    all.equal(ret1[nm], ret2[nm])
}


###################################################
### code chunk number 12: benchmarks
###################################################
LECVxyws <- LinStatExpCov(x, y, weights = weights, subset = subset)
LEVxyws <- LinStatExpCov(x, y, weights = weights, subset = subset, varonly = TRUE)


###################################################
### code chunk number 13: tests
###################################################
### with X given
testit <- function(...) {
    a <- LinStatExpCov(x, y, ...)
    b <- LECV(x, y, ...)
    d <- LinStatExpCov(X = iX2d, ix = ix, Y = iY2d, iy = iy, ...)
    return(cmpr(a, b) && cmpr(d, b))
}
stopifnot(
    testit() && testit(weights = weights) &&
    testit(subset = subset) && testit(weights = weights, subset = subset) &&
    testit(block = block) && testit(weights = weights, block = block) &&
    testit(subset = subset, block = block) && 
    testit(weights = weights, subset = subset, block = block))
### without dummy matrix X
testit <- function(...) {
    a <- LinStatExpCov(X = ix, y, ...)
    b <- LECV(Xfactor, y, ...)
    d <- LinStatExpCov(X = integer(0), ix = ix, Y = iY2d, iy = iy, ...)
    return(cmpr(a, b) && cmpr(d, b))
}
stopifnot(
    testit() && testit(weights = weights) &&
    testit(subset = subset) && testit(weights = weights, subset = subset) &&
    testit(block = block) && testit(weights = weights, block = block) &&
    testit(subset = subset, block = block) && 
    testit(weights = weights, subset = subset, block = block))


###################################################
### code chunk number 14: permutations-2d
###################################################
LinStatExpCov(X = iX2d, ix = ix, Y = iY2d, iy = iy, weights = weights, subset = subset, nperm = 10)$PermutedLinearStatistic
LinStatExpCov(X = iX2d, ix = ix, Y = iY2d, iy = iy, weights = weights, subset = subset, nperm = 10)$PermutedLinearStatistic


###################################################
### code chunk number 15: ExpectationInfluence
###################################################
sumweights <- sum(weights[subset])
expecty <- a0 <- colSums(y[subset, ] * weights[subset]) / sumweights
a1 <- libcoin:::.libcoinCall("R_ExpectationInfluence", y, weights, subset);
a2 <- libcoin:::.libcoinCall("R_ExpectationInfluence", y, as.double(weights), as.double(subset));
a3 <- libcoin:::.libcoinCall("R_ExpectationInfluence", y, weights, as.double(subset));
a4 <- libcoin:::.libcoinCall("R_ExpectationInfluence", y, as.double(weights), subset);
a5 <- LinStatExpCov(x, y, weights = weights, subset = subset)$ExpectationInfluence

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, a5))


###################################################
### code chunk number 16: CovarianceInfluence
###################################################
sumweights <- sum(weights[subset])
yc <- t(t(y) - expecty)
r1y <- rep(1:ncol(y), ncol(y))
r2y <- rep(1:ncol(y), each = ncol(y))
a0 <- colSums(yc[subset, r1y] * yc[subset, r2y] * weights[subset]) / sumweights
a0 <- matrix(a0, ncol = ncol(y))
vary <- diag(a0)
a0 <- a0[lower.tri(a0, diag = TRUE)]
a1 <- libcoin:::.libcoinCall("R_CovarianceInfluence", y, weights, subset, 0L);
a2 <- libcoin:::.libcoinCall("R_CovarianceInfluence", y, as.double(weights), as.double(subset), 0L);
a3 <- libcoin:::.libcoinCall("R_CovarianceInfluence", y, weights, as.double(subset), 0L);
a4 <- libcoin:::.libcoinCall("R_CovarianceInfluence", y, as.double(weights), subset, 0L);
a5 <- LinStatExpCov(x, y, weights = weights, subset = subset)$CovarianceInfluence

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, a5))

a1 <- libcoin:::.libcoinCall("R_CovarianceInfluence", y, weights, subset, 1L);
a2 <- libcoin:::.libcoinCall("R_CovarianceInfluence", y, as.double(weights), as.double(subset), 1L);
a3 <- libcoin:::.libcoinCall("R_CovarianceInfluence", y, weights, as.double(subset), 1L);
a4 <- libcoin:::.libcoinCall("R_CovarianceInfluence", y, as.double(weights), subset, 1L);
a5 <- LinStatExpCov(x, y, weights = weights, subset = subset, varonly = TRUE)$VarianceInfluence
a0 <- vary

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, a5))


###################################################
### code chunk number 17: ExpectationCovarianceX
###################################################
a0 <- colSums(x[subset, ] * weights[subset]) 
a0
a1 <- libcoin:::.libcoinCall("R_ExpectationX", x, P, weights, subset);
a2 <- libcoin:::.libcoinCall("R_ExpectationX", x, P, as.double(weights), as.double(subset));
a3 <- libcoin:::.libcoinCall("R_ExpectationX", x, P, weights, as.double(subset));
a4 <- libcoin:::.libcoinCall("R_ExpectationX", x, P, as.double(weights), subset);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECVxyws$ExpectationX))

a0 <- colSums(x[subset, ]^2 * weights[subset]) 
a1 <- libcoin:::.libcoinCall("R_CovarianceX", x, P, weights, subset, 1L);
a2 <- libcoin:::.libcoinCall("R_CovarianceX", x, P, as.double(weights), as.double(subset), 1L);
a3 <- libcoin:::.libcoinCall("R_CovarianceX", x, P, weights, as.double(subset), 1L);
a4 <- libcoin:::.libcoinCall("R_CovarianceX", x, P, as.double(weights), subset, 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a0 <- as.vector(colSums(Xfactor[subset, ] * weights[subset]))
a0
a1 <- libcoin:::.libcoinCall("R_ExpectationX", ix, Lx, weights, subset);
a2 <- libcoin:::.libcoinCall("R_ExpectationX", ix, Lx, as.double(weights), as.double(subset));
a3 <- libcoin:::.libcoinCall("R_ExpectationX", ix, Lx, weights, as.double(subset));
a4 <- libcoin:::.libcoinCall("R_ExpectationX", ix, Lx, as.double(weights), subset);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a1 <- libcoin:::.libcoinCall("R_CovarianceX", ix, Lx, weights, subset, 1L);
a2 <- libcoin:::.libcoinCall("R_CovarianceX", ix, Lx, as.double(weights), as.double(subset), 1L);
a3 <- libcoin:::.libcoinCall("R_CovarianceX", ix, Lx, weights, as.double(subset), 1L);
a4 <- libcoin:::.libcoinCall("R_CovarianceX", ix, Lx, as.double(weights), subset, 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

r1x <- rep(1:ncol(Xfactor), ncol(Xfactor))
r2x <- rep(1:ncol(Xfactor), each = ncol(Xfactor))
a0 <- colSums(Xfactor[subset, r1x] * Xfactor[subset, r2x] * weights[subset])
a0 <- matrix(a0, ncol = ncol(Xfactor))
vary <- diag(a0)
a0 <- a0[lower.tri(a0, diag = TRUE)]

a1 <- libcoin:::.libcoinCall("R_CovarianceX", ix, Lx, weights, subset, 0L)
a2 <- libcoin:::.libcoinCall("R_CovarianceX", ix, Lx, as.double(weights), as.double(subset), 0L)
a3 <- libcoin:::.libcoinCall("R_CovarianceX", ix, Lx, weights, as.double(subset), 0L)
a4 <- libcoin:::.libcoinCall("R_CovarianceX", ix, Lx, as.double(weights), subset, 0L)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))


###################################################
### code chunk number 18: SimpleSums
###################################################
a0 <- sum(weights[subset])
a1 <- libcoin:::.libcoinCall("R_Sums", N, weights, subset)
a2 <- libcoin:::.libcoinCall("R_Sums", N, as.double(weights), as.double(subset))
a3 <- libcoin:::.libcoinCall("R_Sums", N, weights, as.double(subset))
a4 <- libcoin:::.libcoinCall("R_Sums", N, as.double(weights), subset)
stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))


###################################################
### code chunk number 19: KronSums
###################################################
r1 <- rep(1:ncol(x), ncol(y))
r2 <- rep(1:ncol(y), each = ncol(x))

a0 <- colSums(x[subset,r1] * y[subset,r2] * weights[subset])
a1 <- libcoin:::.libcoinCall("R_KronSums", x, P, y, weights, subset, 0L)
a2 <- libcoin:::.libcoinCall("R_KronSums", x, P, y, as.double(weights),
as.double(subset), 0L)
a3 <- libcoin:::.libcoinCall("R_KronSums", x, P, y, weights,
as.double(subset), 0L)
a4 <- libcoin:::.libcoinCall("R_KronSums", x, P, y, as.double(weights),
subset, 0L)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a0 <- as.vector(colSums(Xfactor[subset,r1Xfactor] * 
                        y[subset,r2Xfactor] * weights[subset]))
a1 <- libcoin:::.libcoinCall("R_KronSums", ix, Lx, y, weights, subset, 0L)
a2 <- libcoin:::.libcoinCall("R_KronSums", ix, Lx, y, as.double(weights),
as.double(subset), 0L)
a3 <- libcoin:::.libcoinCall("R_KronSums", ix, Lx, y, weights,
as.double(subset), 0L)
a4 <- libcoin:::.libcoinCall("R_KronSums", ix, Lx, y, as.double(weights),
subset, 0L)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))



###################################################
### code chunk number 20: KronSums-Permutation
###################################################
a0 <- colSums(x[subset,r1] * y[subsety, r2])
a1 <- libcoin:::.libcoinCall("R_KronSums_Permutation", x, P, y, subset, subsety)
a2 <- libcoin:::.libcoinCall("R_KronSums_Permutation", x, P, y, as.double(subset), as.double(subsety))
stopifnot(all.equal(a0, a1) && all.equal(a0, a1))

a0 <- as.vector(colSums(Xfactor[subset,r1Xfactor] * y[subsety, r2Xfactor]))
a1 <- libcoin:::.libcoinCall("R_KronSums_Permutation", ix, Lx, y, subset, subsety)
a1 <- libcoin:::.libcoinCall("R_KronSums_Permutation", ix, Lx, y, as.double(subset), as.double(subsety))
stopifnot(all.equal(a0, a1))


###################################################
### code chunk number 21: colSums
###################################################
a0 <- colSums(x[subset,] * weights[subset])
a1 <- libcoin:::.libcoinCall("R_colSums", x, weights, subset)
a2 <- libcoin:::.libcoinCall("R_colSums", x, as.double(weights), as.double(subset))
a3 <- libcoin:::.libcoinCall("R_colSums", x, weights, as.double(subset))
a4 <- libcoin:::.libcoinCall("R_colSums", x, as.double(weights), subset)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))


###################################################
### code chunk number 22: OneTableSum
###################################################

a0 <- as.vector(xtabs(weights ~ ixf, subset = subset))
a1 <- ctabs(ix, weights = weights, subset = subset)[-1]
a2 <- ctabs(ix, weights = as.double(weights), subset = as.double(subset))[-1]
a3 <- ctabs(ix, weights = weights, subset = as.double(subset))[-1]
a4 <- ctabs(ix, weights = as.double(weights), subset = subset)[-1]

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))


###################################################
### code chunk number 23: TwoTableSum
###################################################

a0 <- xtabs(weights ~ ixf + iyf, subset = subset)
class(a0) <- "matrix"
dimnames(a0) <- NULL
attributes(a0)$call <- NULL
a1 <- ctabs(ix, iy, weights = weights, subset = subset)[-1, -1]
a2 <- ctabs(ix, iy, weights = as.double(weights), 
            subset = as.double(subset))[-1, -1]
a3 <- ctabs(ix, iy, weights = weights, subset = as.double(subset))[-1, -1]
a4 <- ctabs(ix, iy, weights = as.double(weights), subset = subset)[-1, -1]

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))


###################################################
### code chunk number 24: ThreeTableSum
###################################################
a0 <- xtabs(weights ~ ixf + iyf + block, subset = subset)
class(a0) <- "array"
dimnames(a0) <- NULL
attributes(a0)$call <- NULL
a1 <- ctabs(ix, iy, block, weights, subset)[-1, -1,]
a2 <- ctabs(ix, iy, block, as.double(weights), as.double(subset))[-1,-1,]
a3 <- ctabs(ix, iy, block, weights, as.double(subset))[-1,-1,]
a4 <- ctabs(ix, iy, block, as.double(weights), subset)[-1,-1,]

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))


###################################################
### code chunk number 25: blocks
###################################################
sb <- sample(block)
ns1 <- do.call("c", tapply(subset, sb[subset], function(i) i))
ns2 <- libcoin:::.libcoinCall("R_order_subset_wrt_block", y, integer(0), subset, sb)
all.equal(ns1, ns2, check.attributes = FALSE)


