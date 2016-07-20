
library("libcoin")
library("coin")

set.seed(29)
n <- 100
B <- 1000
x <- round(runif(n), 2)
y <- rnorm(n, mean = x < .5, sd = 2.6)
y2 <- runif(n)
blk <- gl(4, n/4)
ux <- sort(unique(x))
ix <- unclass(cut(x, breaks = c(-Inf, ux[-length(ux)] + diff(ux) / 2, Inf)))

(mt <- maxstat_test(y ~ x , distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, y, B = 1000)
(tst <- Test(lev, xtrafo = "maxstat", type = "max"))
ux[tst$index]
(tst <- Test(lev, xtrafo = "maxstat", type = "quad"))
ux[tst$index]
lev <- LinStatExpCov(ix, y, B = 1000, varonly = TRUE)
(tst <- Test(lev, xtrafo = "maxstat", type = "max"))
ux[tst$index]
(tst <- Test(lev, xtrafo = "maxstat", type = "quad"))
ux[tst$index]

(mt <- maxstat_test(y ~ x | blk, distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, y, block = blk, B = 1000)
(tst <- Test(lev, xtrafo = "maxstat", type = "max"))
ux[tst$index]
(tst <- Test(lev, xtrafo = "maxstat", type = "quad"))
ux[tst$index]
lev <- LinStatExpCov(ix, y, block = blk, B = 1000, varonly = TRUE)
try(tst <- Test(lev, xtrafo = "maxstat", type = "max"))

(mt <- maxstat_test(y + y2 ~ x , distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, cbind(y, y2), B = 1000)
(tst <- Test(lev, xtrafo = "maxstat", type = "max"))
ux[tst$index]
(tst <- Test(lev, xtrafo = "maxstat", type = "quad"))
ux[tst$index]
lev <- LinStatExpCov(ix, cbind(y, y2), B = 1000, varonly = TRUE)
(tst <- Test(lev, xtrafo = "maxstat", type = "max"))
ux[tst$index]
(tst <- Test(lev, xtrafo = "maxstat", type = "quad"))
ux[tst$index]

(mt <- maxstat_test(y + y2 ~ x | blk, distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, cbind(y, y2), block = blk, B = 1000)
(tst <- Test(lev, xtrafo = "maxstat", type = "max"))
ux[tst$index]
(tst <- Test(lev, xtrafo = "maxstat", type = "quad"))
ux[tst$index]
lev <- LinStatExpCov(ix, cbind(y, y2), block = blk, B = 1000, varonly = TRUE)
try(tst <- Test(lev, xtrafo = "maxstat", type = "max"))

x <- sample(gl(5, n))
y <- rnorm(length(x), mean = x %in% levels(x)[c(1, 3, 5)], sd = 4.5)
y2 <- runif(length(x))
ix <- unclass(x)
blk <- gl(5, n)

(mt <- maxstat_test(y ~ x , distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, y, B = 1000)
(tst <- Test(lev, xtrafo = "maxstat", type = "max", ordered = FALSE))
(tst <- Test(lev, xtrafo = "maxstat", type = "quad", ordered = FALSE))
lev <- LinStatExpCov(ix, y, B = 1000, varonly = TRUE)
(tst <- Test(lev, xtrafo = "maxstat", type = "max", ordered = FALSE))
(tst <- Test(lev, xtrafo = "maxstat", type = "quad", ordered = FALSE))

(mt <- maxstat_test(y ~ x | blk, distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, y, block = blk, B = 1000)
(tst <- Test(lev, xtrafo = "maxstat", type = "max", ordered = FALSE))
(tst <- Test(lev, xtrafo = "maxstat", type = "quad", ordered = FALSE))
lev <- LinStatExpCov(ix, y, block = blk, B = 1000, varonly = TRUE)
try(tst <- Test(lev, xtrafo = "maxstat", type = "max", ordered = FALSE))

(mt <- maxstat_test(y + y2 ~ x , distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, cbind(y, y2), B = 1000)
(tst <- Test(lev, xtrafo = "maxstat", type = "max", ordered = FALSE))
(tst <- Test(lev, xtrafo = "maxstat", type = "quad", ordered = FALSE))
lev <- LinStatExpCov(ix, cbind(y, y2), B = 1000, varonly = TRUE)
(tst <- Test(lev, xtrafo = "maxstat", type = "max", ordered = FALSE))
(tst <- Test(lev, xtrafo = "maxstat", type = "quad", ordered = FALSE))

(mt <- maxstat_test(y + y2 ~ x | blk, distrib = approximate(B = 1000)))
lev <- LinStatExpCov(ix, cbind(y, y2), block = blk, B = 50)
(tst <- Test(lev, xtrafo = "maxstat", type = "max", ordered = FALSE))
(tst <- Test(lev, xtrafo = "maxstat", type = "quad", ordered = FALSE))

xx <- factor(x == levels(x)[tst$index == 1])
(it <- independence_test(y + y2 ~ xx | blk, teststat = "quad"))

lev <- LinStatExpCov(ix, cbind(y, y2), block = blk, B = 1000, varonly = TRUE)
try(tst <- Test(lev, xtrafo = "maxstat", type = "max", ordered = FALSE))
