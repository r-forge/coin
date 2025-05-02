
library("lehmann")
library("rms")
library("survival")

set.seed(29)

N <- 500
x <- gl(2, N)
y <- rlogis(length(x), location = c(0, 2)[x])

ci <- confint(m <- orm(y ~ x))
c(rev(coef(m))[1], ci[nrow(ci),])

d <- data.frame(y = y, x = x)

trafo.test(y = y, x = x)
trafo.test(y ~ x, data = d)
trafo.test(y = y, x = x, nbins = 99)
trafo.test(y = y, x = x, nbins = 100)

-coef(m <- coxph(Surv(y, rep(TRUE, length(y))) ~ x))
-rev(confint(m))
### rms: cloglog for P(Y < y) is loglog for P(Y >= y), the latter
### being used in orm()
ci <- confint(m <- orm(y ~ x, family = "loglog"))
c(rev(coef(m))[1], ci[nrow(ci),])

trafo.test(y = y, x = x, link = "cloglog")
trafo.test(y = y, x = x, link = "cloglog", nbins = 99)
trafo.test(y = y, x = x, link = "cloglog", nbins = 100)

ci <- confint(m <- orm(y ~ x, family = "cloglog"))
c(rev(coef(m))[1], ci[nrow(ci),])

trafo.test(y = y, x = x, link = "loglog")
trafo.test(y = y, x = x, link = "loglog", nbins = 99)
trafo.test(y = y, x = x, link = "loglog", nbins = 100)

### with offset
mu <- 1
off <- (x == levels(x)[2]) * mu

ci <- confint(m <- orm(y ~ x + offset(off)))
c(rev(coef(m))[1], ci[nrow(ci),])

trafo.test(y = y, x = x, mu = 1)
trafo.test(y = y, x = x, mu = 1, nbins = 99)
trafo.test(y = y, x = x, mu = 1, nbins = 100)

### permutations
N <- 15
x <- gl(2, N)
y <- rlogis(length(x), location = c(0, 2)[x])

ci <- confint(m <- orm(y ~ x))
c(rev(coef(m))[1], ci[nrow(ci),])

trafo.test(y = y, x = x)
trafo.test(y = y, x = x, inference = "MLScore")
trafo.test(y = y, x = x, inference = "LRatio")
trafo.test(y = y, x = x, inference = "PermScore", B = 0)
trafo.test(y = y, x = x, inference = "PermScore", B = 10000)
trafo.test(y = y, x = x, inference = "PermScore", B = Inf)
