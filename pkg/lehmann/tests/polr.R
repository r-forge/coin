
library("lehmann")
library("MASS")
library("tram")
library("coin")

set.seed(29)

N <- 500
x <- gl(2, N, labels = c("A", "B"))
y <- rlogis(length(x), location = c(0, 2)[x])
y <- cut(y, breaks = c(-Inf, -1, 0, 1, Inf), ordered = TRUE)

ci <- confint(m <- polr(y ~ x, method = "logistic"))
c(coef(m)["xB"], ci)
#score_test(m, parm = "xB")

m <- Polr(y ~ x, method = "logistic")
score_test(m, parm = "xB")
trafo.test(y = y, x = x)
trafo.test(y = y, x = x, inference = "MLScore")
perm_test(m, parm = "xB", distribution = "asymptotic")
trafo.test(y = y, x = x, inference = "PermScore")
perm_test(m, parm = "xB", distribution = approximate(nresample = 10000))
trafo.test(y = y, x = x, inference = "PermScore", B = 10000)

ci <- confint(m <- polr(y ~ x, method = "cloglog"))
c(coef(m)["xB"], ci)

m <- Polr(y ~ x, method = "cloglog")
score_test(m, parm = "xB")

trafo.test(y = y, x = x, link = "cloglog")
trafo.test(y = y, x = x, link = "cloglog", inference = "MLScore")
trafo.test(y = y, x = x, link = "cloglog", inference = "PermScore", B = 10000)

ci <- confint(m <- polr(y ~ x, method = "loglog"))
c(coef(m)["xB"], ci)

m <- Polr(y ~ x, method = "loglog")
score_test(m, parm = "xB")

trafo.test(y = y, x = x, link = "loglog")
trafo.test(y = y, x = x, link = "loglog", inference = "MLScore")
trafo.test(y = y, x = x, link = "loglog", inference = "PermScore", B = 10000)

m <- Polr(y ~ x, method = "probit")
score_test(m, parm = "xB")

trafo.test(y = y, x = x, link = "probit")
trafo.test(y = y, x = x, link = "probit", inference = "MLScore")
trafo.test(y = y, x = x, link = "probit", inference = "PermScore", B = 10000)
