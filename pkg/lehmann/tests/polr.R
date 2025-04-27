
library("lehmann")
library("MASS")
library("tram")

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

lehmann::Lehmann(y = y, x = x)
lehmann::Lehmann(y = y, x = x, B = 10000)

ci <- confint(m <- polr(y ~ x, method = "cloglog"))
c(coef(m)["xB"], ci)

m <- Polr(y ~ x, method = "cloglog")
score_test(m, parm = "xB")

lehmann::Lehmann(y = y, x = x, type = "HazardRatio")
lehmann::Lehmann(y = y, x = x, type = "HazardRatio", B = 10000)

ci <- confint(m <- polr(y ~ x, method = "loglog"))
c(coef(m)["xB"], ci)

m <- Polr(y ~ x, method = "loglog")
score_test(m, parm = "xB")

lehmann::Lehmann(y = y, x = x, type = "Lehmann")
lehmann::Lehmann(y = y, x = x, type = "Lehmann", B = 10000)

