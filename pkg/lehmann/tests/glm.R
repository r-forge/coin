
library("lehmann")

set.seed(29)

N <- 500
x <- gl(2, N, labels = c("A", "B"))
y <- rlogis(length(x), location = c(0, 2)[x])
y <- factor(y < 0)

ci <- confint(m <- glm(y ~ x, family = binomial()))
c(coef(m)["xB"], ci["xB",])

trafo.test(y = y, x = x)
trafo.test(y = y, x = x, inference = "MLScore")
trafo.test(y = y, x = x, inference = "PermScore", B = 10000)

ci <- confint(m <- glm(y ~ x, family = binomial(link = "cloglog")))
c(coef(m)["xB"], ci["xB",])

trafo.test(y = y, x = x, type = "Savage")
trafo.test(y = y, x = x, type = "Savage", inference = "MLScore")
trafo.test(y = y, x = x, type = "Savage", inference = "PermScore", B = 10000)

y <- relevel(y, levels(y)[2])
ci <- confint(m <- glm(y ~ x, family = binomial(link = "cloglog")))
c(coef(m)["xB"], ci["xB",])

trafo.test(y = y, x = x, type = "Lehmann")
trafo.test(y = y, x = x, type = "Lehmann", inference = "MLScore")
trafo.test(y = y, x = x, type = "Lehmann", inference = "PermScore", B = 10000)

ci <- confint(m <- glm(y ~ x, family = binomial(link = "cauchit")))
c(coef(m)["xB"], ci["xB",])

trafo.test(y = y, x = x, type = "Cauchy")
trafo.test(y = y, x = x, type = "Cauchy", inference = "MLScore")
trafo.test(y = y, x = x, type = "Cauchy", inference = "PermScore", B = 10000)
