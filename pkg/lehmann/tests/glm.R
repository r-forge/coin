
library("lehmann")
#library("tram")

set.seed(29)

N <- 500
x <- gl(2, N, labels = c("A", "B"))
y <- rlogis(length(x), location = c(0, 2)[x])
y <- factor(y < 0)

ci <- confint(m <- glm(y ~ x, family = binomial()))
c(coef(m)["xB"], ci["xB",])
#score_test(m, parm = "xB")

#m <- Polr(y ~ x, method = "logistic")
#score_test(m, parm = "xB")

lehmann::trafo.test(y = y, x = x)
lehmann::trafo.test(y = y, x = x, B = 10000)

ci <- confint(m <- glm(y ~ x, family = binomial(link = "cloglog")))
c(coef(m)["xB"], ci["xB",])

trafo.test(y = y, x = x, type = "Savage")
trafo.test(y = y, x = x, type = "Savage", B = 10000)

y <- relevel(y, levels(y)[2])
ci <- confint(m <- glm(y ~ x, family = binomial(link = "cloglog")))
c(coef(m)["xB"], ci["xB",])

trafo.test(y = y, x = x, type = "Lehmann")
trafo.test(y = y, x = x, type = "Lehmann", B = 10000)

