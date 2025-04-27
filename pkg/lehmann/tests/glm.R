
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

lehmann::Lehmann(y = y, x = x)
# lehmann::Lehmann(y = y, x = x, B = 10000)

ci <- confint(m <- glm(y ~ x, family = binomial(link = "cloglog")))
c(coef(m)["xB"], ci["xB",])

Lehmann(y = y, x = x, type = "HazardRatio")
# Lehmann(y = y, x = x, type = "HazardRatio", B = 10000)

y <- relevel(y, levels(y)[2])
ci <- confint(m <- glm(y ~ x, family = binomial(link = "cloglog")))
c(coef(m)["xB"], ci["xB",])

Lehmann(y = y, x = x, type = "Lehmann")
# Lehmann(y = y, x = x, type = "Lehmann", B = 10000)

