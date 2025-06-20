
library("free1way")

set.seed(29)

N <- 500
x <- gl(2, N, labels = c("A", "B"))
y <- rlogis(length(x), location = c(0, 2)[x])
y <- factor(y < 0)

ci <- confint(m <- glm(y ~ x, family = binomial()), test = "LRT")
c(coef(m)["xB"], ci["xB",])
ci <- confint(m <- glm(y ~ x, family = binomial()), test = "Rao")
c(coef(m)["xB"], ci["xB",])

ft <- free1way.test(y = y, x = x, B = 10000)
coef(ft)
print(ft, test = "LRT")
confint(ft, test = "LRT")
print(ft, test = "Rao")
confint(ft, test = "Rao")
print(ft, test = "Wald")
confint(ft, test = "Wald")
print(ft, test = "Permutation")
confint(ft, test = "Permutation")

y <- relevel(y, levels(y)[2])
ci <- confint(m <- glm(y ~ x, family = binomial(link = "cloglog")), test = "LRT")
c(coef(m)["xB"], ci["xB",])
ci <- confint(m <- glm(y ~ x, family = binomial(link = "cloglog")), test = "Rao")
c(coef(m)["xB"], ci["xB",])

ft <- free1way.test(y = y, x = x, link = "loglog", B = 10000)
coef(ft)
print(ft, test = "LRT")
confint(ft, test = "LRT")
print(ft, test = "Rao")
confint(ft, test = "Rao")
print(ft, test = "Wald")
confint(ft, test = "Wald")
print(ft, test = "Permutation")
confint(ft, test = "Permutation")
