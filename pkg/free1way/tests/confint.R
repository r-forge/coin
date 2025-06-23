
library("free1way")

set.seed(29)

level <- .9
N <- 15
w <- gl(2, N)
beta <- -1
y <- rlogis(length(w), location = c(0, beta)[w])

ft <- free1way.test(y ~ w)
coef(ft)
logLik(ft)
(ci <- confint(ft, test = "Wald", level = level))
summary(ft1 <- free1way.test(y ~ w, mu = ci[1]), test = "Wald", alternative =
"greater")$p.value
summary(ft2 <- free1way.test(y ~ w, mu = ci[2]), test = "Wald", alternative =
"less")$p.value
summary(ft12 <- free1way.test(y ~ w, mu = ci[1]), test = "Wald")$p.value
summary(ft22 <- free1way.test(y ~ w, mu = ci[2]), test = "Wald")$p.value
coef(ft1) + ci[1]
coef(ft2) + ci[2]
logLik(ft1)
logLik(ft2)

all.equal(ci, confint(ft1, test = "Wald", level = level) + ci[1])
all.equal(ci, confint(ft2, test = "Wald", level = level) + ci[2])
all.equal(ci, confint(ft12, test = "Wald", level = level) + ci[1])
all.equal(ci, confint(ft22, test = "Wald", level = level) + ci[2])


(ci <- confint(ft, test = "LRT", level = level))
summary(ft12 <- free1way.test(y ~ w, mu = ci[1]), test = "LRT")$p.value
summary(ft22 <- free1way.test(y ~ w, mu = ci[2]), test = "LRT")$p.value

all.equal(ci, confint(ft12, test = "LRT", level = level) + ci[1])
all.equal(ci, confint(ft22, test = "LRT", level = level) + ci[2])


(ci <- confint(ft, test = "Rao", level = level))
summary(ft12 <- free1way.test(y ~ w, mu = ci[1]), test = "Rao")$p.value
summary(ft22 <- free1way.test(y ~ w, mu = ci[2]), test = "Rao")$p.value

all.equal(ci, confint(ft12, test = "Rao", level = level) + ci[1])
all.equal(ci, confint(ft22, test = "Rao", level = level) + ci[2])


(ci <- confint(ft, test = "Permutation", level = .9))
summary(ft1 <- free1way.test(y ~ w, mu = ci[1]), test = "Permutation", alternative = "greater")
summary(ft2 <- free1way.test(y ~ w, mu = ci[2]), test = "Permutation", alternative = "less")
summary(ft12 <- free1way.test(y ~ w, mu = ci[1]), test = "Permutation")
summary(ft22 <- free1way.test(y ~ w, mu = ci[2]), test = "Permutation")
coef(ft1) + ci[1]
coef(ft2) + ci[2]
logLik(ft1)
logLik(ft2)

confint(ft1, test = "Permutation", level = level) + ci[1]
confint(ft2, test = "Permutation", level = level) + ci[2]
confint(ft12, test = "Permutation", level = level) + ci[1]
confint(ft22, test = "Permutation", level = level) + ci[2]

ftc <- free1way.test(y ~ w, mu = coef(ft))
summary(ftc, test = "Permutation", alternative = "greater")
summary(ftc, test = "Permutation", alternative = "less")
summary(ftc, test = "Permutation")

confint(ftc, test = "Permutation", level = level) + coef(ft)

ft <- free1way.test(y ~ w, B = 1000)
(ci <- confint(ft, test = "Permutation", level = .9))
summary(ft1 <- free1way.test(y ~ w, mu = ci[1]), test = "Permutation", alternative = "greater")
summary(ft2 <- free1way.test(y ~ w, mu = ci[2]), test = "Permutation", alternative = "less")
summary(ft12 <- free1way.test(y ~ w, mu = ci[1]), test = "Permutation")
summary(ft22 <- free1way.test(y ~ w, mu = ci[2]), test = "Permutation")
coef(ft1) + ci[1]
coef(ft2) + ci[2]
logLik(ft1)
logLik(ft2)

confint(ft1, test = "Permutation", level = level) + ci[1]
confint(ft2, test = "Permutation", level = level) + ci[2]
confint(ft12, test = "Permutation", level = level) + ci[1]
confint(ft22, test = "Permutation", level = level) + ci[2]

summary(ftg <- free1way.test(y ~ w, mu = coef(ft)), test = "Permutation", alternative = "greater")
summary(ftl <- free1way.test(y ~ w, mu = coef(ft)), test = "Permutation", alternative = "less")
summary(ftt <- free1way.test(y ~ w, mu = coef(ft)), test = "Permutation")

ftc <- free1way.test(y ~ w, mu = coef(ft), B = 1000)
summary(ftc, test = "Permutation", alternative = "greater")
summary(ftc, test = "Permutation", alternative = "less")
summary(ftc, test = "Permutation")

confint(ftc, test = "Permutation", level = level) + coef(ft)
