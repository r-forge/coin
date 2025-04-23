
library("lehmann")
library("rms")
library("survival")

set.seed(29)

N <- 500
x <- gl(2, N)
y <- rlogis(length(x), location = c(0, 2)[x])

ci <- confint(m <- orm(y ~ x))
c(rev(coef(m))[1], ci[nrow(ci),])

Lehmann(y = y, x = x)
Lehmann(y = y, x = x, nbins = 99)
Lehmann(y = y, x = x, nbins = 100)

-coef(m <- coxph(Surv(y, rep(TRUE, length(y))) ~ x))
-rev(confint(m))
### rms: cloglog for P(Y < y) is loglog for P(Y >= y), the latter
### being used in orm()
ci <- confint(m <- orm(y ~ x, family = "loglog"))
c(rev(coef(m))[1], ci[nrow(ci),])

Lehmann(y = y, x = x, type = "HazardRatio")
Lehmann(y = y, x = x, type = "HazardRatio", nbins = 99)
Lehmann(y = y, x = x, type = "HazardRatio", nbins = 100)

ci <- confint(m <- orm(y ~ x, family = "cloglog"))
c(rev(coef(m))[1], ci[nrow(ci),])

Lehmann(y = y, x = x, type = "Lehmann")
Lehmann(y = y, x = x, type = "Lehmann", nbins = 99)
Lehmann(y = y, x = x, type = "Lehmann", nbins = 100)


