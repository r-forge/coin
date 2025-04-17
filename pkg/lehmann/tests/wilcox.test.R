
library("lehmann")
library("rms")

set.seed(29)

N <- 5000
x <- gl(2, N)
y <- rlogis(length(x), location = c(0, 2)[x])

system.time(ci <- confint(m <- orm(y ~ x)))
c(rev(coef(m))[1], ci[nrow(ci),])

system.time(ci <- Lehmann(y = y, x = x))
ci

system.time(ci <- Lehmann(y = y, x = x, nbins = 100))
ci

