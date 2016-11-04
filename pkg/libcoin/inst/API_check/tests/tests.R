
library("dummy")


n <- 100
p <- 4
q <- 2
X <- matrix(runif(p * n), nc = p)
Y <- matrix(runif(q * n), nc = q)

.Call("myR_ExpectationCovarianceStatistic", X, Y, integer(0), integer(0), integer(0), as.integer(FALSE), package = "dummy")
