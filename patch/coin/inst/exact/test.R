
dyn.load("2s.so")
library(coin)
library(exactRankTests)

x <- 1:10 
storage.mode(x) <- "double"
m <- 5
storage.mode(m) <- "integer"
obs <- sum(x) - 20
storage.mode(obs) <- "double"

.Call("R2sample", x, m, obs)
pwilcox(20, 5, 5)

x <- normalfun(rnorm(10))
storage.mode(x) <- "double"
m <- 5
storage.mode(m) <- "integer"
obs <- 2
storage.mode(obs) <- "double"
.Call("R2sample", x, m, obs)

x <- normalfun(rnorm(30))
storage.mode(x) <- "double"
m <- 10
storage.mode(m) <- "integer"
obs <- 2
storage.mode(obs) <- "double"
.Call("R2sample", x, m, obs)

x <- 1:50
storage.mode(x) <- "double"
m <- 10
storage.mode(m) <- "integer"
obs <- 100
storage.mode(obs) <- "double"
.Call("R2sample", x, m, obs)
.Call("R2sample", x, m, as.double(150))
.Call("R2sample", x, m, as.double(210))

pperm(100, x, 10)
pperm(150, x, 10)
pperm(210, x, 10)
