
library("lehmann")
library("tram")

set.seed(29)

N <- 10
w <- gl(2, 1)
s <- gl(5, 1)
y <- gl(4, 1, ordered = TRUE)
d <- expand.grid(w = w, s = s, y = y)
d$W <- sample(1:nrow(d), N)

x <- xtabs(W ~ y + w + s, data = d)
trafo.test(x)
summary(m <- Polr(y | s ~ w, data = d, weights = d$W))

trafo.test(x, inference = "LR")
confint(profile(m))

trafo.test(x, inference = "MLS")
score_test(m, parm = "w2")


