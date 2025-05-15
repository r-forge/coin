
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


example(mantelhaen.test)

mantelhaen.test(Rabbits, correct = FALSE)
trafo.test(as.table(Rabbits))

x <- as.table(matrix(c(1, 1, 0, 0, 1, 1,
                       1, 0, 0, 0, 2, 1), ncol = 2))
trafo.test(x)

trafo.test(x[-c(3, 4),])

library("coin")

N <- 50
s <- gl(3, 1)
d <- expand.grid(1:N, w = w, s = s)
d$y <- as.ordered(d$yi <- sample(1:3, nrow(d), replace = TRUE))

logLik(m <- Polr(y | 0 + s ~ w, data = d))
confint(profile(m))

trafo.test(y ~ w + s, data = d, inference = "LR")


score_test(m)
trafo.test(y ~ w + s, data = d, inference = "MLScore")

perm_test(m)
(tt <- trafo.test(y ~ w + s, data = d, inference = "PermScore"))
-tt$neglogLik
wilcox_test(yi ~ w | s, data = d)

perm_test(m, distribution = approximate(10000))

(tt <- trafo.test(y ~ w + s, data = d, inference = "PermScore", B = 10000))
-tt$neglogLik
wilcox_test(yi ~ w | s, data = d, distribution = approximate(10000))

