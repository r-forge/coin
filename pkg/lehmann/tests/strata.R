
library("lehmann")
library("tram")
library("coin")

set.seed(29)

w <- gl(2, 1)
s <- gl(3, 1)
y <- gl(2, 1)
d <- expand.grid(y = y, s = s, w = w)
d$Freq <- sample(nrow(d)) + 100

x <- xtabs(Freq ~ y + w + s, data = d)

m <- glm(y ~ w + s, data = d, weights = d$Freq, family = binomial())
summary(m)
exp(confint(profile(m, type = "LRT"))["w2",])
exp(confint(profile(m, type = "Rao"))["w2",])
logLik(m)

mantelhaen.test(x, correct = FALSE)

tt <- trafo.test(x)
(tt$statistic)^2
tt$statistic
tt$p.value
exp(tt$conf.int)
exp(tt$estimate)
-tt$neglogLik

tt <- trafo.test(x, inference = "LR")
(tt$statistic)^2
tt$statistic
tt$p.value
exp(tt$conf.int)
exp(tt$estimate)


tt <- trafo.test(x, inference = "MLS")
(tt$statistic)^2
tt$statistic
tt$p.value
exp(tt$conf.int)
exp(tt$estimate)



independence_test(y ~ w | s, data = d, weights = ~ Freq)


tt <- trafo.test(x, inference = "PermScore")
(tt$statistic)^2
(tt$statistic)
tt$p.value
exp(tt$conf.int)
exp(tt$estimate)

mantelhaen.test(x, correct = FALSE, exact = TRUE)


independence_test(y ~ w | s, data = d, weights = ~ Freq, distribution = approximate(100000))


tt <- trafo.test(x, inference = "PermScore", B = 100000)
(tt$statistic)^2
tt$statistic
tt$p.value
exp(tt$conf.int)
exp(tt$estimate)

independence_test(y ~ w | s, data = d, weights = ~ Freq, distribution = "exact")

w <- gl(2, 1)
s <- gl(3, 1)
N <- 10
d <- expand.grid(1:N, s = s, w = w)
d$y <- rlogis(nrow(d))

(it <- independence_test(y ~ w | s, data = d,
    ytrafo = function(data) 
        trafo(data, numeric_trafo = rank_trafo, block = d$s)))

(tt <- trafo.test(y ~ w + s, data = d, inference = "PermScore"))
-tt$neglogLik

independence_test(y ~ w | s, data = d,
    distribution = approximate(10000),
    ytrafo = function(data) 
        trafo(data, numeric_trafo = rank_trafo, block = d$s))


(tt <- trafo.test(y ~ w + s, data = d, inference = "PermScore", 
    B = 10000))

statistic(it, "linear")

r <- unlist(tapply(d$y, d$s, rank))
sum(r * (1:0)[d$w])

independence_test(y ~ w | s, data = d, distribution = "exact",
    ytrafo = function(data) 
        trafo(data, numeric_trafo = rank_trafo, block = d$s))

