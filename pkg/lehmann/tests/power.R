
library("lehmann")

set.seed(29)
alt <- "two"

N <- 25
delta <- .759

(p <- power.trafo.test(n = N, delta = delta, link = "logit", alt = alt))

power.trafo.test(power = p$power, delta = delta, link = "logit", alt = alt)

power.trafo.test(n = p$n, power = p$power, link = "logit", alt = alt)

power.trafo.test(n = p$n, power = p$power, delta = p$delta, sig.level =
NULL, link = "logit", alt = alt)

power.trafo.test(n = N, delta = -delta, link = "logit", alt = alt)

K <- 5
prob <- runif(K)
prob <- prob / sum(prob)

(p <- power.trafo.test(n = N, prob = prob, delta = delta, link = "logit", alt = alt))

K <- 10
prob <- runif(K)
prob <- prob / sum(prob)

(p <- power.trafo.test(n = N, prob = prob, delta = delta, link = "logit", alt = alt))

(p <- power.trafo.test(n = N, aratio = 2, prob = prob, delta = delta, link = "logit", alt = alt))


power.prop.test(n = 25, p1 = .5, p2 = plogis(qlogis(.5) - delta))
power.trafo.test(n = 25, prob = c(.5, .5), delta = delta)
