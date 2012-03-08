

flies <- data.frame(wl = c(1.72, 1.64, 1.74, 1.70, 1.82, 1.82, 1.90, 1.82,
                           2.08, 1.78, 1.86, 1.96, 2.00, 2.00, 1.96),
                    al = c(1.24, 1.38, 1.36, 1.40, 1.38, 1.48, 1.38, 1.54, 1.56,
                           1.14, 1.20, 1.30, 1.26, 1.28, 1.18),
                    species = factor(c(rep("Af", 9), rep("Apf", 6))))

library("coin")

it <- independence_test(wl ~ species, data = flies, distr = exact())
it <- independence_test(al ~ species, data = flies, distr = exact())

sit <- support(it)
a <- sapply(sit, function(i) pperm(it, i))
plot(sit, a, type = "S")
curve(pnorm, min(sit), max(sit), add = TRUE, col = "red")

it <- independence_test(al + wl ~ species, data = flies, teststat = "quad",
distribution = approximate(B = 100000))


typeface <- data.frame(speed = c(135, 91, 111, 87, 122, 175, 130, 514, 283,
105, 147, 159, 107, 194), stype = factor(c(rep(1, 5), rep(2, 4), rep(3,
5))))

it <- oneway_test(speed ~ stype, data = typeface, teststat = "quad", dist = approximate(B = 100000))

sit <- unique(round(support(it), 1))
a <- sapply(sit, function(i) pperm(it, i))          
plot(sit, a, type = "S")                
lines(sit, sapply(sit, function(i) pchisq(i, df = 2)), col = "red")

summary(aov(speed ~ stype, data = typeface))
