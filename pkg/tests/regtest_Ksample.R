
### Regression tests for the K sample problem, i.e.,
### testing the independence of a numeric variable
### `y' and a factor `x' (possibly blocked)

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### generate data
dat <- data.frame(x = gl(4, 25), y = rnorm(100), block = gl(5, 20))[sample(1:100, 50),]


### Kruskal-Wallis Test

### asymptotic distribution
ptwo <- kruskal.test(y ~ x, data = dat)$p.value

stopifnot(isequal(pvalue(kruskal_test(y ~ x, data = dat)), ptwo))

stopifnot(isequal(pvalue(oneway_test(y ~ x, data = dat, distribution = "asympt",
                                   teststat = "quad",
    ytrafo = function(data) trafo(data, numeric_trafo = rank_trafo))), ptwo))

### approximated distribution
rtwo <- pvalue(kruskal_test(y ~ x, data = dat, distribution = "approx")) / ptwo
stopifnot(all(rtwo > 0.9 &
              rtwo < 1.1))

### <FIXME> add block examples </FIXME>

### sanity checks
try(kruskal_test(x ~ y, data = dat))
try(kruskal_test(x ~ y | y, data = dat))


### Fligner Test

### asymptotic distribution
ptwo <- fligner.test(y ~ x, data = dat)$p.value

stopifnot(isequal(pvalue(fligner_test(y ~ x, data = dat)), ptwo))

dat$yy <- dat$y - tapply(dat$y, dat$x, median)[dat$x]
stopifnot(isequal(pvalue(oneway_test(yy ~ x, data = dat, distribution = "asympt",
                                   teststat = "quad",
    ytrafo = function(data) trafo(data, numeric_trafo = fligner_trafo))), ptwo))

### approximated distribution
rtwo <- pvalue(fligner_test(y ~ x, data = dat, distribution = "approx")) / ptwo
stopifnot(all(rtwo > 0.9 &
              rtwo < 1.1))

### <FIXME> add block examples </FIXME>

### sanity checks
try(fligner_test(x ~ y, data = dat))
try(fligner_test(x ~ y | y, data = dat))


### One-way Test
oneway_test(y ~ x, data = dat)

oneway_test(y ~ ordered(x), data = dat)
oneway_test(y ~ ordered(x), data = dat,
            alternative = "less")
oneway_test(y ~ ordered(x), data = dat,
            alternative = "greater")

oneway_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
oneway_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "less")
oneway_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "greater")


### Normal Scores Test
normal_test(y ~ x, data = dat)

normal_test(y ~ ordered(x), data = dat)
normal_test(y ~ ordered(x), data = dat,
            alternative = "less")
normal_test(y ~ ordered(x), data = dat,
            alternative = "greater")

normal_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
normal_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "less")
normal_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "greater")


### Median Test
median_test(y ~ x, data = dat)

median_test(y ~ ordered(x), data = dat)
median_test(y ~ ordered(x), data = dat,
            alternative = "less")
median_test(y ~ ordered(x), data = dat,
            alternative = "greater")

median_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
median_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "less")
median_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "greater")


### Savage Test
savage_test(y ~ x, data = dat)

savage_test(y ~ ordered(x), data = dat)
savage_test(y ~ ordered(x), data = dat,
            alternative = "less")
savage_test(y ~ ordered(x), data = dat,
            alternative = "greater")

savage_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
savage_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "less")
savage_test(y ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
            alternative = "greater")


### Logrank Test
surv_test(Surv(y) ~ x, data = dat)

surv_test(Surv(y) ~ ordered(x), data = dat)
surv_test(Surv(y) ~ ordered(x), data = dat,
          alternative = "less")
surv_test(Surv(y) ~ ordered(x), data = dat,
          alternative = "greater")

surv_test(Surv(y) ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)))
surv_test(Surv(y) ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
          alternative = "less")
surv_test(Surv(y) ~ x, data = dat, scores = list(x = c(2, 4, 6, 8)),
          alternative = "greater")
