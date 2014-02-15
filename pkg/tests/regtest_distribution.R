### Regression tests for the distribution functions

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)


### generate independent two-sample data
dta <- data.frame(y = rnorm(20), x = gl(2, 10), b = factor(rep(1:4, 5)),
                  w = rep(1:3, length = 20))
dta$y5 <- round(dta$y, 5)
dta$y3 <- round(dta$y, 3)

### generate paired two-sample data
dta2 <- data.frame(y = as.vector(rbind(abs(dta$y) * (dta$y >= 0),
                                       abs(dta$y) * (dta$y < 0))),
                   x = factor(rep(0:1, length(dta$y)),
                              labels = c("pos", "neg")),
                   b = gl(length(dta$y), 2))
dta2$y5 <- round(dta2$y, 5)
dta2$y3 <- round(dta2$y, 3)


### check 'algorithm = "auto"'

### scores that can't be mapped into integers

### two-sample with block
try(independence_test(y ~ x | b, data = dta,
                      distribution = exact(algo = "auto")))

try(independence_test(y ~ x | b, data = dta,
                      distribution = exact(algo = "shift")))

try(independence_test(y ~ x | b, data = dta,
                      distribution = exact(algo = "split-up")))

### two-sample without block
it <- independence_test(y ~ x, data = dta,
                        distribution = exact(algo = "auto"))
it@distribution@name
pvalue(it)

try(independence_test(y ~ x, data = dta,
                      distribution = exact(algo = "shift")))

it <- independence_test(y ~ x, data = dta,
                        distribution = exact(algo = "split-up"))
it@distribution@name
pvalue(it)

### paired two-sample
try(independence_test(y ~ x | b, data = dta2, paired = TRUE,
                      distribution = exact(algo = "auto")))

try(independence_test(y ~ x | b, data = dta2, paired = TRUE,
                      distribution = exact(algo = "shift")))

try(independence_test(y ~ x | b, data = dta2, paired = TRUE,
                      distribution = exact(algo = "split-up")))

### mapped into integers using 'fact'

### two-sample with block
it <- independence_test(y5 ~ x | b, data = dta,
                        distribution = exact(algo = "auto", fact = 1e5))
it@distribution@name
pvalue(it)

it <- independence_test(y5 ~ x | b, data = dta,
                        distribution = exact(algo = "shift", fact= 1e5))
it@distribution@name
pvalue(it)

try(independence_test(y5 ~ x | b, data = dta,
                        distribution = exact(algo = "split-up")))

### two-sample without block
it <- independence_test(y5 ~ x, data = dta,
                        distribution = exact(algo = "auto", fact = 1e5))
it@distribution@name
pvalue(it)

it <- independence_test(y5 ~ x, data = dta,
                        distribution = exact(algo = "shift", fact= 1e5))
it@distribution@name
pvalue(it)

it <- independence_test(y5 ~ x, data = dta,
                        distribution = exact(algo = "split-up"))
it@distribution@name
pvalue(it)

### paired two-sample
it <- independence_test(y5 ~ x | b, data = dta2, paired = TRUE,
                        distribution = exact(algo = "auto", fact = 1e5))
it@distribution@name
pvalue(it)

it <- independence_test(y5 ~ x | b, data = dta2, paired = TRUE,
                        distribution = exact(algo = "shift", fact = 1e5))
it@distribution@name
pvalue(it)

try(independence_test(y5 ~ x | b, data = dta2, paired = TRUE,
                      distribution = exact(algo = "split-up")))

### automatically mapped into integers

### two-sample with block
it <- independence_test(y3 ~ x | b, data = dta,
                        distribution = exact(algo = "auto"))
it@distribution@name
pvalue(it)

it <- independence_test(y3 ~ x | b, data = dta,
                        distribution = exact(algo = "shift"))
it@distribution@name
pvalue(it)

try(independence_test(y3 ~ x | b, data = dta,
                      distribution = exact(algo = "split-up")))

### two-sample without block
it <- independence_test(y3 ~ x, data = dta,
                        distribution = exact(algo = "auto"))
it@distribution@name
pvalue(it)

it <- independence_test(y3 ~ x, data = dta,
                        distribution = exact(algo = "shift"))
it@distribution@name
pvalue(it)

it <- independence_test(y3 ~ x, data = dta,
                        distribution = exact(algo = "split-up"))
it@distribution@name
pvalue(it)

### paired two-sample
it <- independence_test(y3 ~ x | b, data = dta2, paired = TRUE,
                        distribution = exact(algo = "auto"))
it@distribution@name
pvalue(it)

it <- independence_test(y3 ~ x | b, data = dta2, paired = TRUE,
                        distribution = exact(algo = "shift"))
it@distribution@name
pvalue(it)

try(independence_test(y3 ~ x | b, data = dta2, paired = TRUE,
                      distribution = exact(algo = "split-up")))


### check exact tests with weights
itw1 <- independence_test(y3 ~ x, data = dta, weights = ~ w,
                          distribution = exact(algorithm = "shift"))
itw2 <- independence_test(y3 ~ x, data = dta, weights = ~ w,
                          distribution = exact(algorithm = "split-up"))
y3w <- with(dta, rep(y3, w))
xw <- with(dta, rep(x, w))
it1 <- independence_test(y3w ~ xw, distribution = exact(algorithm = "shift"))
it2 <- independence_test(y3w ~ xw, distribution = exact(algorithm = "split-up"))
isequal(pvalue(itw1), pvalue(it1))
isequal(pvalue(itw1), pvalue(it2))
isequal(pvalue(itw2), pvalue(it1))
isequal(pvalue(itw2), pvalue(it2))

Convictions <-
    matrix(c(2, 10, 15, 3), nrow = 2,
           dimnames = list(c("Dizygotic", "Monozygotic"),
                           c("Convicted", "Not convicted")))
itw1 <- independence_test(as.table(Convictions), alternative = "less",
                          distribution = exact(algorithm = "shift"))
itw2 <- independence_test(as.table(Convictions), alternative = "less",
                          distribution = exact(algorithm = "split-up"))
it1 <- independence_test(Var2 ~ Var1, alternative = "less",
                         data = coin:::table2df(as.table(Convictions)),
                         distribution = exact(algorithm = "shift"))
it2 <- independence_test(Var2 ~ Var1, alternative = "less",
                         data = coin:::table2df(as.table(Convictions)),
                         distribution = exact(algorithm = "split-up"))
isequal(pvalue(itw1), pvalue(it1))
isequal(pvalue(itw1), pvalue(it2))
isequal(pvalue(itw2), pvalue(it1))
isequal(pvalue(itw2), pvalue(it2))
