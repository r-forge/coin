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
                        distribution = exact(algo = "shift", fact = 1e5))
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
                        distribution = exact(algo = "shift", fact = 1e5))
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
stopifnot(isequal(pvalue(itw1), pvalue(it1)))
stopifnot(isequal(pvalue(itw1), pvalue(it2)))
stopifnot(isequal(pvalue(itw2), pvalue(it1)))
stopifnot(isequal(pvalue(itw2), pvalue(it2)))

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
stopifnot(isequal(pvalue(itw1), pvalue(it1)))
stopifnot(isequal(pvalue(itw1), pvalue(it2)))
stopifnot(isequal(pvalue(itw2), pvalue(it1)))
stopifnot(isequal(pvalue(itw2), pvalue(it2)))


### check support, pperm, dperm, qperm,
y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
dta1 <- data.frame(
    y = c(y1, y2),
    x = gl(2, length(y1)),
    b = factor(rep(seq_along(y1), 2)))

diff <- y1 - y2
dta2 <- data.frame(
    y = as.vector(rbind(abs(diff) * (diff >= 0), abs(diff) * (diff < 0))),
    x = factor(rep(0:1, length(diff)), labels = c("pos", "neg")),
    b <- gl(length(diff), 2))

### shift without block
it1_SR <- independence_test(y ~ x, data = dta1,
                            distribution = exact(algorithm = "shift"))
supp_it1_SR <- support(it1_SR)
stopifnot(!is.unsorted(supp_it1_SR))
stopifnot(all(supp_it1_SR == unique(supp_it1_SR)))
round(pp_it1_SR <- pperm(it1_SR, supp_it1_SR), 7)
round(dp_it1_SR <- dperm(it1_SR, supp_it1_SR), 7)
round(qp_it1_SR <- qperm(it1_SR, seq(0, 1, 0.01)), 7)

### split-up without block
it1_vdW <- independence_test(y ~ x, data = dta1,
                             distribution = exact(algorithm = "split-up"))
round(pp_it1_vdW <- pperm(it1_vdW, supp_it1_SR), 7)
round(qp_it1_vdW <- qperm(it1_vdW, seq(0, 1, 0.01)), 7)

### should be equal
stopifnot(isequal(pp_it1_SR, pp_it1_vdW))
stopifnot(isequal(qp_it1_SR[-c(1, 101)], qp_it1_vdW[-c(1, 101)]))
stopifnot(isequal(pvalue(it1_SR), pvalue(it1_vdW)))

### shift with block
it2_SR <- independence_test(y ~ x | b, data = dta1,
                            distribution = exact(algorithm = "shift"))
supp_it2_SR <- support(it2_SR)
stopifnot(!is.unsorted(supp_it2_SR))                         # failed in < 1.1-0
stopifnot(all(supp_it2_SR == unique(supp_it2_SR)))           # failed in < 1.1-0

round(pp_it2_SR <- pperm(it2_SR, supp_it2_SR), 7)
round(dp_it2_SR <- dperm(it2_SR, supp_it2_SR), 7)            # failed in < 1.1-0
round(qp_it2_SR <- qperm(it2_SR, seq(0, 1, 0.01)), 7)

### paired shift with block
it3_SR <- independence_test(y ~ x | b, data = dta2,
                            distribution = exact(algorithm = "shift"),
                            paired = TRUE)
supp_it3_SR <- support(it3_SR)
stopifnot(!is.unsorted(supp_it3_SR))
stopifnot(all(supp_it3_SR == unique(supp_it3_SR)))
round(pp_it3_SR <- pperm(it3_SR, supp_it3_SR), 7)
round(dp_it3_SR <- dperm(it3_SR, supp_it3_SR), 7)
round(qp_it3_SR <- qperm(it3_SR, seq(0, 1, 0.01)), 7)

### should be equal
stopifnot(isequal(pp_it2_SR, pp_it3_SR))                     # failed in < 1.1-0
stopifnot(isequal(qp_it2_SR, qp_it3_SR))                     # failed in < 1.1-0
stopifnot(isequal(pvalue(it2_SR), pvalue(it3_SR)))


### exact test based on quadratic forms

### shift with block
itq1 <- independence_test(y ~ x | b, data = dta1,
                          distribution = exact(algorithm = "shift"),
                          teststat = "quad")
its1 <- independence_test(y ~ x | b, data = dta1,
                          distribution = exact(algorithm = "shift"),
                          teststat = "scalar")
stopifnot(isequal(statistic(itq1), statistic(its1)^2))
stopifnot(isequal(pvalue(itq1), pvalue(its1)))
stopifnot(isequal(support(itq1), support(its1)[support(its1) >= 0]^2))

### paired shift with block
its2 <- independence_test(y ~ x | b, data = dta2,
                          distribution = exact(algorithm = "shift"),
                          paired = TRUE, teststat = "scalar")
itq2 <- independence_test(y ~ x | b, data = dta2,
                          distribution = exact(algorithm = "shift"),
                          paired = TRUE, teststat = "quad")
stopifnot(isequal(statistic(itq2), statistic(its2)^2))
stopifnot(isequal(pvalue(itq2), pvalue(its2)))
stopifnot(isequal(support(itq2), support(its2)[support(its2) >= 0]^2))

### should be equal
stopifnot(isequal(statistic(itq1), statistic(itq2)))
stopifnot(isequal(pvalue(itq1), pvalue(itq2)))
stopifnot(isequal(support(itq1), support(itq2)))
