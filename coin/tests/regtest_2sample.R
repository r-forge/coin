
PKG <- "coin"

### Regression tests for the 2 sample problem, i.e.,
### testing the independence of a numeric variable
### `y' and a binary factor `x' (possibly blocked)

set.seed(290875)
library(PKG, character.only = TRUE)

### generate data
dat <- data.frame(x = gl(2, 50), y = rnorm(100), block = gl(5, 20))[sample(1:100,
75),]

### Wilcoxon Mann-Whitney Rank Sum Test

### asymptotic distribution
ptwo <- wilcox.test(y ~ x, data = dat, correct = FALSE, exact = FALSE)$p.value
pless <- wilcox.test(y ~ x, data = dat, alternative = "less", 
                     correct = FALSE, exact = FALSE)$p.value
pgreater <- wilcox.test(y ~ x, data = dat, alternative = "greater", 
                        correct = FALSE, exact = FALSE)$p.value

stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "less")),
                  pless))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater")), 
                  pgreater))

stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank))), ptwo))
stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
    alternative = "less")), pless))
stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
    alternative = "greater")), pgreater))

### exact distribution
ptwo <- wilcox.test(y ~ x, data = dat, exact = TRUE)$p.value
pless <- wilcox.test(y ~ x, data = dat, alternative = "less", exact = TRUE)$p.value
pgreater <- wilcox.test(y ~ x, data = dat, alternative = "greater", 
                        exact = TRUE)$p.value

stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, distribution = "exact")), 
                  ptwo))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "less", 
                             distribution = "exact")), pless))
stopifnot(isequal(pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater",
                             distribution = "exact")), pgreater))

stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank))), ptwo))
stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
    alternative = "less")), pless))
stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = rank), 
    alternative = "greater")), pgreater))

### approximated distribution
rtwo <- pvalue(wilcox_test(y ~ x, data = dat, distribution = "approx")) / ptwo
rless <- pvalue(wilcox_test(y ~ x, data = dat, alternative = "less",
                   distribution = "approx")) / pless
rgreater <- pvalue(wilcox_test(y ~ x, data = dat, alternative = "greater",
                   distribution = "approx")) / pgreater
stopifnot(all(c(rtwo, rless, rgreater) > 0.9 & 
              c(rtwo, rless, rgreater) < 1.1))

### <FIXME> add block examples </FIXME>

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact"))

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt", 
                   alternative = "less"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx", 
                   alternative = "less"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact", 
                   alternative = "less"))

pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "asympt", 
                   alternative = "greater"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "approx", 
                   alternative = "greater"))
pvalue(wilcox_test(y ~ x | block, data = dat, distribution = "exact", 
                   alternative = "greater"))

### StatXact 6 manual, page 345-346
load("employment.rda")
stopifnot(isequal(round(pvalue(wilcox_test(Salary ~ Gender | Year, data = employment, 
                         distribution = "exact")), 4), 0.04))
stopifnot(isequal(round(pvalue(wilcox_test(Salary ~ Gender | Year, data = employment, 
                         distribution = "exact", alternative = "less")), 4), 0.04))
stopifnot(isequal(round(pvalue(wilcox_test(Salary ~ Gender | Year, data = employment, 
                         distribution = "exact", alternative = "greater")), 4), 1))
stopifnot(isequal(round(pvalue(wilcox_test(Salary ~ Gender | Year, data = employment, 
                         distribution = "asympt")), 4), 0.0578))
stopifnot(isequal(round(pvalue(wilcox_test(Salary ~ Gender | Year, data = employment, 
                         distribution = "asympt", alternative = "less")), 4), 0.0289))

### sanify checks
try(wilcox_test(x ~ y, data = dat))
try(wilcox_test(x ~ y | y, data = dat))

### Ansari-Bradley Test 

### asymptotic distribution
### <FIXME> alternative is defined in another way here </FIXME>
ptwo <- ansari.test(y ~ x, data = dat, correct = FALSE, exact = FALSE)$p.value
pless <- ansari.test(y ~ x, data = dat, alternative = "greater", 
                     correct = FALSE, exact = FALSE)$p.value
pgreater <- ansari.test(y ~ x, data = dat, alternative = "less", 
                        correct = FALSE, exact = FALSE)$p.value

stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "less")),
                  pless))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "greater")), 
                  pgreater))

stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo))), ptwo))
stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
    alternative = "less")), pless))
stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "asympt", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
    alternative = "greater")), pgreater))

### exact distribution
ptwo <- ansari.test(y ~ x, data = dat, exact = TRUE)$p.value
pless <- ansari.test(y ~ x, data = dat, alternative = "greater", exact = TRUE)$p.value
pgreater <- ansari.test(y ~ x, data = dat, alternative = "less", 
                        exact = TRUE)$p.value

### <FIXME>: Definition of two-sided P-values! </FIXME>
(isequal(pvalue(ansari_test(y ~ x, data = dat, distribution = "exact")), 
                  ptwo))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "less", 
                             distribution = "exact")), pless))
stopifnot(isequal(pvalue(ansari_test(y ~ x, data = dat, alternative = "greater",
                             distribution = "exact")), pgreater))

### <FIXME>: Definition of two-sided P-values! </FIXME>
(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo))), ptwo))
stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
    alternative = "less")), pless))
stopifnot(isequal(pvalue(perm_test(y ~ x, data = dat, distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = ansari_trafo), 
    alternative = "greater")), pgreater))

### approximated distribution
rtwo <- pvalue(ansari_test(y ~ x, data = dat, distribution = "approx")) / ptwo
rless <- pvalue(ansari_test(y ~ x, data = dat, alternative = "less",
                   distribution = "approx")) / pless
rgreater <- pvalue(ansari_test(y ~ x, data = dat, alternative = "greater",
                   distribution = "approx")) / pgreater
### <FIXME> ??? </FIXME>
(all(c(rtwo, rless, rgreater) > 0.9 & 
              c(rtwo, rless, rgreater) < 1.1))

### <FIXME> add block examples </FIXME>

### sanify checks
try(ansari_test(x ~ y, data = dat))
try(ansari_test(x ~ y | y, data = dat))

### the remaining three candidates
### asymptotic distribution
ptwo <- wilcox.test(y ~ x, data = dat, correct = FALSE)$p.value
pless <- wilcox.test(y ~ x, data = dat, alternative = "less",
                     correct = FALSE)$p.value
pgreater <- wilcox.test(y ~ x, data = dat, alternative = "greater", 
                        correct = FALSE)$p.value

perm_test(y ~ x, dat = dat)
perm_test(y ~ x, dat = dat, alternative = "less")
perm_test(y ~ x, dat = dat, alternative = "greater")

normal_test(y ~ x, dat = dat)
normal_test(y ~ x, dat = dat, alternative = "less")
normal_test(y ~ x, dat = dat, alternative = "greater")

median_test(y ~ x, dat = dat)
median_test(y ~ x, dat = dat, alternative = "less")
median_test(y ~ x, dat = dat, alternative = "greater")

