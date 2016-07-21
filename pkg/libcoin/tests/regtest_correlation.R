
### Regression tests for the correlation problem, i.e.,
### testing the independence of two numeric variables
### `x' and `y' (possibly blocked)

set.seed(290875)
library("coin")
library("libcoin")
source("check_vs_coin.R")

isequal <- coin:::isequal
options(useFancyQuotes = FALSE)

### generate data
dat <- data.frame(x = rnorm(100), y = rnorm(100), block = gl(10, 10))

### not really the same, T = (rank(x) - rank(y))^2 is used here
cor.test(~ x + y, data = dat, method = "spearman")$p.value
cor.test(~ x + y, data = dat, alternative = "less", method = "spearman")$p.value
cor.test(~ x + y, data = dat, alternative = "greater", method = "spearman")$p.value

### without blocks
pvalue(lc("spearman_test",y ~ x, data = dat))
pvalue(lc("spearman_test",x ~ y, data = dat))
pvalue(lc("spearman_test", ~ y + x, data = dat))
pvalue(lc("spearman_test", ~ x + y, data = dat))

pvalue(lc("fisyat_test",y ~ x, data = dat))
pvalue(lc("fisyat_test",x ~ y, data = dat))
pvalue(lc("fisyat_test", ~ y + x, data = dat))
pvalue(lc("fisyat_test", ~ x + y, data = dat))

pvalue(lc("quadrant_test",y ~ x, data = dat))
pvalue(lc("quadrant_test",x ~ y, data = dat))
pvalue(lc("quadrant_test", ~ y + x, data = dat))
pvalue(lc("quadrant_test", ~ x + y, data = dat))

pvalue(lc("koziol_test",y ~ x, data = dat))
pvalue(lc("koziol_test",x ~ y, data = dat))
pvalue(lc("koziol_test", ~ y + x, data = dat))
pvalue(lc("koziol_test", ~ x + y, data = dat))

### with blocks
pvalue(lc("spearman_test",y ~ x | block, data = dat))
pvalue(lc("spearman_test",x ~ y | block, data = dat))
pvalue(lc("spearman_test", ~ y + x | block, data = dat))
pvalue(lc("spearman_test", ~ x + y | block, data = dat))

pvalue(lc("fisyat_test",y ~ x | block, data = dat))
pvalue(lc("fisyat_test",x ~ y | block, data = dat))
pvalue(lc("fisyat_test", ~ y + x | block, data = dat))
pvalue(lc("fisyat_test", ~ x + y | block, data = dat))

pvalue(lc("quadrant_test",y ~ x | block, data = dat))
pvalue(lc("quadrant_test",x ~ y | block, data = dat))
pvalue(lc("quadrant_test", ~ y + x | block, data = dat))
pvalue(lc("quadrant_test", ~ x + y | block, data = dat))

pvalue(lc("koziol_test",y ~ x | block, data = dat))
pvalue(lc("koziol_test",x ~ y | block, data = dat))
pvalue(lc("koziol_test", ~ y + x | block, data = dat))
pvalue(lc("koziol_test", ~ x + y | block, data = dat))

### sanity checks, those should be errors
dat <- data.frame(x = gl(2, 50), y = rnorm(100), block = rnorm(100))

try(lc("spearman_test",y ~ x, data = dat))
try(lc("spearman_test",y ~ x | block, data = dat))

try(lc("fisyat_test",y ~ x, data = dat))
try(lc("fisyat_test",y ~ x | block, data = dat))

try(lc("quadrant_test",y ~ x, data = dat))
try(lc("quadrant_test",y ~ x | block, data = dat))

try(lc("koziol_test",y ~ x, data = dat))
try(lc("koziol_test",y ~ x | block, data = dat))
