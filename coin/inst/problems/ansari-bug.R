
set.seed(290875)
library(exactRankTests)

### generate data
dat <- data.frame(x = gl(2, 50), y = rnorm(100), block = gl(5,
20))[sample(1:100, 75),]

ansari.test(y ~ x, data = dat, exact = TRUE)$p.value
ansari.exact(y ~ x, dat = dat, exact = TRUE)$p.value

ansari.test(y ~ x, data = dat, exact = TRUE)$p.value
ansari.exact(y ~ x, dat = dat, exact = TRUE)$p.value
a <- ansari.test(y ~ x, data = dat, exact = TRUE)

.C("pansari", as.integer(1), as.double(a$statistic), 
    as.integer(38), as.integer(37), PACKAGE = "stats")
