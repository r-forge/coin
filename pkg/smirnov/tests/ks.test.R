
library("smirnov")

set.seed(29)
x <- runif(15) * 10
y <- runif(20) * 10
B <- 1e5

### two-sided asymptotic
stats::ks.test(x = x, y = y, exact = FALSE)
smirnov::ks.test(x = x, y = y, exact = FALSE)

### two-sided exact
stats::ks.test(x = x, y = y, exact = TRUE)
smirnov::ks.test(x = x, y = y, exact = TRUE)
smirnov::ks.test(x = x, y = y, exact = FALSE, simulate.p.value = TRUE, B = B)

### less asymptotic
stats::ks.test(x = x, y = y, exact = FALSE, alternative = "less")
smirnov::ks.test(x = x, y = y, exact = FALSE, alternative = "less")

### less exact
### NOTE: stats::ks.test reports asymptotic p-value SILENTLY
stats::ks.test(x = x, y = y, exact = TRUE, alternative = "less")
smirnov::ks.test(x = x, y = y, exact = TRUE, alternative = "less")
smirnov::ks.test(x = x, y = y, exact = FALSE, simulate.p.value = TRUE, 
                 alternative = "less", B = B)

### greater asymptotic
stats::ks.test(x = x, y = y, exact = FALSE, alternative = "greater")
smirnov::ks.test(x = x, y = y, exact = FALSE, alternative = "greater")

### greater exact
### NOTE: stats::ks.test reports asymptotic p-value SILENTLY
stats::ks.test(x = x, y = y, exact = TRUE, alternative = "greater")
smirnov::ks.test(x = x, y = y, exact = TRUE, alternative = "greater")
smirnov::ks.test(x = x, y = y, exact = FALSE, simulate.p.value = TRUE, 
                 B = B, alternative = "greater")

x <- round(x)
length(unique(x))
y <- round(y)
length(unique(y))

### two-sided asymptotic
stats::ks.test(x = x, y = y, exact = FALSE)
smirnov::ks.test(x = x, y = y, exact = FALSE)

### two-sided exact
stats::ks.test(x = x, y = y, exact = TRUE)
smirnov::ks.test(x = x, y = y, exact = TRUE)
smirnov::ks.test(x = x, y = y, exact = FALSE, simulate.p.value = TRUE, B = B)

### less asymptotic
stats::ks.test(x = x, y = y, exact = FALSE, alternative = "less")
smirnov::ks.test(x = x, y = y, exact = FALSE, alternative = "less")

### less exact
### NOTE: stats::ks.test reports asymptotic p-value SILENTLY
stats::ks.test(x = x, y = y, exact = TRUE, alternative = "less")
smirnov::ks.test(x = x, y = y, exact = TRUE, alternative = "less")
smirnov::ks.test(x = x, y = y, exact = FALSE, simulate.p.value = TRUE, 
                 B = B, alternative = "less")

### greater asymptotic
stats::ks.test(x = x, y = y, exact = FALSE, alternative = "greater")
smirnov::ks.test(x = x, y = y, exact = FALSE, alternative = "greater")

### greater exact
### NOTE: stats::ks.test reports asymptotic p-value SILENTLY
stats::ks.test(x = x, y = y, exact = TRUE, alternative = "greater")
smirnov::ks.test(x = x, y = y, exact = TRUE, alternative = "greater")
smirnov::ks.test(x = x, y = y, exact = FALSE, simulate.p.value = TRUE, 
                 B = B, alternative = "greater")

### an extreme example
x <- runif(5) * 10
y <- runif(500) * 10
### two-sided asymptotic
smirnov::ks.test(x = x, y = y, exact = FALSE)
smirnov::ks.test(x = x, y = y, exact = TRUE)
smirnov::ks.test(x = x, y = y, exact = FALSE, simulate.p.value = TRUE, B = B)

### even more extreme example
x <- 1:5
y <- round(runif(500) * 10)
### two-sided asymptotic
smirnov::ks.test(x = x, y = y, exact = FALSE)
smirnov::ks.test(x = x, y = y, exact = TRUE)
smirnov::ks.test(x = x, y = y, exact = FALSE, simulate.p.value = TRUE, B = B)
