
library("wilcox")

x <- rnorm(50)
y <- rnorm(50)

.wilcox_score(x, y)

x <- rnorm(500)
y <- rnorm(500)

.wilcox_score(x, y, exact = FALSE)

x <- rnorm(5000)
y <- rnorm(5000)

.wilcox_score(x, y, exact = FALSE)

