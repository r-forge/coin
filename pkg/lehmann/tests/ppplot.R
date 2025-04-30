
library("lehmann")

N <- 50
x <- gl(2, N, labels = LETTERS[1:2])
y <- rlogis(length(x), location = c(0, .5)[x])
y <- rnorm(length(x), sd = c(1, 2)[x])

ppplot(y ~ x, pch = 19)
ppplot(y ~ x, pch = 19, conf.level = 0.95, 
       conf.args = list(type = "Wilcoxon", col = "lightgrey"))

