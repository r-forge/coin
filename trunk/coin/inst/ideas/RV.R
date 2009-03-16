
library("FactoMineR")
library("coin")

example(coeffRV)

X <- as.matrix(X)
Y <- as.matrix(Y)

fm <- paste(paste(colnames(Y), collapse = "+"), "~", 
            paste(colnames(X), collapse = "+"), collapse = "")
fm <- as.formula(fm)

it <- independence_test(fm, data = wine, teststat = "quad")

x <- it@statistic@xtrans
y <- it@statistic@ytrans
w <- it@statistic@weights
b <- as.integer(it@statistic@weights)

mc <- coin:::MCfun(x, y, w, b, B = 20000)
s1 <- apply(mc, 2, function(x) crossprod(x))

mean(s1 > sum(statistic(it, "linear")^2))
