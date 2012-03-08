
library("coin")
### get rid of namespace
load(file.path(.find.package("coin"), "R", "all.rda"))
source("bootstrap.R")

library("multcomp")

     YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44, 
                                  42, 60, 32, 42, 45, 58, 27, 51, 42, 52, 
                                  38, 33, 26, 25, 28, 28, 26, 27, 27, 27, 
                                  31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                       site = factor(c(rep("I", 10), rep("II", 10),
                                       rep("III", 10), rep("IV", 10))))


NDWD <- oneway_test(length ~ site, data = YOY,
#             ytrafo = function(data) trafo(data, numeric_trafo = rank),
             xtrafo = function(data) trafo(data, factor_trafo = function(x)
                 model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
             teststat = "max")

mt.bootstrap(NDWD, B = 10000, xcontrast = FALSE)
mt.bootstrap(NDWD, B = 10000, xcontrast = TRUE)

mtest <- function(mu = 0.5) {

    y <- rnorm(100)
    x <- gl(4, 25)
    y[1:25] <- y[1:25] + mu
    boxplot(y ~ x)
    mydf <- data.frame(y = y, x = x)
    it <- oneway_test(y ~ x, data = mydf,
             xtrafo = function(data) trafo(data, factor_trafo = function(x)
                 model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
             teststat = "max")
    # print(drop(pvalue(it, "single")))
    mt.bootstrap(it, B = 1000, xcontrast = TRUE)
}

n <- 1000
p <- matrix(0, nrow = n, ncol = 6)
for (i in 1:nrow(p)) {
    print(i)
    p[i,] <- mtest(mu = 0.5)
}

colMeans(p < 0.05)
mean(apply(p[,4:6], 1, min) < 0.05)

