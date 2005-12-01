
library("coin")

plin <- function(x, y, s)
    coin:::LinearStatistic(x, y[s, , drop = FALSE], rep(1, nrow(x)))

foo <- function(group, set1, set2, s, x, y, p) {
    i1 <- s[group == "1"]
    if (length(set1) > 0)
        i1 <- i1[i1 > max(set1)]
    i2 <- s[group == "2"]
    if (length(set2) > 0)
        i2 <- i2[i2 > max(set2)]

    sd <- s
    sd[set2] <- s[set1]
    sd[set1] <- s[set2]

    for (j in i2) {
        for (i in i1) {
            tmp <- sd
            tmp[i] <- sd[j]
            tmp[j] <- sd[i]
            T <<- rbind(T, c(plin(x, y, tmp), p))
            foo(group, c(set1, i), c(set2, j), s, x, y, p)
        }
    }
}

x <- gl(2, 5)
xx <- model.matrix(~ x - 1)
y <- 1:10
yy <- matrix(y)
indx <- xx
storage.mode(indx) <- "logical"
s <- 1:10
p <- prod(factorial(colSums(xx)))

T <<- c()
foo(x, set1 = numeric(0), set2 = numeric(0), s, xx, yy, p = p)
T <<- rbind(c(plin(xx, yy, s), p), T)

T <- T[order(T[,1]),]
sum(T[,3]) == factorial(10)


it <- wilcox_test(y ~ x, dist = exact())
support(it) * sqrt(variance(it)) + expectation(it)
