
library("inum")

### there was a warning; reported by Fabian Scheipl
x <- 1:2 + .1
inum(data.frame(x = x))


### by Susanne Dandl
sepallen <- iris[, "Sepal.Length", drop = FALSE]
sepallen$Sepal.Length[c(1, 10)] <- NA

a <- inum(sepallen, nmax = 5, as.interval = "Sepal.Length")
b <- inum(sepallen, nmax = 5, total = TRUE)
c <- inum(sepallen, nmax = 5, total = TRUE, complete.cases.only = TRUE)
all.equal(length(a), length(b), length(c))

cbind(sepallen, a, as.numeric(b), as.numeric(c))

stopifnot(length(attr(b, "levels")[unclass(b),"Sepal.Length"]) == 150)
stopifnot(length(attr(c, "levels")[unclass(c),"Sepal.Length"]) == 148)
