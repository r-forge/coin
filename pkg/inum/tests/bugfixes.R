
library("inum")
set.seed(29)

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

### by Susanne Dandl
## mini data frame with some missings
d <- data.frame(
  y = rep(1:5, each = 2),
  x = factor(rep(0:1, 5), labels = c("a", "b")),
  z = 1:10,
  w = 0:9/9
)
d$y[c(1, 10)] <- NA

i <- inum(d, total = TRUE, complete = FALSE)
attr(i, "levels")[i,]

i <- inum(d, total = TRUE, complete = TRUE)
rbind(NA, attr(i, "levels"))[i + 1,]

d <- expand.grid(y = 1:5, z = 1:10)
d$y[c(1, nrow(d))] <- NA
d$w <- rpois(nrow(d), lambda = 3)

i1 <- inum(d, total = TRUE, complete = FALSE)
attr(i1, "levels")[i1,]

i2 <- inum(d, total = TRUE, complete = TRUE)
rbind(NA, attr(i2, "levels"))[i2 + 1,]
