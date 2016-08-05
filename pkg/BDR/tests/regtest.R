
library("BDR")
library("survival")
data("iris")

iris[3, "Sepal.Width"] <- NA

iris1 <- BDR(iris, nmax = 5)

iris2 <- BDR(iris, nmax = 5, total = TRUE)

all.equal(as.data.frame(iris1), as.data.frame(iris1))

iris$y <- Surv(runif(nrow(iris)), sample(c(TRUE, FALSE), nrow(iris), 
                                         replace = TRUE))

iris1 <- BDR(iris, nmax = 5)

iris2 <- BDR(iris, nmax = 5, total = TRUE)

all.equal(as.data.frame(iris1), as.data.frame(iris1))



