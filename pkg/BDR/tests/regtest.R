
library("BDR")
data("iris")

iris[3, "Sepal.Width"] <- NA

iris1 <- BDR(iris, nmax = 5, as.interval = "Sepal.Width")

iris2 <- BDR(iris, nmax = 5, total = TRUE, as.interval = "Sepal.Width")

all.equal(x1 <- as.data.frame(iris1), x2 <- as.data.frame(iris2))

table(x1$Species, iris$Species)

tapply(iris$Sepal.Width, x1$Sepal.Width, range)
levels(x1$Sepal.Width)

(w <- weights(iris2))
sum(w)

x <- runif(100)
x[1:3] <- NA   
ix <- interval(x, breaks = 0:10/10)

levels(ix)
nlevels(ix)
ix
  
table(ix)
ix[1:10]

enum(gl(3, 3))
enum(gl(3, 3, ordered = TRUE))
enum(c(TRUE, FALSE))
enum(c(1:3, 20L, 30L))

x <- sample(c(1:3, 10L, 20L), 100, replace = TRUE)
x[1:3] <- NA
ix <- enum(x)
levels(ix)   
nlevels(ix)  
ix
  
table(ix)

