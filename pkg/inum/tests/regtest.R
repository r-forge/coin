
library("inum")
data("iris")
set.seed(29)

iris[3, "Sepal.Width"] <- NA

iris1 <- inum(iris, nmax = 5, as.interval = "Sepal.Width")

iris1a <- inum(iris, nmax = 5, as.interval = c("Sepal.Width", "Sepal.Length"))

iris2 <- inum(iris, nmax = 5, total = TRUE, as.interval = "Sepal.Width")
iris2cc <- inum(iris, nmax = 5, total = TRUE, as.interval = "Sepal.Width", complete.cases.only = TRUE)

x1 <- as.data.frame(iris1)

table(x1$Species, iris$Species)

tapply(iris$Sepal.Width, x1$Sepal.Width, range)
levels(x1$Sepal.Width)

as.data.frame(iris2)
(w <- weights(iris2))
sum(w)

as.data.frame(iris2cc)
(w <- weights(iris2cc))
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

is.na(enum(c(NA, 1:3)))
is.na(interval(c(NA, runif(100))))

