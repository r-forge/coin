
library("libcoin")

w <- 1:nrow(iris)

y <- model.matrix(~ Species - 1, data = iris)
storage.mode(y) <- "double"

.Call("R_LinstatExpCov", iris, colnames(iris) != "Species", y, w)


