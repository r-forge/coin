
library("FactoMineR")
library("coin")

example(coeffRV)

X <- as.matrix(X)
Y <- as.matrix(Y)

fm <- paste(paste(names(Y), collapse = "+"), "~", 
            paste(names(X), collapse = "+"), collapse = "")
fm <- as.formula(fm)

it <- independence_test(fm, data = wine, teststat = "quad")

    coefficientRV <- function(X, Y) {
        if (dim(X)[[1]] != dim(Y)[[1]]) 
            stop("no the same dimension for X and Y")
        if (dim(X)[[1]] == 1) {
            print("1 configuration RV is  NA")
            rv = NA
        }
        else {
            Y <- scale(Y, scale = FALSE)
            X <- scale(X, scale = FALSE)
            W1 <- X %*% t(X)
            W2 <- Y %*% t(Y)
            rv <- sum(diag(W1 %*% W2))/(sum(diag(W1 %*% W1)) * 
                sum(diag(W2 %*% W2)))^0.5
        }
        return(rv)
    }


p <- q <- 2
n <- 100

pval <- numeric(1000)
for (i in 1:length(pval)) {
    x <- as.data.frame(matrix(runif(n * p * q), nrow = n))
    fm <- paste(paste(names(x)[1:p], collapse = "+"), "~", 
                paste(names(x)[(p+1):(p + q)], collapse = "+") , collapse = "")
    fm <- as.formula(fm)
    pval[i] <- pvalue(independence_test(fm, data = x, teststat = "quad"))
}
