
library(coin)
library(multcomp)

d <- data.frame(y = rnorm(100), x = gl(4, 25))
d$y[d$x %in% c(2,4)] <- d$y[d$x %in% c(2,4)] + 0.5

sp <- new("SymmetryProblem", d["x"], d["y"])
sp@y[[1]] <- sp@y[[1]] + as.numeric(sp@block)
d$b <- sp@block

cm <- function(x) {
    model.matrix(~ x - 1) %*% t(contrMat(table(x), type = "Tukey") )
}

pt <- perm_test(sp, xfun = function(data) tfun(data, factor_fun = cm))

pvalue(pt, adjust = TRUE)
pvalue(pt)
friedman.test(y ~ x | b, data = d)$p.value

