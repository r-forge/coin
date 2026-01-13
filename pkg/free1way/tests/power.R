
library("free1way")
options(digits = 4)

args <- expand.grid(n = c(15, 50),
                    delta = c(0, 1),
                    alloc_ratio = c(1, 3))
args$seed <- 290875
args$pwr <- 0
p <- which(names(args) == "pwr")

for (i in 1:nrow(args))
    args$pwr[i] <- do.call("power.free1way.test", args[i,-p])$power

args
