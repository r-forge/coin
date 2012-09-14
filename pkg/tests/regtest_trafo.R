
### Regression tests for the transformation functions

set.seed(290875)
library(coin)
isequal <- coin:::isequal


### NA handling: continuous
x <- c(1, 2, NA, 3, 3, NA, 4, 5, NA)
cc <- complete.cases(x)

id_trafo(x)
id_trafo(x[cc])

rank_trafo(x)
rank_trafo(x[cc])
rank_trafo(x, ties.method = "random")
rank_trafo(x[cc], ties.method = "random")

ansari_trafo(x)
ansari_trafo(x[cc])
ansari_trafo(x, ties.method = "average-scores")
ansari_trafo(x[cc], ties.method = "average-scores")

fligner_trafo(x)
fligner_trafo(x[cc])
fligner_trafo(x, ties.method = "average-scores")
fligner_trafo(x[cc], ties.method = "average-scores")

normal_trafo(x)
normal_trafo(x[cc])
normal_trafo(x, ties.method = "average-scores")
normal_trafo(x[cc], ties.method = "average-scores")

median_trafo(x)
median_trafo(x[cc])

consal_trafo(x)
consal_trafo(x[cc])
consal_trafo(x, ties.method = "average-scores")
consal_trafo(x[cc], ties.method = "average-scores")

maxstat_trafo(x)
maxstat_trafo(x[cc])
maxstat_trafo(x, maxprob = 0.3)
maxstat_trafo(x[cc], maxprob = 0.3)


### NA handling: survival
x <- c(1, 2, NA, 3, 3, NA, 4, 5, NA)
cc <- complete.cases(x)

logrank_trafo(Surv(x))
logrank_trafo(Surv(x[cc]))
logrank_trafo(Surv(x), ties.method = "HL")
logrank_trafo(Surv(x[cc]), ties.method = "HL")
logrank_trafo(Surv(x), ties.method = "average-scores")
logrank_trafo(Surv(x[cc]), ties.method = "average-scores")

x <- c(1, 2, 3, 3, 3, 4, 4, 5, 5)
e <- rep(c(0, NA, 1, 1), length.out = 9)
cc <- complete.cases(x, e)

logrank_trafo(Surv(x, e))
logrank_trafo(Surv(x[cc], e[cc]))
logrank_trafo(Surv(x, e), ties.method = "HL")
logrank_trafo(Surv(x[cc], e[cc]), ties.method = "HL")
logrank_trafo(Surv(x, e), ties.method = "average-scores")
logrank_trafo(Surv(x[cc], e[cc]), ties.method = "average-scores")

x <- c(1, 2, NA, 3, 3, NA, 4, 5, NA)
e <- rep(c(0, NA, 1, 1), length.out = 9)
cc <- complete.cases(x, e)

logrank_trafo(Surv(x, e))
logrank_trafo(Surv(x[cc], e[cc]))
logrank_trafo(Surv(x, e), ties.method = "HL")
logrank_trafo(Surv(x[cc], e[cc]), ties.method = "HL")
logrank_trafo(Surv(x, e), ties.method = "average-scores")
logrank_trafo(Surv(x[cc], e[cc]), ties.method = "average-scores")


### NA handling: factor
x <- factor(c(1, 1, NA, 2, NA, 3, 3, NA, 4))
cc <- complete.cases(x)

f_trafo(x)
f_trafo(x[cc])

of_trafo(x)
of_trafo(x[cc])

fmaxstat_trafo(x)
fmaxstat_trafo(x[cc])
fmaxstat_trafo(x, maxprob = 0.4)
fmaxstat_trafo(x[cc], maxprob = 0.4)

mcp_trafo(x = "Tukey")(data.frame(x))
mcp_trafo(x = "Tukey")(data.frame(x = x[cc]))

x[9] <- NA
cc <- complete.cases(x)

f_trafo(x)
f_trafo(x[cc])

of_trafo(x)
of_trafo(x[cc])

fmaxstat_trafo(x)
fmaxstat_trafo(x[cc])
fmaxstat_trafo(x, maxprob = 0.4)
fmaxstat_trafo(x[cc], maxprob = 0.4)

mcp_trafo(x = "Tukey")(data.frame(x))
mcp_trafo(x = "Tukey")(data.frame(x = x[cc]))
