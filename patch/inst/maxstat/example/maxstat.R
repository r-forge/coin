
library("coin")
set.seed(290875)

load("preOP_maxstat.rda")

preOP <- subset(preOP, complete.cases(preOP))

## potential combinations of levels
lev1 <- c("pT0", "pT1", "pT2", "pT3", "pT4")
lev2 <- c("pN0", "pN1", "pN2 + pN3")
lev <- expand.grid(t = lev1, n = lev2)

## potential partitions
part <- coin:::fsplits(length(lev1) * length(lev2))
nrow(part) ## 2^(15-1) - 1 not 2^(13-1) - 1

## convenience function for displaying
## and checking partitions
show_part <- function(x) {
  table(lev[x,][,1], lev[x,][,2])
}

is_marg_ordered <- function(x) {
  x <- table(lev[x,][,1], lev[x,][,2])

  x1 <- apply(x, 1, function(z) sum(diff(z)))
  x2 <- apply(x, 2, function(z) sum(diff(z)))
  x1abs <- apply(x, 1, function(z) sum(abs(diff(z))))
  x2abs <- apply(x, 2, function(z) sum(abs(diff(z))))

  rval1 <- TRUE ## all(x1 >= 0) | all(x1 <= 0)
  rval2 <- TRUE ## all(x2 >= 0) | all(x2 <= 0)
  rval3 <- all(c(x1abs, x2abs) < 2)

  all(c(rval1, rval2, rval3))
}

## example for ordered and non-ordered partitions
show_part(part[31,])
is_marg_ordered(part[31,])

show_part(part[32,])
is_marg_ordered(part[32,])

## obtain ordered partitions
wi <- sapply(1:nrow(part), function(i) is_marg_ordered(part[i,]))
part2 <- part[wi,]

tn <- preOP$tclass:preOP$nclass
colnames(part2) <- gsub(":", ".", levels(tn))
mymaxstat_trafo <- function(x, minprob = 0.1, maxprob = 0.9) {

    sp <- part2 ###fsplits(nlevels(x))
    x <- x[[1]]
    lev <- levels(x)
    sp <- sp[,colnames(sp) %in% lev]
    tr <- matrix(0, nrow = length(x), ncol = nrow(sp))
    cn <- vector(mode = "character", length = nrow(sp))
    for (i in 1:nrow(sp)) {
        tr[ ,i] <- x %in% lev[sp[i, ]]
        cn[i] <- paste("{", paste(lev[sp[i, ]], collapse = ", "), "} vs. {",
                       paste(lev[!sp[i, ]], collapse = ", "), "}", sep = "")
    }
    rownames(tr) <- 1:length(x)
    colnames(tr) <- cn
    tr <- tr[, colMeans(tr) >= minprob & colMeans(tr) <= maxprob]
    tr
}


mt <- independence_test(Surv(time, event) ~ tn, data = preOP,
    xtrafo = mymaxstat_trafo, distr = approximate(B = 50000))

teststat <- statistic(mt)
stat <- statistic(mt, "standardized")
cutpoint <- mt@statistic@xtrans[,which.max(abs(stat))]
pval <- pvalue(mt)

save(preOP, teststat, cutpoint, pval, stat, file = "maxstat.rda")

risk <- as.factor(cutpoint) ###preOP$tn %in% cutpoint
table(risk)
pvalue(mt)

table(risk, preOP$stadium)

library("survival")
plot(survfit(Surv(time, event) ~ risk, data = preOP))

