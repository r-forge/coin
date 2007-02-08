
library("coin")
set.seed(290875)

load("preOP_maxstat.rda")

mt <- maxstat_test(Surv(time, event) ~ tn, data = preOP, 
    distribution = approximate(B = 1000))



teststat <- statistic(mt)
stat <- statistic(mt, "standardized")
cutpoint <- mt@statistic@estimates$estimate$cutpoint
pval <- pvalue(mt)

save(teststat, cutpoint, pval, stat, file = "maxstat.rda")

risk <- preOP$tn %in% cutpoint
table(risk)
pvalue(mt)

table(risk, preOP$stadium)

library("survival")
plot(survfit(Surv(time, event) ~ risk, data = preOP))

