
library("coin")
set.seed(290875)

load("preOP_maxstat.rda")

mt <- maxstat_test(Surv(time, event) ~ tn, data = preOP, 
    distribution = approximate(B = 1000))

save(mt, file = "maxstat.rda")

risk <- preOP$tn %in% mt@statistic@estimates$estimate$cutpoint
table(risk)
pvalue(mt)

table(risk, preOP$stadium)

library("survival")
plot(survfit(Surv(time, event) ~ risk, data = preOP))

