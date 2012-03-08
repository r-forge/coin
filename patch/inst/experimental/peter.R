
### example from your email
peter <- read.table("peter.csv", header = FALSE, sep = " ")[,-29]
names(peter)[29] <- "group"

library("coin")

fm <- paste(paste("V", 1:28, sep = "", collapse = " + "), " ~ group")
fm <- as.formula(fm)

### univariate p-value for V1
pvalue(independence_test(V1 ~ group, data = peter, 
                  distr = approximate(B = 100000)))

it <- independence_test(fm, data = peter, 
                        distr = approximate(B = 10000))

### global p-value
pvalue(it)

### permutation single-step
pvalue(it, adjusted = TRUE)

### permutation step-down
pvalue(it, adjusted = TRUE, method = "step-down")

### permutation Bonferroni
pvalue(it, adjusted = TRUE, method = "discrete")


### example from Westfall & Wolfinger (1997), Table 2 (two-sided)
df <- data.frame(group = factor(c(rep("Control", 50), rep("Treatment", 48))),
                 V1 = c(rep(0, 50), rep(1, 0), rep(0, 48 - 5), rep(1, 5)),
                 V2 = c(rep(0, 50 - 4), rep(1, 4), rep(0, 48 - 3), rep(1, 3)),
                 V3 = c(rep(0, 50), rep(1, 0), rep(0, 48 - 4), rep(1, 4)),
                 V4 = c(rep(0, 50 - 6), rep(1, 6), rep(0, 48 - 4), rep(1, 4)))

### check Fisher-test p-values: OK
sapply(df[,2:5], function(x) fisher.test(table(x, df$group), 
                                         alternative = "two")$p.value)

it <- independence_test(V1 + V2 + V3 + V4 ~ group, data = df, 
                        distr = approximate(B = 90000))

### page 7: 0.05261 for V1
pvalue(it, adjusted = TRUE, method = "discrete")

### example from Westfall & Wolfinger (1997), Table 4
tab <- as.table(matrix(c(12, 1, 3, 8, 17, 12, 9, 9, 16, 24), nrow = 2,
    dimnames = list(group = c("Placebo", "Active"), 
                    response = c("Very Poor", "Poor", "Fair", "Good", "Excellent"))))
df <- coin:::table2df(tab)

it <- independence_test(response ~ group, data = df, 
                        distr = approximate(B = 10000))

### Table 5, last column: OK
pvalue(it, adjusted = TRUE, method = "step-down")

### Table 5, 5th column: hm?
pvalue(it, adjusted = TRUE, method = "discrete")
