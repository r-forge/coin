
### Regression tests for multiple adjustments

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### example from Westfall & Wolfinger (1997), Table 4
tab <- as.table(matrix(c(12, 1, 3, 8, 17, 12, 9, 9, 16, 24), nrow = 2,
    dimnames = list(group = c("Placebo", "Active"),
                    response = c("Very Poor", "Poor", "Fair", "Good",
"Excellent"))))
df <- coin:::table2df(tab)

it <- independence_test(response ~ group, data = df,
                        distr = approximate(B = 100000))

### Table 5, last column: OK
pvalue(it, adjusted = TRUE, method = "step-down")
