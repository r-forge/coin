
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

### Table 5, first column: OK
pvalue(it, method = "unadjusted")

### Table 5, next-to-last column: OK
pvalue(it, method = "Sidak-Holm")

### Table 5, last column: OK
pvalue(it, method = "step-down")

### example from Westfall & Wolfinger (1997), Table 1
df <- data.frame(group = factor(c(rep("Control", 50), rep("Treatment", 48))),
                 V1 = c(rep(0, 50), rep(1, 0), rep(0, 48 - 5), rep(1, 5)),
                 V2 = c(rep(0, 50 - 4), rep(1, 4), rep(0, 48 - 3), rep(1, 3)),
                 V3 = c(rep(0, 50), rep(1, 0), rep(0, 48 - 4), rep(1, 4)),
                 V4 = c(rep(0, 50 - 6), rep(1, 6), rep(0, 48 - 4), rep(1, 4)))

### alternative: less
it <- independence_test(V1 + V2 + V3 + V4 ~ group, data = df,
                        distr = approximate(B = 100000), alt = "less")

### page 4, 2nd column: adjusted p-value = 0.03665 for V1
pvalue(it, method = "Sidak")

### page 4, 2nd column: adjusted p-value = 0.03698 for V1
### Note: 0.02521 + 0.00532 + 0 + 0.00645 = 0.03698
pvalue(it, method = "Bonferroni")

### alternative: two
it <- independence_test(V1 + V2 + V3 + V4 ~ group, data = df,
                        distr = approximate(B = 100000), alt = "two")

### page 5, 1st column: adjusted p-value = 0.05261 for V1
pvalue(it, method = "Sidak")

### page 5, 1st column: adjusted p-value = 0.05352 for V1
### Note: 0.02521 + 0.01254 + 0 + 0.01577 = 0.05352
pvalue(it, method = "Bonferroni")

### artificial example, checked against `multtest:mt.maxT'

set.seed(290875)

gr <- gl(2, 50) 
x1 <- rnorm(100) + (as.numeric(gr) - 1) * 0.5
x2 <- rnorm(100) - (as.numeric(gr) - 1) * 0.5

pvalue(independence_test(x1 + x2 ~ gr, alt = "two.sided"), "single-step")
pvalue(independence_test(x1 + x2 ~ gr, alt = "less"), "single-step")
pvalue(independence_test(x1 + x2 ~ gr, alt = "greater"), "single-step")

pvalue(independence_test(x1 + x2 ~ gr, alt = "two.sided"), "step-down")
pvalue(independence_test(x1 + x2 ~ gr, alt = "less"), "step-down")
pvalue(independence_test(x1 + x2 ~ gr, alt = "greater"), "step-down")

pvalue(independence_test(x1 + x2 ~ gr, alt = "two.sided", 
                         dist = approximate(B = 10000)), "single-step")
pvalue(independence_test(x1 + x2 ~ gr, alt = "less", 
                         dist = approximate(B = 10000)), "single-step")
pvalue(independence_test(x1 + x2 ~ gr, alt = "greater", 
                         dist = approximate(B = 10000)), "single-step")

pvalue(independence_test(x1 + x2 ~ gr, alt = "two.sided", 
                         dist = approximate(B = 10000)), "step-down")
pvalue(independence_test(x1 + x2 ~ gr, alt = "less", 
                         dist = approximate(B = 10000)), "step-down")
pvalue(independence_test(x1 + x2 ~ gr, alt = "greater", 
                         dist = approximate(B = 10000)), "step-down")

if (FALSE) {
    #library("multtest")
    #a <- mt.maxT(t(cbind(x1, x2)), as.numeric(gr) - 1)
    #a[order(a$index),]
    #a <- mt.maxT(t(cbind(x1, x2)), as.numeric(gr) - 1, side = "upper")
    #a[order(a$index),]
    #a <- mt.maxT(t(cbind(x1, x2)), as.numeric(gr) - 1, side = "lower")
    #a[order(a$index),]
}

### Monte-Carlo distribution

y <- rnorm(20)
x <- runif(20)

mt <- maxstat_test(y ~ x, distribution = approximate())
pvalue(mt)
pperm(mt, 1)
qperm(mt, 0.9)
dperm(mt, qperm(mt, 0.9))
support(mt)

mt <- maxstat_test(y ~ x, distribution = approximate(), alternative = "greater")
pvalue(mt)
pperm(mt, 1)
qperm(mt, 0.9)
dperm(mt, qperm(mt, 0.9))
support(mt)

mt <- maxstat_test(y ~ x, distribution = approximate(), alternative = "less")
pvalue(mt)
pperm(mt, 1)
qperm(mt, 0.9)
dperm(mt, qperm(mt, 0.9))
support(mt)

### unadjusted

set.seed(290875)

gr <- gl(3, 50) 
x1 <- rnorm(150) + (as.numeric(gr) - 1) * 0.5
x2 <- rnorm(150) - (as.numeric(gr) - 1) * 0.5

pvalue(it1 <- independence_test(x1 + x2 ~ gr, alt = "two.sided"), "unadjusted")
pvalue(it2 <- independence_test(x1 + x2 ~ gr, alt = "less"), "unadjusted")
pvalue(it3 <- independence_test(x1 + x2 ~ gr, alt = "greater"), "unadjusted")

pvalue(it4 <- independence_test(x1 + x2 ~ gr, alt = "two.sided", 
                                dist = approximate(B = 10000)), "unadjusted")
pvalue(it5 <- independence_test(x1 + x2 ~ gr, alt = "less", 
                                dist = approximate(B = 10000)), "unadjusted")
pvalue(it6 <- independence_test(x1 + x2 ~ gr, alt = "greater", 
                                dist = approximate(B = 10000)), "unadjusted")

### consistency of minimum p-value for "global"/"single-step"/"step-down"

set.seed(290875); pg1 <- pvalue(it1, "global")[1]
set.seed(290875); pss1 <- pvalue(it1, "single-step")
set.seed(290875); psd1 <- pvalue(it1, "step-down")
identical(pg1, min(pss1))
identical(pg1, min(psd1))

set.seed(290875); pg2 <- pvalue(it2, "global")[1]
set.seed(290875); pss2 <- pvalue(it2, "single-step")
set.seed(290875); psd2 <- pvalue(it2, "step-down")
identical(pg2, min(pss2))
identical(pg2, min(psd2))

set.seed(290875); pg3 <- pvalue(it3, "global")[1]
set.seed(290875); pss3 <- pvalue(it3, "single-step")
set.seed(290875); psd3 <- pvalue(it3, "step-down")
identical(pg3, min(pss3))
identical(pg3, min(psd3))

pg4 <- pvalue(it4, "global")[1]
pss4 <- pvalue(it4, "single-step")
psd4 <- pvalue(it4, "step-down")
identical(pg4, min(pss4))
identical(pg4, min(psd4))

pg5 <- pvalue(it5, "global")[1]
pss5 <- pvalue(it5, "single-step")
psd5 <- pvalue(it5, "step-down")
identical(pg5, min(pss5))
identical(pg5, min(psd5))

pg6 <- pvalue(it6, "global")[1]
pss6 <- pvalue(it6, "single-step")
psd6 <- pvalue(it6, "step-down")
identical(pg6, min(pss6))
identical(pg6, min(psd6))

### adjusted marginal asymptotic p-values

pvalue(it1, "Bonferroni")
pvalue(it1, "Sidak")
pvalue(it1, "Bonferroni-Holm")
pvalue(it1, "Sidak-Holm")

pvalue(it2, "Bonferroni")
pvalue(it2, "Sidak")
pvalue(it2, "Bonferroni-Holm")
pvalue(it2, "Sidak-Holm")

pvalue(it3, "Bonferroni")
pvalue(it3, "Sidak")
pvalue(it3, "Bonferroni-Holm")
pvalue(it3, "Sidak-Holm")

### mcp

YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
                             42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
                             38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
                             31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                  site = factor(c(rep("I", 10), rep("II", 10),
                                  rep("III", 10), rep("IV", 10))))

### permutation based Dunnett
it <- independence_test(length ~ site, data = YOY,
                        xtrafo = mcp_trafo(site = "Dunnett"),
                        distribution = approximate(10000),
                        alternative = "two.sided")
pvalue(it, method = "npmcp")

### asymptotic Dunnett
it <- independence_test(length ~ site, data = YOY,
                        xtrafo = mcp_trafo(site = "Dunnett"),
                        alternative = "two.sided")
pvalue(it, method = "npmcp")

### asymptotic Dunnett, user-defined w/o column names
cm <- rbind("II  vs I" = c(-1, 1, 0, 0),
            "III vs I" = c(-1, 0, 1, 0),
            "IV  vs I" = c(-1, 0, 0, 1))
it <- independence_test(length ~ site, data = YOY,
                        xtrafo = mcp_trafo(site = cm),
                        alternative = "two.sided")
pvalue(it, method = "npmcp")
