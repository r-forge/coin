################################################################################
### Comparisons with examples and output given in the StaXact 9 manual
################################################################################

library("coin")
set.seed(290875)
isequal <- coin:::isequal
table2df <- coin:::table2df
table2df_sym <- coin:::table2df_sym
### <FIXME>  Remove this
od <- setwd("e:/research/coin/pkg/coin_update_StatXact_tests/tests/")
### <\FIXME>

################################################################################
### 6      Two-sample Inference: Related Samples
################################################################################

################################################################################
### 6.4    Sign Test
################################################################################

################################################################################
### 6.4.1  Example: AZT for AIDS, p. 106
### <FIXME>  Rename to azt1, or ???
load("AIDS.rda")
### </FIXME>

diff <- with(AIDS, post - pre)
diff0 <- diff[abs(diff) > 0]
y <- as.vector(t(cbind(as.numeric(diff0 < 0), as.numeric(diff0 > 0))))
lvls <- c("pos", "neg")
x <- factor(rep.int(c("neg", "pos"), length(diff0)), levels = lvls)
b <- gl(length(diff0), 2)

##
## 'symmetry.test.formula':
##

## One-sided asymptotic, p. 107
stao <- symmetry_test(y ~ x | b,
                      alternative = "less")
stopifnot(isequal(round(statistic(stao), 4), -3)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0013)) # p-value

## Two-sided asymptotic, p. 107
sta <- symmetry_test(y ~ x | b)
stopifnot(isequal(round(pvalue(sta), 4), 0.0027)) # p-value

## One-sided exact, p. 107
steo <- symmetry_test(y ~ x | b,
                      distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(steo), 4), 0.0021)) # p-value

## Two-sided exact, p. 107
ste <- symmetry_test(y ~ x | b,
                     distribution = "exact")
stopifnot(isequal(round(pvalue(ste), 4), 0.0042)) # p-value
stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0018)) # point prob

## One-sided approximative, p. 107
stxo <- symmetry_test(y ~ x | b,
                      distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0021 & pci[2] > 0.0021) # p-value

## Two-sided approximative, p. 107
stx <- symmetry_test(y ~ x | b,
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0042 & pci[2] > 0.0042) # p-value

## Clean-up
rm(AIDS, diff, diff0, y, lvls, x, b, stao, sta, steo, ste,
   stxo, stx, pci)
################################################################################


################################################################################
### 6.5    Wilcoxon Signed Rank Test
################################################################################

################################################################################
### 6.5.1  Example: Wilcoxon Signed Rank Test, p. 113
### <FIXME>  Rename to azt1, or ???
load("AIDS.rda")
### </FIXME>

y <- c(AIDS$pre, AIDS$post)
x <- gl(2, nrow(AIDS))
b <- factor(rep.int(1:nrow(AIDS), 2))

##
## 'wilcoxsign_test.formula', without block:
##

## One-sided asymptotic, p. 115
wtao <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Wilcoxon",
                        alternative = "less")
stopifnot(isequal(round(statistic(wtao), 4), -2.8957)) # test statistic
stopifnot(isequal(round(pvalue(wtao), 4), 0.0019)) # p-value

## Two-sided asymptotic, p. 115
wta <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Wilcoxon")
stopifnot(isequal(round(pvalue(wta), 4), 0.0038)) # p-value

## One-sided exact, p. 115
wteo <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Wilcoxon",
                        distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(wteo), 4), 0.0011)) # p-value

## Two-sided exact, p. 115
wte <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Wilcoxon",
                       distribution = "exact")
stopifnot(isequal(round(pvalue(wte), 4), 0.0021)) # p-value
stopifnot(isequal(round(dperm(wte, statistic(wte)), 4), 0.0002)) # point prob

## One-sided approximative, p. 115
wtxo <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Wilcoxon",
                        distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(wtxo), "conf.int")
stopifnot(pci[1] < 0.0011 & pci[2] > 0.0011) # p-value

## Two-sided approximative, p. 115
wtx <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Wilcoxon",
                       distribution = approximate(B = 10000))
pci <- attr(pvalue(wtx), "conf.int")
stopifnot(pci[1] < 0.0021 & pci[2] > 0.0021) # p-value

##
## 'wilcoxsign_test.formula', with block:
##

## One-sided asymptotic, p. 115
wtao2 <- wilcoxsign_test(y ~ x | b, zero.method = "Wilcoxon",
                        alternative = "less")
stopifnot(isequal(statistic(wtao2), statistic(wtao))) # test statistic
stopifnot(isequal(pvalue(wtao2), pvalue(wtao))) # p-value

## Two-sided asymptotic, p. 115
wta2 <- wilcoxsign_test(y ~ x | b, zero.method = "Wilcoxon")
stopifnot(isequal(pvalue(wta2), pvalue(wta))) # p-value

## One-sided exact, p. 115
wteo2 <- wilcoxsign_test(y ~ x | b, zero.method = "Wilcoxon",
                         distribution = "exact", alternative = "less")
stopifnot(isequal(pvalue(wteo2), pvalue(wteo))) # p-value

## Two-sided exact, p. 115
wte2 <- wilcoxsign_test(y ~ x | b, zero.method = "Wilcoxon",
                        distribution = "exact")
stopifnot(isequal(pvalue(wte2), pvalue(wte))) # p-value
stopifnot(isequal(dperm(wte2, statistic(wte2)), dperm(wte, statistic(wte)))) # point prob

## One-sided approximative, p. 115
wtxo2 <- wilcoxsign_test(y ~ x | b, zero.method = "Wilcoxon",
                         distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(wtxo2), "conf.int")
stopifnot(pci[1] < 0.0011 & pci[2] > 0.0011) # p-value

## Two-sided approximative, p. 115
wtx2 <- wilcoxsign_test(y ~ x | b, zero.method = "Wilcoxon",
                        distribution = approximate(B = 10000))
pci <- attr(pvalue(wtx2), "conf.int")
stopifnot(pci[1] < 0.0021 & pci[2] > 0.0021) # p-value

## Clean-up
rm(AIDS, y, x, b, wtao, wta, wteo, wte,
   wtxo, wtx, wtao2, wta2, wteo2, wte2,
   wtxo2, wtx2, pci)
################################################################################


################################################################################
### 6.7    Permutation Test with Arbitrary Scores
################################################################################

################################################################################
### 6.7.1  Example 1: AIDS Data Including the Zeros, p.123
### <FIXME>  Rename to azt1, or ???
load("AIDS.rda")
### </FIXME>

y <- c(AIDS$pre, AIDS$post)
x <- gl(2, nrow(AIDS))
b <- factor(rep.int(1:nrow(AIDS), 2))

##
## 'wilcoxsign_test.formula', without block:
##

## One-sided asymptotic, p. 125
wtao <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Pratt",
                        alternative = "less")
stopifnot(isequal(round(statistic(wtao), 4), -3.0023)) # test statistic
stopifnot(isequal(round(pvalue(wtao), 4), 0.0013)) # p-value

## Two-sided asymptotic, p. 125
wta <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Pratt")
stopifnot(isequal(round(pvalue(wta), 4), 0.0027)) # p-value

## One-sided exact, p. 125
wteo <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Pratt",
                        distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(wteo), 4), 0.0008)) # p-value

## Two-sided exact, p. 125
wte <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Pratt",
                       distribution = "exact")
stopifnot(isequal(round(pvalue(wte), 4), 0.0016)) # p-value
stopifnot(isequal(round(dperm(wte, statistic(wte)), 4), 0.0001)) # point prob

## One-sided approximative, p. 125
wtxo <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Pratt",
                        distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(wtxo), "conf.int")
stopifnot(pci[1] < 0.0008 & pci[2] > 0.0008) # p-value

## Two-sided approximativet, p. 125
wtx <- wilcoxsign_test(post ~ pre, data = AIDS, zero.method = "Pratt",
                       distribution = approximate(B = 10000))
pci <- attr(pvalue(wtx), "conf.int")
stopifnot(pci[1] < 0.0016 & pci[2] > 0.0016) # p-value

##
## 'wilcoxsign_test.formula', with block:
##

## One-sided asymptotic, p. 125
wtao2 <- wilcoxsign_test(y ~ x | b, zero.method = "Pratt",
                        alternative = "less")
stopifnot(isequal(statistic(wtao2), statistic(wtao))) # test statistic
stopifnot(isequal(pvalue(wtao2), pvalue(wtao))) # p-value

## Two-sided asymptotic, p. 125
wta2 <- wilcoxsign_test(y ~ x | b, zero.method = "Pratt")
stopifnot(isequal(pvalue(wta2), pvalue(wta))) # p-value

## One-sided exact, p. 125
wteo2 <- wilcoxsign_test(y ~ x | b, zero.method = "Pratt",
                         distribution = "exact", alternative = "less")
stopifnot(isequal(pvalue(wteo2), pvalue(wteo))) # p-value

## Two-sided exact, p. 125
wte2 <- wilcoxsign_test(y ~ x | b, zero.method = "Pratt",
                        distribution = "exact")
stopifnot(isequal(pvalue(wte2), pvalue(wte))) # p-value
stopifnot(isequal(dperm(wte2, statistic(wte2)), dperm(wte, statistic(wte)))) # point prob

## One-sided approximative, p. 125
wtxo2 <- wilcoxsign_test(y ~ x | b, zero.method = "Pratt",
                         distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(wtxo2), "conf.int")
stopifnot(pci[1] < 0.0011 & pci[2] > 0.0011) # p-value

## Two-sided approximative, p. 125
wtx2 <- wilcoxsign_test(y ~ x | b, zero.method = "Pratt",
                        distribution = approximate(B = 10000))
pci <- attr(pvalue(wtx2), "conf.int")
stopifnot(pci[1] < 0.0021 & pci[2] > 0.0021) # p-value

## Clean-up
rm(AIDS, y, x, b, wtao, wta, wteo, wte,
   wtxo, wtx, wtao2, wta2, wteo2, wte2,
   wtxo2, wtx2, pci)
################################################################################

################################################################################
### 6.7.2  Example 2: AIDS Data with Raw Differences, p. 125
### <FIXME>  Rename to azt1, or ???
load("AIDS.rda")
### </FIXME>

diff <- with(AIDS, post - pre)
y <- as.vector(t(cbind(abs(diff) * (diff < 0), abs(diff) * (diff >= 0))))
x <- factor(rep(c("neg", "pos"), length(diff)), levels = c("pos", "neg"))
b <- gl(length(diff), 2)

##
## symmetry_test.formula:
##

## One-sided asymptotic, p. 126
stao <- symmetry_test(y ~ x | b,
                      alternative = "less")
stopifnot(isequal(round(statistic(stao), 4), -1.7072)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0439)) # p-value

## Two-sided asymptotic, p. 126
sta <- symmetry_test(y ~ x | b)
stopifnot(isequal(round(pvalue(sta), 4), 0.0878)) # p-value

## One-sided exact, p. 126
steo <- symmetry_test(y ~ x | b,
                      distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(steo), 4), 0.0011)) # p-value

## Two-sided exact, p. 126
ste <- symmetry_test(y ~ x | b,
                     distribution = "exact")
stopifnot(isequal(round(pvalue(ste), 4), 0.0021)) # p-value
stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0000)) # point prob

## One-sided approximative, p. 126
stxo <- symmetry_test(y ~ x | b,
                      distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0011 & pci[2] > 0.0011) # p-value

## Two-sided approximative, p. 126
stx <- symmetry_test(y ~ x | b,
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0021 & pci[2] > 0.0021) # p-value

## Clean-up
rm(AIDS, diff, y, x, b, stao, sta, steo, ste,
   stxo, stx, pci)
################################################################################


################################################################################
### 6.8    McNemar's Test
################################################################################

################################################################################
### 6.8.3  Example: Voters' Preference, p. 129
load("VOTE.rda")
diag(vote) <- 0 # no contribution from the diagonal elements

dta <- table2df_sym(vote)
y <- dta$response
lvls <- c("Preference.Before.TV.Debate", "Preference.After.TV.Debate")
x <- factor(dta$groups, levels = lvls)
b <- factor(rep.int(seq_len(sum(vote)), 2))

##
## 'synmmetry_test.formula':
##

## One-sided asymptotic, p. 130
stao <- symmetry_test(y ~ x | b,
                      alternative = "greater")
stopifnot(isequal(round(statistic(stao), 4), 1.3416)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0899)) # p-value

## Two-sided asymptotic, p. 130
sta <- symmetry_test(y ~ x | b)
stopifnot(isequal(round(pvalue(sta), 4), 0.1797)) # p-value

## ## One-sided exact, p. 130
## steo <- symmetry_test(y ~ x | b,
##                       distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(steo), 4), 0.1316)) # p-value

## ## Two-sided exact, p. 130
## ste <- symmetry_test(y ~ x | b,
##                      distribution = "exact")
## stopifnot(isequal(round(pvalue(ste), 4), 0.2632)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0739)) # point prob

## One-sided approximative, p. 130
stxo <- symmetry_test(y ~ x | b,
                      distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.1316 & pci[2] > 0.1316) # p-value

## Two-sided approximative, p. 130
stx <- symmetry_test(y ~ x | b,
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.2632 & pci[2] > 0.2632) # p-value

##
## 'symmetry_test.table':
##

## One-sided asymptotic, p. 130
stao2 <- symmetry_test(vote,
                       xtrafo = function(data)
                           trafo(data, factor_trafo = function(x)
                               f_trafo(factor(x, levels = lvls))),
                       alternative = "greater")
stopifnot(isequal(statistic(stao2), statistic(stao))) # test statistic
stopifnot(isequal(pvalue(stao2), pvalue(stao))) # p-value

## Two-sided asymptotic, p. 130
sta2 <- symmetry_test(vote,
                      xtrafo = function(data)
                          trafo(data, factor_trafo = function(x)
                              f_trafo(factor(x, levels = lvls))))
stopifnot(isequal(pvalue(sta2), pvalue(sta))) # p-value

## ## One-sided exact, p. 130
## steo2 <- symmetry_test(vote,
##                        xtrafo = function(data)
##                            trafo(data, factor_trafo = function(x)
##                                f_trafo(factor(x, levels = lvls))),
##                        distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(steo2), pvalue(steo))) # p-value

## ## Two-sided exact, p. 130
## ste2 <- symmetry_test(vote,
##                       xtrafo = function(data)
##                           trafo(data, factor_trafo = function(x)
##                               f_trafo(factor(x, levels = lvls))),
##                       distribution = "exact")
## stopifnot(isequal(pvalue(ste2), pvalue(ste))) # p-value
## stopifnot(isequal(dperm(ste2, statistic(ste2)), dperm(ste, statistic(ste)))) # point prob

## One-sided approximative, p. 130
stxo2 <- symmetry_test(vote,
                       xtrafo = function(data)
                           trafo(data, factor_trafo = function(x)
                               f_trafo(factor(x, levels = lvls))),
                       distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo2), "conf.int")
stopifnot(pci[1] < 0.1316 & pci[2] > 0.1316) # p-value

## Two-sided approximative, p. 130
stx2 <- symmetry_test(vote,
                      xtrafo = function(data)
                          trafo(data, factor_trafo = function(x)
                              f_trafo(factor(x, levels = lvls))),
                      distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.2632 & pci[2] > 0.2632) # p-value

##
## 'mh_test.formula':
##

## Two-sided asymptotic, p. 130
mta <- mh_test(y ~ x | b)
stopifnot(isequal(sqrt(statistic(mta)), unname(statistic(sta)))) # test statistic
stopifnot(isequal(pvalue(mta), pvalue(sta))) # test statistic

## Two-sided approximative, p. 130
mtx <- mh_test(y ~ x | b,
               distribution = approximate(B = 10000))
pci <- attr(pvalue(mtx), "conf.int")
stopifnot(pci[1] < 0.2632 & pci[2] > 0.2632) # p-value

##
## 'mh_test.table':
##

## Two-sided asymptotic, p. 130
mta2 <- mh_test(vote)
stopifnot(isequal(sqrt(statistic(mta2)), unname(statistic(sta)))) # test statistic

## Two-sided approximative, p. 130
mtx2 <- mh_test(vote,
               distribution = approximate(B = 10000))
pci <- attr(pvalue(mtx2), "conf.int")
stopifnot(pci[1] < 0.2632 & pci[2] > 0.2632) # p-value

## Clean-up
rm(vote, dta, y, lvls, x, b, stao, sta, #steo, ste,
   stxo, stx, stao2, sta2, #steo2, ste2,
   stxo2, stx2, mta, mtx, mta2, mtx2, pci)
################################################################################


################################################################################
### 6.9    Marginal Homogeneity Test
################################################################################

################################################################################
### 6.9.1  Example 1: Matched Case-Control Study of Endometrial Cancer, p. 131
load("ENDOMET.rda")
diag(endomet) <- 0 # no contribution from the diagonal elements
scrs <- c(0, 0.2, 0.5125, 0.7)

dta <- table2df_sym(endomet)
y <- dta$response
x <- factor(dta$groups)
b <- factor(rep.int(seq_len(sum(endomet)), 2))

##
## 'symmetry_test.formula':
##

## One-sided asymptotic, p. 134
stao <- symmetry_test(y ~ x | b, scores = list(y = scrs),
                      alternative = "greater")
stopifnot(isequal(round(statistic(stao), 4), 3.7346)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0001)) # p-value

## Two-sided asymptotic, p. 134
sta <- symmetry_test(y ~ x | b, scores = list(y = scrs))
stopifnot(isequal(round(pvalue(sta), 4), 0.0002)) # p-value

## ## One-sided exact, p. 134
## steo <- symmetry_test(y ~ x | b, scores = list(y = scrs),
##                       distribution = exact(fact = 1e5), alternative = "greater")
## stopifnot(isequal(round(pvalue(steo), 4), 0.0001)) # p-value

## ## Two-sided exact, p. 134
## ste <- symmetry_test(y ~ x | b, scores = list(y = scrs),
##                      distribution = exact(fact = 1e5))
## stopifnot(isequal(round(pvalue(ste), 4), 0.0001)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0000)) # point prob

## One-sided approximative, p. 134
stxo <- symmetry_test(y ~ x | b, scores = list(y = scrs),
                      distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0001 & pci[2] > 0.0001) # p-value

## Two-sided approximative, p. 134
stx <- symmetry_test(y ~ x | b, scores = list(y = scrs),
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0001 & pci[2] > 0.0001) # p-value

##
## 'symmetry_test.table':
##

## One-sided asymptotic, p. 134
stao2 <- symmetry_test(endomet, scores = list(response = scrs),
                       alternative = "greater")
stopifnot(isequal(statistic(stao2), statistic(stao))) # test statistic
stopifnot(isequal(pvalue(stao2), pvalue(stao))) # p-value

## Two-sided asymptotic, p. 134
sta2 <- symmetry_test(endomet, scores = list(response = scrs))
stopifnot(isequal(pvalue(sta2), pvalue(sta))) # p-value

## ## One-sided exact, p. 134
## steo2 <- symmetry_test(endomet, scores = list(response = scrs),
##                        distribution = exact(fact = 1e5), alternative = "greater")
## stopifnot(isequal(pvalue(steo2), pvalue(steo))) # p-value

## ## Two-sided exact, p. 134
## ste2 <- symmetry_test(endomet, scores = list(response = scrs),
##                       distribution = exact(fact = 1e5))
## stopifnot(isequal(pvalue(ste2), pvalue(ste))) # p-value
## stopifnot(isequal(dperm(ste2, statistic(ste2)), dperm(ste, statistic(ste)))) # point prob

## One-sided approximative, p. 134
stxo2 <- symmetry_test(endomet, scores = list(response = scrs),
                       distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo2), "conf.int")
stopifnot(pci[1] < 0.0001 & pci[2] > 0.0001) # p-value

## Two-sided approximative, p. 134
stx2 <- symmetry_test(endomet, scores = list(response = scrs),
                      distribution = approximate(B = 10000))
pci <- attr(pvalue(stx2), "conf.int")
stopifnot(pci[1] < 0.0001 & pci[2] > 0.0001) # p-value

##
## 'mh_test.formula':
##

## Two-sided asymptotic, p. 134
mta <- mh_test(y ~ x | b, scores = list(y = scrs))
stopifnot(isequal(sqrt(statistic(mta)), unname(statistic(sta)))) # test statistic
## <FIXME>  Why is this ???  </FIXME>

## Two-sided approximative, p. 134
mtx <- mh_test(y ~ x | b, scores = list(y = scrs),
               distribution = approximate(B = 10000))
pci <- attr(pvalue(mtx), "conf.int")
stopifnot(pci[1] < 0.0001 & pci[2] > 0.0001) # p-value
## <FIXME>  Why is this ???  </FIXME>

##
## 'mh_test.table':
##

## Two-sided asymptotic, p. 134
mta2 <- mh_test(endomet, scores = list(response = scrs))
stopifnot(isequal(sqrt(statistic(mta2)), unname(statistic(sta)))) # test statistic

## Two-sided approximative, p. 134
mtx2 <- mh_test(endomet, scores = list(response = scrs),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(mtx2), "conf.int")
stopifnot(pci[1] < 0.0001 & pci[2] > 0.0001) # p-value

## Clean-up
rm(endomet, scrs, dta, y, x, b, stao, sta, #steo, ste,
   stxo, stx, stao2, sta2, #steo2, ste2,
   stxo2, stx2, mta, mtx, mta2, mtx2, pci)
################################################################################

################################################################################
### 6.9.2  Example 2: Pap-Smear Classification by Two Pathologists, p. 134
load("PAP-SMR.rda")
diag(pap_smr) <- 0 # no contribution from the diagonal elements

dta <- table2df_sym(pap_smr)
y <- dta$response
x <- factor(dta$groups)
b <- factor(rep.int(seq_len(sum(pap_smr)), 2))

##
## 'symmetry_test.formula':
##

## One-sided asymptotic, p. 136
stao <- symmetry_test(y ~ x | b, scores = list(y = 1:5),
                      alternative = "greater")
stopifnot(isequal(round(statistic(stao), 4), 1.1523)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.1246)) # p-value

## Two-sided asymptotic, p. 136
sta <- symmetry_test(y ~ x | b, scores = list(y = 1:5))
stopifnot(isequal(round(pvalue(sta), 4), 0.2492)) # p-value

## ## One-sided exact, p. 136
## steo <- symmetry_test(y ~ x | b, scores = list(y = 1:5),
##                       alternative = "greater", distribution = "exact")
## stopifnot(isequal(round(pvalue(steo), 4), 0.1536)) # p-value

## ## Two-sided exact, p. 136
## ste <- symmetry_test(y ~ x | b, scores = list(y = 1:5),
##                      distribution = "exact")
## stopifnot(isequal(round(pvalue(ste), 4), 0.3073)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0531)) # point prob

## One-sided approximative, p. 136
stxo <- symmetry_test(y ~ x | b, scores = list(y = 1:5),
                      distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.1536 & pci[2] > 0.1536)

## Two-sided approximative, p. 136
stx <- symmetry_test(y ~ x | b, scores = list(y = 1:5),
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.3073 & pci[2] > 0.3073)

##
## 'symmetry_test.table':
##

## One-sided asymptotic, p. 136
stao2 <- symmetry_test(pap_smr, scores = list(response = 1:5),
                       alternative = "greater")
stopifnot(isequal(statistic(stao2), statistic(stao))) # test statistic
stopifnot(isequal(pvalue(stao2), pvalue(stao))) # p-value

## Two-sided asymptotic, p. 136
sta2 <- symmetry_test(pap_smr, scores = list(response = 1:5))
stopifnot(isequal(pvalue(sta2), pvalue(sta))) # p-value

## ## One-sided exact, p. 136
## steo2 <- symmetry_test(pap_smr, scores = list(response = 1:5),
##                        alternative = "greater", distribution = "exact")
## stopifnot(isequal(pvalue(steo2), pvalue(steo))) # p-value

## ## Two-sided exact, p. 136
## ste2 <- symmetry_test(pap_smr, scores = list(response = 1:5),
##                       distribution = "exact")
## stopifnot(isequal(pvalue(ste2), pvalue(ste))) # p-value
## stopifnot(isequal(dperm(ste2, statistic(ste2)), dperm(ste, statistic(ste)))) # point prob

## One-sided approximative, p. 136
stxo2 <- symmetry_test(pap_smr, scores = list(response = 1:5),
                       distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo2), "conf.int")
stopifnot(pci[1] < 0.1536 & pci[2] > 0.1536)

## Two-sided approximative, p. 136
stx2 <- symmetry_test(pap_smr, scores = list(response = 1:5),
                      distribution = approximate(B = 10000))
pci <- attr(pvalue(stx2), "conf.int")
stopifnot(pci[1] < 0.3073 & pci[2] > 0.3073)

##
## 'mh_test.formula':
##

## Two-sided asymptotic, p. 136
mta <- mh_test(y ~ x | b, scores = list(y = 1:5))
stopifnot(isequal(sqrt(statistic(mta)), unname(statistic(sta)))) # test statistic
## <FIXME>  Why is this ???  </FIXME>

## Two-sided approximative, p. 136
mtx <- mh_test(y ~ x | b, scores = list(y = 1:5),
               distribution = approximate(B = 10000))
pci <- attr(pvalue(mtx), "conf.int")
stopifnot(pci[1] < 0.3073 & pci[2] > 0.3073) # p-value
## <FIXME>  Why is this ???  </FIXME>

##
## 'mh_test.table':
##

## Two-sided asymptotic, p. 136
mta2 <- mh_test(pap_smr, scores = list(response = 1:5))
stopifnot(isequal(sqrt(statistic(mta2)), unname(statistic(sta)))) # test statistic

## Two-sided approximative, p. 136
mtx2 <- mh_test(pap_smr, scores = list(response = 1:5),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(mtx2), "conf.int")
stopifnot(pci[1] < 0.3073 & pci[2] > 0.3073) # p-value

## Clean-up
rm(pap_smr, dta, y, x, b, stao, sta, #steo, ste,
   stxo, stx, stao2, sta2, #steo2, ste2,
   stxo2, stx2, mta, mtx, mta2, mtx2, pci)
################################################################################


################################################################################
### 7      Two-sample Inference: Independent Samples
################################################################################

################################################################################
### 7.4    Wilcoxon-Mann-Whitney Test
################################################################################

################################################################################
### 7.4.3  Example with Unstratified Data: Diastolic Blod Pressure, p.154
load("BLOODPR.rda")

##
## 'wilcox_test.formula':
##

## One-sided asymptotic, p. 156
wtao <- wilcox_test(response ~ group, data = bloodpr,
                    alternative = "greater")
stopifnot(isequal(round(statistic(wtao), 4), 1.7205)) # test statistic
stopifnot(isequal(round(pvalue(wtao), 4), 0.0427)) # p-value

## Two-sided asymptotic, p. 156
wta <- wilcox_test(response ~ group, data = bloodpr)
stopifnot(isequal(round(pvalue(wta), 4), 0.0853)) # p-value

## One-sided exact, p. 156
wteo <- wilcox_test(response ~ group, data = bloodpr,
                    distribution = "exact", alternative = "greater")
stopifnot(isequal(round(pvalue(wteo), 4), 0.0542)) # p-value

## Two-sided exact, p. 156
wte <- wilcox_test(response ~ group, data = bloodpr,
                    distribution = "exact")
stopifnot(isequal(round(pvalue(wte), 4), 0.0989)) # p-value
stopifnot(isequal(round(dperm(wte, statistic(wte)), 4), 0.0190)) # point prob

## One-sided approximative, p. 156, (157--158)
wtxo <- wilcox_test(response ~ group, data = bloodpr,
                    distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(wtxo), "conf.int")
stopifnot(pci[1] < 0.0542 & pci[2] > 0.0542) # p-value

## Two-sided approximative, p. 156, (157--158)
wtx <- wilcox_test(response ~ group, data = bloodpr,
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(wtx), "conf.int")
stopifnot(pci[1] < 0.0989 & pci[2] > 0.0989) # p-value

## Clean-up
rm(bloodpr, wtao, wta, wteo, wte, wtxo, wtx, pci)
################################################################################

################################################################################
### 7.4.4  Example with Stratified Data: Employment Discrimination, p.158
load("DISCRIM.rda")

##
## 'wilcox_test.formula':
##

## One-sided asymptotic, p. 160
wtao <- wilcox_test(salary ~ gender | year, data = discrim,
                    alternative = "less")
stopifnot(isequal(round(statistic(wtao), 4), -1.8970)) # test statistic
stopifnot(isequal(round(pvalue(wtao), 4), 0.0289)) # p-value

## Two-sided asymptotic, p. 160
wta <- wilcox_test(salary ~ gender | year, data = discrim)
stopifnot(isequal(round(pvalue(wta), 4), 0.0578)) # p-value

## ## One-sided exact, p. 160
## wteo <- wilcox_test(salary ~ gender | year, data = discrim,
##                     distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(wteo), 4), 0.0400)) # p-value

## ## Two-sided exact, p. 160
## wte <- wilcox_test(salary ~ gender | year, data = discrim,
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(wte), 4), 0.0400)) # p-value
## stopifnot(isequal(round(dperm(wte, statistic(wte)), 4), 0.0400)) # point prob

## One-sided approximative, p. 160
wtxo <- wilcox_test(salary ~ gender | year, data = discrim,
                    distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(wtxo), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Two-sided approximative, p. 160
wtx <- wilcox_test(salary ~ gender | year, data = discrim,
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(wtx), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Clean-up
rm(discrim, wtao, wta, #wteo, wte,
   wtxo, wtx, pci)
################################################################################


################################################################################
### 7.6    Normal Scores Test
################################################################################

################################################################################
### 7.6.3  Example with Unstratified Data: Diastolic Blod Pressure, p.165
load("BLOODPR.rda")

##
## 'normal_test.formula':
##

## One-sided asymptotic, p. 167
ntao <- normal_test(response ~ group, data = bloodpr, ties.method = "average",
                    alternative = "greater")
stopifnot(isequal(round(statistic(ntao), 4), 1.7885)) # test statistic
stopifnot(isequal(round(pvalue(ntao), 4), 0.0368)) # p-value

## Two-sided asymptotic, p. 167
nta <- normal_test(response ~ group, data = bloodpr, ties.method = "average")
stopifnot(isequal(round(pvalue(nta), 4), 0.0737)) # p-value

## One-sided exact, p. 167
nteo <- normal_test(response ~ group, data = bloodpr, ties.method = "average",
                    distribution = "exact", alternative = "greater")
stopifnot(isequal(round(pvalue(nteo), 4), 0.0462)) # p-value

## Two-sided exact, p. 167
nte <- normal_test(response ~ group, data = bloodpr, ties.method = "average",
                   distribution = "exact")
stopifnot(isequal(round(pvalue(nte), 4), 0.0799)) # p-value
stopifnot(isequal(round(dperm(nte, statistic(nte)), 4), 0.0176)) # point prob
## <FIXME>  Why is this???  </FIXME>

## One-sided approximative, p. 167, (168--169)
ntxo <- normal_test(response ~ group, data = bloodpr, ties.method = "average",
                    distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(ntxo), "conf.int")
stopifnot(pci[1] < 0.0462 & pci[2] > 0.0462) # p-value

## Two-sided approximative, p. 167, (168--169)
ntx <- normal_test(response ~ group, data = bloodpr, ties.method = "average",
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(ntx), "conf.int")
stopifnot(pci[1] < 0.0799 & pci[2] > 0.0799) # p-value

## Clean-up
rm(bloodpr, ntao, nta, nteo, nte, ntxo, ntx, pci)
################################################################################

################################################################################
### 7.6.4  Example with Stratified Data: Employment Discrimination, p.169
load("DISCRIM.rda")

##
## 'normal_test.formula':
##

## One-sided asymptotic, p. 170
ntao <- normal_test(salary ~ gender | year, data = discrim, ties.method = "average",
                    alternative = "less")
stopifnot(isequal(round(statistic(ntao), 4), -1.8018)) # test statistic
stopifnot(isequal(round(pvalue(ntao), 4), 0.0358)) # p-value

## Two-sided asymptotic, p. 170
nta <- normal_test(salary ~ gender | year, data = discrim, ties.method = "average")
stopifnot(isequal(round(pvalue(nta), 4), 0.0716)) # p-value

## ## One-sided exact, p. 170
## nteo <- normal_test(salary ~ gender | year, data = discrim, ties.method = "average",
##                     distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(nteo), 4), 0.0400)) # p-value

## ## Two-sided exact, p. 170
## nte <- normal_test(salary ~ gender | year, data = discrim, ties.method = "average",
##                     distribution = "exact")
## stopifnot(isequal(round(pvalue(nte), 4), 0.0400)) # p-value
## stopifnot(isequal(round(dperm(nte, statistic(nte)), 4), 0.0400)) # point prob

## One-sided approximative, p. 170
ntxo <- normal_test(salary ~ gender | year, data = discrim, ties.method = "average",
                    distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(ntxo), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Two-sided approximative, p. 170
ntx <- normal_test(salary ~ gender | year, data = discrim, ties.method = "average",
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(ntx), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Clean-up
rm(discrim, ntao, nta, #nteo, nte,
   ntxo, ntx, pci)
################################################################################


################################################################################
### 7.7    Savage Scores Test
################################################################################

################################################################################
### 7.7.3  Unstratified Data: Lung Cancer Example, p. 171
load("CANCER.rda")

##
## 'surv_test.formula':
##

## One-sided asymptotic, p. 172
stao <- surv_test(Surv(response) ~ drug, data = cancer, ties.method = "average",
                  alternative = "greater")
stopifnot(isequal(round(statistic(stao), 4), 2.0437)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0205)) # p-value

## Two-sided asymptotic, p. 172
sta <- surv_test(Surv(response) ~ drug, data = cancer, ties.method = "average")
stopifnot(isequal(round(pvalue(sta), 4), 0.0410)) # p-value

## One-sided exact, p. 172
steo <- surv_test(Surv(response) ~ drug, data = cancer, ties.method = "average",
                  distribution = "exact", alternative = "greater")
stopifnot(isequal(round(pvalue(steo), 4), 0.0220)) # p-value

## Two-sided exact, p. 172
ste <- surv_test(Surv(response) ~ drug, data = cancer, ties.method = "average",
                 distribution = "exact")
stopifnot(isequal(round(pvalue(ste), 4), 0.0260)) # p-value
stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0010)) # point prob
## <FIXME>  Why is this???  </FIXME>

## One-sided approximative, p. 172
stxo <- surv_test(Surv(response) ~ drug, data = cancer, ties.method = "average",
                  distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0220 & pci[2] > 0.0220) # p-value

## Two-sided approximative, p. 172
stx <- surv_test(Surv(response) ~ drug, data = cancer, ties.method = "average",
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0260 & pci[2] > 0.0260) # p-value

## Clean-up
rm(cancer, stao, sta, steo, ste, stxo, stx, pci)
################################################################################

################################################################################
### 7.7.4  Stratified Data Example, p.172
load("DISCRIM.rda")

##
## 'surv_test.formula':
##

## One-sided asymptotic, p. 174
stao <- surv_test(Surv(salary) ~ gender | year, data = discrim, ties.method = "average",
                  alternative = "less")
stopifnot(isequal(round(statistic(stao), 4), -1.5357)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0623)) # p-value

## Two-sided asymptotic, p. 174
sta <- surv_test(Surv(salary) ~ gender | year, data = discrim, ties.method = "average")
stopifnot(isequal(round(pvalue(sta), 4), 0.1246)) # p-value

## ## One-sided exact, p. 174
## steo <- surv_test(Surv(salary) ~ gender | year, data = discrim, ties.method = "average",
##                   distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(steo), 4), 0.0400)) # p-value

## ## Two-sided exact, p. 174
## ste <- surv_test(Surv(salary) ~ gender | year, data = discrim, ties.method = "average",
##                  distribution = "exact")
## stopifnot(isequal(round(pvalue(ste), 4), 0.0400)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0400)) # point prob

## One-sided approximative, p. 174
stxo <- surv_test(Surv(salary) ~ gender | year, data = discrim, ties.method = "average",
                  distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Two-sided approximative, p. 174
stx <- surv_test(Surv(salary) ~ gender | year, data = discrim, ties.method = "average",
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Clean-up
rm(discrim, stao, sta, #steo, ste,
   stxo, stx, pci)
################################################################################


## ################################################################################
## ### 7.8    Siegel-Tukey Test
## ################################################################################

## ################################################################################
## ### 7.8.3  Example: Cereal Packaging Machine, p. 176
## load("CER-PKG.rda")

## ##
## ## 'siegel_test.formula':
## ##

## ## One-sided asymptotic, p. 177
## stao <- siegel_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                     alternative = "less")
## stopifnot(isequal(round(statistic(stao), 4), -1.9838)) # test statistic
## stopifnot(isequal(round(pvalue(stao), 4), 0.0236)) # p-value

## ## Two-sided asymptotic, p. 177
## sta <- siegel_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average")
## stopifnot(isequal(round(pvalue(sta), 4), 0.0473)) # p-value

## ## One-sided exact, p. 177
## steo <- siegel_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                     distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(steo), 4), 0.0253)) # p-value

## ## Two-sided exact, p. 177
## ste <- siegel_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(ste), 4), 0.0530)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0051)) # point prob
## ## <FIXME>  Why is this???  </FIXME>

## ## One-sided approximative, p. 177
## stxo <- siegel_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                     distribution = approximate(B = 10000), alternative = "less")
## pci <- attr(pvalue(stxo), "conf.int")
## stopifnot(pci[1] < 0.0253 & pci[2] > 0.0253) # p-value

## ## Two-sided approximative, p. 177
## stx <- siegel_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                    distribution = approximate(B = 10000))
## pci <- attr(pvalue(stx), "conf.int")
## stopifnot(pci[1] < 0.0530 & pci[2] > 0.0530) # p-value

## ## Clean-up
## rm(cer_pkg, stao, sta, steo, ste, stxo, stx, pci)
## ################################################################################


################################################################################
### 7.9    Ansari-Bradley Test
################################################################################

################################################################################
### 7.9.3  Example: Cereal Packaging Machine, p. 180
load("CER-PKG.rda")

##
## 'ansari_test.formula':
##

## One-sided asymptotic, p. 182
atao <- ansari_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
                    alternative = "less")
stopifnot(isequal(round(statistic(atao), 4), -1.9983)) # test statistic
stopifnot(isequal(round(pvalue(atao), 4), 0.0228)) # p-value
## <FIXME>  Why is this???  </FIXME>

## Two-sided asymptotic, p. 182
ata <- ansari_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average")
stopifnot(isequal(round(pvalue(ata), 4), 0.0457)) # p-value

## One-sided exact, p. 182
ateo <- ansari_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
                    distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(ateo), 4), 0.0253)) # p-value
## <FIXME>  Why is this???  </FIXME>

## Two-sided exact, p. 182
ate <- ansari_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
                   distribution = "exact")
stopifnot(isequal(round(pvalue(ate), 4), 0.0581)) # p-value
stopifnot(isequal(round(dperm(ate, statistic(ate)), 4), 0.0051)) # point prob

## One-sided approximative, p. 182
atxo <- ansari_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
                    distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(atxo), "conf.int")
stopifnot(pci[1] < 0.0253 & pci[2] > 0.0253) # p-value
## <FIXME>  Why is this???  </FIXME>

## Two-sided approximative, p. 182
atx <- ansari_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(atx), "conf.int")
stopifnot(pci[1] < 0.0581 & pci[2] > 0.0581) # p-value

## Clean-up
rm(cer_pkg, atao, ata, ateo, ate, atxo, atx, pci)
################################################################################


## ################################################################################
## ### 7.10    Klotz Test
## ################################################################################

## ################################################################################
## ### 7.10.3  Example: Cereal Packaging Machine, p. 185
## load("CER-PKG.rda")

## ##
## ## 'klotz_test.formula':
## ##

## ## One-sided asymptotic, p. 186
## ktao <- klotz_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                    alternative = "greater")
## stopifnot(isequal(round(statistic(ktao), 4), 2.3082)) # test statistic
## stopifnot(isequal(round(pvalue(ktao), 4), 0.0105)) # p-value

## ## Two-sided asymptotic, p. 186
## kta <- klotz_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average")
## stopifnot(isequal(round(pvalue(kta), 4), 0.0210)) # p-value

## ## One-sided exact, p. 186
## kteo <- klotz_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                    distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(kteo), 4), 0.0101)) # p-value

## ## Two-sided exact, p. 186
## kte <- klotz_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                   distribution = "exact")
## stopifnot(isequal(round(pvalue(kte), 4), 0.0101)) # p-value
## stopifnot(isequal(round(dperm(kte, statistic(kte)), 4), 0.0051)) # point prob
## ## <FIXME>  Why is this???  </FIXME>

## ## One-sided approximative, p. 186
## ktxo <- klotz_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                    distribution = approximate(B = 10000), alternative = "greater")
## pci <- attr(pvalue(ktxo), "conf.int")
## stopifnot(pci[1] < 0.0101 & pci[2] > 0.0101) # p-value

## ## Two-sided approximative, p. 186
## ktx <- klotz_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                   distribution = approximate(B = 10000))
## pci <- attr(pvalue(ktx), "conf.int")
## stopifnot(pci[1] < 0.0101 & pci[2] > 0.0101) # p-value

## ## Clean-up
## rm(cer_pkg, ktao, kta, kteo, kte, ktxo, ktx, pci)
## ################################################################################


## ################################################################################
## ### 7.11    Mood Test
## ################################################################################

## ################################################################################
## ### 7.11.3  Example: Cereal Packaging Machine, p. 189
## load("CER-PKG.rda")

## ##
## ## 'mood_test.formula':
## ##

## ## One-sided asymptotic, p. 190
## mtao <- mood_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                    alternative = "greater")
## stopifnot(isequal(round(statistic(mtao), 4), 2.2715)) # test statistic
## stopifnot(isequal(round(pvalue(mtao), 4), 0.0116)) # p-value

## ## Two-sided asymptotic, p. 190
## mta <- mood_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average")
## stopifnot(isequal(round(pvalue(mta), 4), 0.0231)) # p-value

## ## One-sided exact, p. 190
## mteo <- mood_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                    distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(mteo), 4), 0.0126)) # p-value

## ## Two-sided exact, p. 190
## mte <- mood_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                   distribution = "exact")
## stopifnot(isequal(round(pvalue(mte), 4), 0.0202)) # p-value
## stopifnot(isequal(round(dperm(mte, statistic(mte)), 4), 0.0051)) # point prob
## ## <FIXME>  Why is this???  </FIXME>

## ## One-sided approximative, p. 190
## mtxo <- mood_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                    distribution = approximate(B = 10000), alternative = "greater")
## pci <- attr(pvalue(mtxo), "conf.int")
## stopifnot(pci[1] < 0.0126 & pci[2] > 0.0126) # p-value

## ## Two-sided approximative, p. 190
## mtx <- mood_test(cerealwt ~ machin, data = cer_pkg, ties.method = "average",
##                   distribution = approximate(B = 10000))
## pci <- attr(pvalue(mtx), "conf.int")
## stopifnot(pci[1] < 0.0202 & pci[2] > 0.0202) # p-value

## ## Clean-up
## rm(cer_pkg, mtao, mta, mteo, mte, mtxo, mtx, pci)
## ################################################################################


## ################################################################################
## ### 7.12    Conover Test
## ################################################################################

## ################################################################################
## ### 7.12.3  Example: Tire Failure Data, p. 193
## load("FAILURE.rda")

## ##
## ## 'conover_test.formula':
## ##

## ## One-sided asymptotic, p. 195
## ctao <- conover_test(ftime ~ tyretype, data = failure, ties.method = "average",
##                      alternative = "less")
## stopifnot(isequal(round(statistic(ctao), 4), -1.5274)) # test statistic
## stopifnot(isequal(round(pvalue(ctao), 4), 0.0633)) # p-value

## ## Two-sided asymptotic, p. 195
## cta <- conover_test(ftime ~ tyretype, data = failure, ties.method = "average")
## stopifnot(isequal(round(pvalue(cta), 4), 0.1267)) # p-value

## ## One-sided exact, p. 195
## cteo <- conover_test(ftime ~ tyretype, data = failure, ties.method = "average",
##                      distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(cteo), 4), 0.0634)) # p-value

## ## Two-sided exact, p. 195
## cte <- conover_test(ftime ~ tyretype, data = failure, ties.method = "average",
##                     distribution = "exact")
## stopifnot(isequal(round(pvalue(cte), 4), 0.1300)) # p-value
## stopifnot(isequal(round(dperm(cte, statistic(cte)), 4), 0.0007)) # point prob
## ## <FIXME>  Why is this???  </FIXME>

## ## One-sided approximative, p. 195
## ctxo <- conover_test(ftime ~ tyretype, data = failure, ties.method = "average",
##                      distribution = approximate(B = 10000), alternative = "less")
## pci <- attr(pvalue(ctxo), "conf.int")
## stopifnot(pci[1] < 0.0634 & pci[2] > 0.0634) # p-value

## ## Two-sided approximative, p. 195
## ctx <- conover_test(ftime ~ tyretype, data = failure, ties.method = "average",
##                     distribution = approximate(B = 10000))
## pci <- attr(pvalue(ctx), "conf.int")
## stopifnot(pci[1] < 0.1300 & pci[2] > 0.1300) # p-value

## ## Clean-up
## rm(cer_pkg, ctao, cta, cteo, cte, ctxo, ctx, pci)
## ################################################################################


################################################################################
### 7.13    Permutation Test with General Scores
################################################################################

################################################################################
### 7.13.2  Example with Unstratified Data: Diastolic Blod Pressure, p.198
load("BLOODPR.rda")

##
## 'oneway_test.formula':
##

## One-sided asymptotic, p. 199
otao <- oneway_test(response ~ group, data = bloodpr,
                    alternative = "greater")
stopifnot(isequal(round(statistic(otao), 4), 1.6117)) # test statistic
stopifnot(isequal(round(pvalue(otao), 4), 0.0535)) # p-value

## Two-sided asymptotic, p. 199
ota <- oneway_test(response ~ group, data = bloodpr)
stopifnot(isequal(round(pvalue(ota), 4), 0.1070)) # p-value

## One-sided exact, p. 199
oteo <- oneway_test(response ~ group, data = bloodpr,
                    distribution = "exact", alternative = "greater")
stopifnot(isequal(round(pvalue(oteo), 4), 0.0564)) # p-value

## Two-sided exact, p. 199
ote <- oneway_test(response ~ group, data = bloodpr,
                   distribution = "exact")
stopifnot(isequal(round(pvalue(ote), 4), 0.1040)) # p-value
stopifnot(isequal(round(dperm(ote, statistic(ote)), 4), 0.0176)) # point prob

## One-sided approximative, p. 199
otxo <- oneway_test(response ~ group, data = bloodpr,
                    distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(otxo), "conf.int")
stopifnot(pci[1] < 0.0564 & pci[2] > 0.0564) # p-value

## Two-sided approximative, p. 199
otx <- oneway_test(response ~ group, data = bloodpr,
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(otx), "conf.int")
stopifnot(pci[1] < 0.1040 & pci[2] > 0.1040) # p-value

## Clean-up
rm(bloodpr, otao, ota, oteo, ote, otxo, otx, pci)
################################################################################

################################################################################
### 7.13.3  Example with Stratified Data: Employment Discrimination, p.199
load("DISCRIM.rda")

##
## 'oneway_test.formula':
##

## One-sided asymptotic, p. 160
otao <- oneway_test(salary ~ gender | year, data = discrim,
                    alternative = "less")
stopifnot(isequal(round(statistic(otao), 4), -1.8734)) # test statistic
stopifnot(isequal(round(pvalue(otao), 4), 0.0305)) # p-value

## Two-sided asymptotic, p. 160
ota <- oneway_test(salary ~ gender | year, data = discrim)
stopifnot(isequal(round(pvalue(ota), 4), 0.0610)) # p-value

## ## One-sided exact, p. 160
## oteo <- oneway_test(salary ~ gender | year, data = discrim,
##                     distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(oteo), 4), 0.0400)) # p-value

## ## Two-sided exact, p. 160
## ote <- oneway_test(salary ~ gender | year, data = discrim,
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(ote), 4), 0.0400)) # p-value
## stopifnot(isequal(round(dperm(ote, statistic(ote)), 4), 0.0400)) # point prob

## One-sided approximative, p. 160
otxo <- oneway_test(salary ~ gender | year, data = discrim,
                    distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(otxo), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Two-sided approximative, p. 160
otx <- oneway_test(salary ~ gender | year, data = discrim,
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(otx), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Clean-up
rm(discrim, otao, ota, #oteo, ote,
   otxo, otx, pci)
################################################################################

################################################################################
### 7.13.4  Example using Stratum Specific Scores, p. 201
load("DISCRIM.rda")

##
## 'oneway_test.formula':
##

## One-sided asymptotic, p. 207
otao <- oneway_test(salary ~ gender | year, data = discrim,
                    ytrafo = function(data)
                        trafo(data, numeric_trafo = rank, block = discrim$year),
                    alternative = "less")
stopifnot(isequal(round(statistic(otao), 3), -2.012)) # test statistic
stopifnot(isequal(round(pvalue(otao), 5), 0.02209)) # p-value

## Two-sided asymptotic, p. 207
ota <- oneway_test(salary ~ gender | year, data = discrim,
                   ytrafo = function(data)
                        trafo(data, numeric_trafo = rank, block = discrim$year))
stopifnot(isequal(round(pvalue(ota), 5), 0.04417)) # p-value

## ## One-sided exact, p. 207
## oteo <- oneway_test(salary ~ gender | year, data = discrim,
##                     ytrafo = function(data)
##                         trafo(data, numeric_trafo = rank, block = discrim$year),
##                     distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(oteo), 5), 0.04)) # p-value

## ## Two-sided exact, p. 207
## ote <- oneway_test(salary ~ gender | year, data = discrim,
##                    ytrafo = function(data)
##                         trafo(data, numeric_trafo = rank, block = discrim$year),
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(ote), 5), 0.06667)) # p-value
## stopifnot(isequal(round(dperm(ote, statistic(ote)), 4), 0.0400)) # point prob

## One-sided approximative, p. 207
otxo <- oneway_test(salary ~ gender | year, data = discrim,
                    ytrafo = function(data)
                        trafo(data, numeric_trafo = rank, block = discrim$year),
                    distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(otxo), "conf.int")
stopifnot(pci[1] < 0.0400 & pci[2] > 0.0400) # p-value

## Two-sided approximative, p. 207
otx <- oneway_test(salary ~ gender | year, data = discrim,
                   ytrafo = function(data)
                        trafo(data, numeric_trafo = rank, block = discrim$year),
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(otx), "conf.int")
stopifnot(pci[1] < 0.06667 & pci[2] > 0.06667) # p-value

## Clean-up
rm(discrim, otao, ota, #oteo, ote,
   otxo, otx, pci)
################################################################################

################################################################################
### 7.13.5  Example Using MERT Scores, p. 207
################################################################################


################################################################################
### 7.14    Logrank Test
################################################################################

################################################################################
### 7.14.3  Example with Unstratified Data: Lung Cancer Study, p. 213
load("CANCER.rda")

##
## 'surv_test.formula':
##

## One-sided asymptotic, p. 215
stao <- surv_test(Surv(response, censored) ~ drug, data = cancer,
                  ties.method = "average",
                  alternative = "greater")
stopifnot(isequal(round(statistic(stao), 4), 2.9492)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0016)) # p-value

## Two-sided asymptotic, p. 215
sta <- surv_test(Surv(response, censored) ~ drug, data = cancer,
                 ties.method = "average")
stopifnot(isequal(round(pvalue(sta), 4), 0.0032)) # p-value

## One-sided exact, p. 215
steo <- surv_test(Surv(response, censored) ~ drug, data = cancer,
                  ties.method = "average",
                  distribution = "exact", alternative = "greater")
stopifnot(isequal(round(pvalue(steo), 4), 0.0010)) # p-value

## Two-sided exact, p. 215
ste <- surv_test(Surv(response, censored) ~ drug, data = cancer,
                 ties.method = "average",
                 distribution = "exact")
stopifnot(isequal(round(pvalue(ste), 4), 0.0010)) # p-value
stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0005)) # point prob
## <FIXME>  Why is this???  </FIXME>

## One-sided approximative, p. 215
stxo <- surv_test(Surv(response, censored) ~ drug, data = cancer,
                  ties.method = "average",
                  distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0010 & pci[2] > 0.0010) # p-value

## Two-sided approximative, p. 215
stx <- surv_test(Surv(response, censored) ~ drug, data = cancer,
                 ties.method = "average",
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0010 & pci[2] > 0.0010) # p-value

## Clean-up
rm(cancer, stao, sta, steo, ste, stxo, stx, pci)
################################################################################

################################################################################
### 7.14.4  Example with Stratified Data: Survival by Treatment and Gender, p.218
load("SUR-TRT.rda")

##
## 'surv_test.formula':
##

## One-sided asymptotic, p. 219
stao <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
                  ties.method = "average",
                  alternative = "less")
stopifnot(isequal(round(statistic(stao), 4), -1.5446)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0612)) # p-value

## Two-sided asymptotic, p. 219
sta <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
                 ties.method = "average")
stopifnot(isequal(round(pvalue(sta), 4), 0.1224)) # p-value

## ## One-sided exact, p. 219
## steo <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
##                   ties.method = "average",
##                   distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(steo), 4), 0.0900)) # p-value

## ## Two-sided exact, p. 219
## ste <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
##                  ties.method = "average",
##                  distribution = "exact")
## stopifnot(isequal(round(pvalue(ste), 4), 0.1600)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0100)) # point prob

## One-sided approximative, p. 219
stxo <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
                  ties.method = "average",
                  distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0900 & pci[2] > 0.0900) # p-value

## Two-sided approximative, p. 219
stx <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
                 ties.method = "average",
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.1600 & pci[2] > 0.1600) # p-value

## Clean-up
rm(sur_trt, stao, sta, #steo, ste,
   stxo, stx, pci)
################################################################################


## ################################################################################
## ### 7.15    Generalized Wilcoxon-Gehan Test
## ################################################################################

## ################################################################################
## ### 7.15.3  Example with Unstratified Data: Lung Cancer Study, p. 221
## load("CANCER.rda")

## ##
## ## 'surv_test.formula':
## ##

## ## One-sided asymptotic, p. 222
## stao <- surv_test(Surv(response, censored) ~ drug, data = cancer,
##                   ties.method = "average", type = "Prentice",
##                   alternative = "greater")
## stopifnot(isequal(round(statistic(stao), 4), 2.7813)) # test statistic
## stopifnot(isequal(round(pvalue(stao), 4), 0.0027)) # p-value

## ## Two-sided asymptotic, p. 222
## sta <- surv_test(Surv(response, censored) ~ drug, data = cancer,
##                  ties.method = "average", type = "Prentice")
## stopifnot(isequal(round(pvalue(sta), 4), 0.0054)) # p-value

## ## One-sided exact, p. 222
## steo <- surv_test(Surv(response, censored) ~ drug, data = cancer,
##                   ties.method = "average", type = "Prentice",
##                   distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(steo), 4), 0.0015)) # p-value

## ## Two-sided exact, p. 222
## ste <- surv_test(Surv(response, censored) ~ drug, data = cancer,
##                  ties.method = "average", type = "Prentice",
##                  distribution = "exact")
## stopifnot(isequal(round(pvalue(ste), 4), 0.0030)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0005)) # point prob
## ## <FIXME>  Why is this???  </FIXME>

## ## One-sided approximative, p. 222
## stxo <- surv_test(Surv(response, censored) ~ drug, data = cancer,
##                   ties.method = "average", type = "Prentice",
##                   distribution = approximate(B = 10000), alternative = "greater")
## pci <- attr(pvalue(stxo), "conf.int")
## stopifnot(pci[1] < 0.0015 & pci[2] > 0.0015) # p-value

## ## Two-sided approximative, p. 222
## stx <- surv_test(Surv(response, censored) ~ drug, data = cancer,
##                  ties.method = "average", type = "Prentice",
##                  distribution = approximate(B = 10000))
## pci <- attr(pvalue(stx), "conf.int")
## stopifnot(pci[1] < 0.0030 & pci[2] > 0.0030) # p-value

## ## Clean-up
## rm(cancer, stao, sta, steo, ste, stxo, stx, pci)
## ################################################################################

## ################################################################################
## ### 7.15.4  Example with Stratified Data: Survival by Treatment and Gender, p.223
## load("SUR-TRT.rda")

## ##
## ## 'surv_test.formula':
## ##

## ## One-sided asymptotic, p. 224
## stao <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
##                   ties.method = "average", type = "Prentice",
##                   alternative = "less")
## stopifnot(isequal(round(statistic(stao), 4), -1.7765)) # test statistic
## stopifnot(isequal(round(pvalue(stao), 4), 0.0378)) # p-value

## ## Two-sided asymptotic, p. 224
## sta <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
##                  ties.method = "average", type = "Prentice")
## stopifnot(isequal(round(pvalue(sta), 4), 0.0756)) # p-value

## ## ## One-sided exact, p. 224
## ## steo <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
## ##                  ties.method = "average", type = "Prentice",
## ##                   distribution = "exact", alternative = "less")
## ## stopifnot(isequal(round(pvalue(steo), 4), 0.0700)) # p-value

## ## ## Two-sided exact, p. 224
## ## ste <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
## ##                  ties.method = "average", type = "Prentice",
## ##                  distribution = "exact")
## ## stopifnot(isequal(round(pvalue(ste), 4), 0.1000)) # p-value
## ## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0100)) # point prob

## ## One-sided approximative, p. 224
## stxo <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
##                   ties.method = "average", type = "Prentice",
##                   distribution = approximate(B = 10000), alternative = "less")
## pci <- attr(pvalue(stxo), "conf.int")
## stopifnot(pci[1] < 0.0700 & pci[2] > 0.0700) # p-value

## ## Two-sided approximative, p. 224
## stx <- surv_test(Surv(response, censor) ~ trtmnt | gender, data = sur_trt,
##                  ties.method = "average", type = "Prentice",
##                  distribution = approximate(B = 10000))
## pci <- attr(pvalue(stx), "conf.int")
## stopifnot(pci[1] < 0.1000 & pci[2] > 0.1000) # p-value

## ## Clean-up
## rm(sur_trt, stao, sta, #steo, ste,
##    stxo, stx, pci)
## ################################################################################


################################################################################
### 8      K-sample Inference: Related (Blocked) Samples
################################################################################

################################################################################
### 8.4    Friedman Test
################################################################################

################################################################################
### 8.4.1  Example 1: Effect of Hypnosis on Skin Potential, p. 245
load("HYPNO.rda")

##
## 'friedman_test.formula':
##

## Two-sided asymptotic, p. 246
fta <- friedman_test(potential ~ treatment | subject, data = hypno)
stopifnot(isequal(round(statistic(fta), 4), 9.1525)) # test statistic
stopifnot(isequal(round(pvalue(fta), 4), 0.0574)) # p-value

## ## Two-sided exact, p. 246
## fte <- friedman_test(potential ~ treatment | subject, data = hypno,
##                      distribution = "exact")
## stopifnot(isequal(round(pvalue(fte), 4), 0.0268)) # p-value
## stopifnot(isequal(round(dperm(fte, statistic(fte)), 4), 0.0025)) # point prob

## Two-sided approximative, p. 246 (247)
ftx <- friedman_test(potential ~ treatment | subject, data = hypno,
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(ftx), "conf.int")
stopifnot(pci[1] < 0.0268 & pci[2] > 0.0268) # p-value

## Clean-up
rm(hypno, fta, #fte,
   ftx, pci)
################################################################################

################################################################################
### 8.4.2  Example 2: Relationsship of Friedman's Tests to the Sign, p. 247
### <FIXME>  Rename to azt1, or ???
load("AIDS.rda")
### </FIXME>

diff <- with(AIDS[1:8, ], post - pre)
diff0 <- diff[abs(diff) > 0]
y <- as.vector(t(cbind(abs(diff0) * (diff0 < 0), abs(diff0) * (diff0 >= 0))))
x <- factor(rep(c("neg", "pos"), length(diff0)), levels = c("pos", "neg"))
b <- gl(length(diff0), 2)

##
## 'friedman_test.formula':
##

## Two-sided asymptotic, p. 248
fta <- friedman_test(y ~ x | b)
stopifnot(isequal(round(statistic(fta), 4), 0.6667)) # test statistic
stopifnot(isequal(round(pvalue(fta), 4), 0.4142)) # p-value

## ## Two-sided exact, p. 248
## fte <- friedman_test(y ~ x | b,
##                      distribution = "exact")
## stopifnot(isequal(round(pvalue(fte), 4), 0.6875)) # p-value
## stopifnot(isequal(round(dperm(fte, statistic(fte)), 4), 0.4688)) # point prob

## Two-sided approximative, p. 248
ftx <- friedman_test(y ~ x | b,
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(ftx), "conf.int")
stopifnot(pci[1] < 0.6875 & pci[2] > 0.6875) # p-value

## Clean-up
rm(AIDS, diff, diff0, y, x, b, fta, #fte,
   ftx, pci)
################################################################################


################################################################################
### 8.5    Kendall's W or Coefficient of Concordance
################################################################################

################################################################################
### 8.5.1  Example 1: Choice of Site for an Annual Professional Meeting, p. 249
load("MEET.rda")

##
## 'friedman_test.formula':
##

## Two-sided asymptotic, p. 251
fta <- friedman_test(rating ~ factor | rater, data = meet)
stopifnot(isequal(round(statistic(fta), 4), 13.7778)) # test statistic
stopifnot(isequal(round(pvalue(fta), 4), 0.0553)) # p-value

## Clean-up
rm(meet, fta)
################################################################################


################################################################################
### 8.6    Cochran's Q Test
################################################################################

################################################################################
### 8.6.1  Example 1: Cross-Over Clinical Trial of Analgesic Efficacy, p. 253
load("CLNTRL.rda")
dta <- droplevels(subset(clntrl, treatment != "Aspirin"))

##
## 'mh_test.formula':
##

## Two-sided asymptotic, p. 254
mta <- mh_test(outcome ~ treatment | subject, data = dta)
stopifnot(isequal(round(statistic(mta), 4), 7)) # test statistic
stopifnot(isequal(round(pvalue(mta), 4), 0.0082)) # p-value

## ## Two-sided exact, p. 254
## mte <- mh_test(outcome ~ treatment | subject, data = dta,
##                distribution = "exact")
## stopifnot(isequal(round(pvalue(mte), 4), 0.0156)) # p-value
## stopifnot(isequal(round(dperm(mte, statistic(mte)), 4), 0.0156)) # point prob

## Two-sided approximative, p. 254
mtx <- mh_test(outcome ~ treatment | subject, data = dta,
               distribution = approximate(B = 10000))
pci <- attr(pvalue(mtx), "conf.int")
stopifnot(pci[1] < 0.0156 & pci[2] > 0.0156) # p-value

## Clean-up
rm(clntrl, dta, mta, #mte,
   mtx, pci)
################################################################################


################################################################################
### 8.7    Friedman Aligned Rank Test
################################################################################

################################################################################
### 8.7.1  Example 3: Effect of Hypnosis on Skin Potential, p. 256
load("HYPNO.rda")

align_trafo <-
    function(x, numeric_trafo = id_trafo, block = NULL, align.fun, ...) {
        x[] <- trafo(x, numeric_trafo = function(x) align.fun(x, ...),
                     block = block)
        trafo(x, numeric_trafo = numeric_trafo)
    }

##
## 'oneway_test.formula':
##

## Two-sided asymptotic, p. 257
ota <- oneway_test(potential ~ treatment | subject, data = hypno,
                   ytrafo = function(data)
                       align_trafo(data, numeric_trafo = rank,
                                   block = hypno$subject,
                                   align.fun = function(y) y - mean(y)),
                   teststat = "quad")
stopifnot(isequal(round(statistic(ota), 1), 9.1)) # test statistic
stopifnot(isequal(round(pvalue(ota), 5), 0.05865)) # p-value

## ## Two-sided exact, p. 257
## ote <- oneway_test(potential ~ treatment | subject, data = hypno,
##                    ytrafo = function(data)
##                        align_trafo(data, numeric_trafo = rank,
##                                    block = hypno$subject,
##                                    align.fun = function(y) y - mean(y)),
##                    teststat = "quad",
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(ote), 4), 0.02306)) # p-value
## stopifnot(isequal(round(dperm(ote, statistic(ote)), 4), 0.0004167)) # point prob

## Two-sided approximative, p. 257 (247)
otx <- oneway_test(potential ~ treatment | subject, data = hypno,
                   ytrafo = function(data)
                       align_trafo(data, numeric_trafo = rank,
                                   block = hypno$subject,
                                   align.fun = function(y) y - mean(y)),
                   teststat = "quad",
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(otx), "conf.int")
stopifnot(pci[1] < 0.02306 & pci[2] > 0.02306) # p-value

## Clean-up
rm(hypno, align_trafo, ota, #ote,
   otx, pci)
################################################################################

################################################################################
### 8.7.2  Example 4: Newborn Behaviour Levels, p. 257
load("BEHAVIOUR.rda")

align_trafo <-
    function(x, numeric_trafo = id_trafo, block = NULL, align.fun, ...) {
        x[] <- trafo(x, numeric_trafo = function(x) align.fun(x, ...),
                     block = block)
        trafo(x, numeric_trafo = numeric_trafo)
    }

##
## 'oneway_test.formula':
##

## Two-sided asymptotic, p. 257
ota <- oneway_test(score ~ condition | subject, data = behaviour,
                   ytrafo = function(data)
                       align_trafo(data, numeric_trafo = rank,
                                   block = behaviour$subject,
                                   align.fun = function(y) y - mean(y)),
                   teststat = "quad")
stopifnot(isequal(round(statistic(ota), 3), 2.399)) # test statistic
stopifnot(isequal(round(pvalue(ota), 4), 0.4938)) # p-value
## <FIXME>  Why is this???  </FIXME>

## Clean-up
rm(behaviour, align_trafo, ota)
################################################################################


## ################################################################################
## ### 8.8    Quade Test
## ################################################################################

## ################################################################################
## ### 8.8.1  Example: Effect of Hypnosis on Skin Potential, p. 259
## load("HYPNO.rda")

## ##
## ## 'quade_test.formula':
## ##

## ## Two-sided asymptotic, p. 260
## qta <- quade_test(potential ~ treatment | subject, data = hypno)
## stopifnot(isequal(round(statistic(qta), 4), 4.7500)) # test statistic
## stopifnot(isequal(round(pvalue(qta), 4), 0.0294)) # p-value

## ## ## Two-sided exact, p. 260
## ## qte <- quade_test(potential ~ treatment | subject, data = hypno,
## ##                      distribution = "exact")
## ## stopifnot(isequal(round(pvalue(qte), 4), 0.0171)) # p-value
## ## stopifnot(isequal(round(dperm(qte, statistic(qte)), 4), 0.0022)) # point prob

## ## Two-sided approximative, p. 260 (261)
## qtx <- quade_test(potential ~ treatment | subject, data = hypno,
##                      distribution = approximate(B = 10000))
## pci <- attr(pvalue(qtx), "conf.int")
## stopifnot(pci[1] < 0.0171 & pci[2] > 0.0171) # p-value

## ## Clean-up
## rm(hypno, qta, #qte,
##    qtx, pci)
## ################################################################################


################################################################################
### 8.9    Page Test
################################################################################

################################################################################
### 8.9.1  Example: Breaking Strength of Cotton Fibers, p. 263
load("COTTON.rda")

##
## 'symmetry_test.formula':
##

## One-sided asymptotic, p. 264
stao <- symmetry_test(strength ~ potash | replication, data = cotton,
                      ytrafo = function(data)
                          trafo(data, numeric_trafo = rank,
                                block = cotton$replication),
                      alternative = "greater")
stopifnot(isequal(round(statistic(stao), 3), 2.656)) # test statistic
stopifnot(isequal(round(pvalue(stao), 6), 0.003956)) # p-value

## Two-sided asymptotic, p. 264
sta <- symmetry_test(strength ~ potash | replication, data = cotton,
                     ytrafo = function(data)
                         trafo(data, numeric_trafo = rank,
                               block = cotton$replication))
stopifnot(isequal(round(pvalue(sta), 6), 0.007912)) # p-value

## ## One-sided exact, p. 264
## steo <- symmetry_test(strength ~ potash | replication, data = cotton,
##                       ytrafo = function(data)
##                       trafo(data, numeric_trafo = rank,
##                             block = cotton$replication),
##                       distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(steo), 6), 0.002492)) # p-value

## ## Two-sided exact, p. 264
## ste <- symmetry_test(strength ~ potash | replication, data = cotton,
##                      ytrafo = function(data)
##                          trafo(data, numeric_trafo = rank,
##                                block = cotton$replication),
##                      distribution = "exact")
## stopifnot(isequal(round(pvalue(ste), 6), 0.004985)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 6), 0.001083)) # point prob

## One-sided approximative, p. 264
stxo <- symmetry_test(strength ~ potash | replication, data = cotton,
                     ytrafo = function(data)
                         trafo(data, numeric_trafo = rank,
                               block = cotton$replication),
                      distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.002492 & pci[2] > 0.002492) # p-value

## Two-sided approximative, p. 264
stx <- symmetry_test(strength ~ potash | replication, data = cotton,
                     ytrafo = function(data)
                         trafo(data, numeric_trafo = rank,
                               block = cotton$replication),
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.004985 & pci[2] > 0.004985) # p-value

##
## 'friedman_test.formula':
##

## Two-sided asymptotic, p. 264
fta <- friedman_test(strength ~ potash | replication, data = cotton)
stopifnot(isequal(sqrt(statistic(fta)), unname(statistic(sta)))) # test statistic

## Two-sided approximative, p. 264
ftx <- friedman_test(strength ~ potash | replication, data = cotton,
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.004985 & pci[2] > 0.004985) # p-value

## Clean-up
rm(cotton, stao, sta, #steo, ste,
   stxo, stx, fta, ftx, pci)
################################################################################


################################################################################
### 9      K-sample Inference: Independent Samples
################################################################################

################################################################################
### 9.4    Median Test
################################################################################

################################################################################
### 9.4.1  Example: Hematologic Toxicity Data, p. 278
load("TOXIC.rda")

##
## 'oneway_test.formula':
##

## Two-sided asymptotic, p. 280, 282
ota <- oneway_test(days ~ drug, data = toxic,
                   ytrafo = function(data)
                       trafo(data, numeric_trafo = median_trafo),
                   teststat = "quad")
stopifnot(isequal(round(statistic(ota), 4), 4.1625)) # test statistic
stopifnot(isequal(round(pvalue(ota), 4), 0.3845)) # p-value

## ## Two-sided exact, p. 280, 282
## ote <- oneway_test(days ~ drug, data = toxic,
##                    ytrafo = function(data)
##                        trafo(data, numeric_trafo = median_trafo),
##                    teststat = "quad",
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(ote), 4), 0.4289)) # p-value
## stopifnot(isequal(round(dperm(ote, statistic(ote)), 4), 0.0373)) # point prob

## Two-sided approximative, p. 280, 282
set.seed(711109)
otx <- oneway_test(days ~ drug, data = toxic,
                   ytrafo = function(data)
                       trafo(data, numeric_trafo = median_trafo),
                   teststat = "quad",
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(otx), "conf.int")
stopifnot(pci[1] < 0.4289 & pci[2] > 0.4289) # p-value

##
## 'chisq_test.formula':
##

## StatXact, however, uses the Pearson statistic...
cdays <- factor(toxic$days <= 7)

## Two-sided asymptotic, p. 280, 282
cta <- chisq_test(cdays ~ drug, data = toxic)
stopifnot(isequal(round(statistic(cta), 4), 4.3167)) # test statistic
stopifnot(isequal(round(pvalue(cta), 4), 0.3648)) # p-value

## ## Two-sided exact, p. 280, 282
## cte <- chisq_test(cdays ~ drug, data = toxic,
##                   ytrafo = function(data)
##                       trafo(data, numeric_trafo = median_trafo),
##                   distribution = "exact")
## stopifnot(isequal(round(pvalue(cte), 4), 0.4289)) # p-value
## stopifnot(isequal(round(dperm(cte, statistic(ote)), 4), 0.0373)) # point prob

## Two-sided approximative, p. 280, 282
set.seed(711109)
ctx <- chisq_test(cdays ~ drug, data = toxic,
                  ytrafo = function(data)
                      trafo(data, numeric_trafo = median_trafo),
                  distribution = approximate(B = 10000))
pci <- attr(pvalue(ctx), "conf.int")
stopifnot(pci[1] < 0.4289 & pci[2] > 0.4289) # p-value

## Extra check
stopifnot(isequal(pvalue(ctx), pvalue(otx)))

## Clean-up
rm(toxic, ota, #ote,
   otx, cdays, cta, #cte,
   ctx, pci)
################################################################################


################################################################################
### 9.5    Kruskal-Wallis Test
################################################################################

################################################################################
### 9.5.1  Example: Hematologic Toxicity Data, p. 284
load("TOXIC.rda")

##
## 'kruskal_test.formula':
##

## Two-sided asymptotic, p. 285
kta <- kruskal_test(days ~ drug, data = toxic)
stopifnot(isequal(round(statistic(kta), 4), 9.4147)) # test statistic
stopifnot(isequal(round(pvalue(kta), 4), 0.0515)) # p-value

## Clean-up
rm(toxic, kta)
################################################################################


################################################################################
### 9.6    Normal Scores Test
################################################################################

################################################################################
### 9.6.1  Example: Hematologic Toxicity Data, p. 287
load("TOXIC.rda")

##
## 'oneway_test.formula':
##

## Two-sided asymptotic, p. 288
ota <- oneway_test(days ~ drug, data = toxic,
                   ytrafo = function(data)
                       trafo(data, numeric_trafo = function(y)
                           normal_trafo(y, ties = "average")),
                   teststat = "quad")
stopifnot(isequal(round(statistic(ota), 4), 9.8562)) # test statistic
stopifnot(isequal(round(pvalue(ota), 4), 0.0429)) # p-value

## Clean-up
rm(toxic, ota)
################################################################################


################################################################################
### 9.7    Savage Scores Test
################################################################################

################################################################################
### 9.7.1  Example: Hematologic Toxicity Data, p. 290
load("TOXIC.rda")

##
## 'surv_test.formula':
##

## Two-sided asymptotic, p. 291
sta <- surv_test(Surv(days) ~ drug, data = toxic, ties.method = "average")
stopifnot(isequal(round(statistic(sta), 4), 9.0699)) # test statistic
stopifnot(isequal(round(pvalue(sta), 4), 0.0594)) # p-value

## Clean-up
rm(toxic, sta)
################################################################################


################################################################################
### 9.8    Permutation One-Way ANOVA with General Scores
################################################################################

################################################################################
### 9.8.1  Example: Kruskal-Wallis Testa, p. 293
load("TOXIC.rda")

##
## 'oneway_test.formula':
##

## Two-sided asymptotic, p. 296
ota <- oneway_test(days ~ drug, data = toxic,
                   ytrafo = function(data)
                       trafo(data, numeric_trafo = rank),
                   teststat = "quad")
stopifnot(isequal(round(statistic(ota), 4), 9.4147)) # test statistic
stopifnot(isequal(round(pvalue(ota), 4), 0.0515)) # p-value

## Clean-up
rm(toxic, ota)
################################################################################

################################################################################
### 9.8.2  Example: ANOVA on the Raw Data, p. 297
load("TOXIC.rda")

##
## 'oneway_test.formula':
##

## Two-sided asymptotic, p. 298
ota <- oneway_test(days ~ drug, data = toxic, teststat = "quad")
stopifnot(isequal(round(statistic(ota), 4), 10.9755)) # test statistic
stopifnot(isequal(round(pvalue(ota), 4), 0.0268)) # p-value

## Clean-up
rm(toxic, ota)
################################################################################

################################################################################
### 9.8.3  Example: ANOVA on Censored Survival Data, p. 298
load("BRAIN.rda")

## ##
## ## 'oneway_test.formula':
## ##

## ## Two-sided asymptotic, p. 302
## ota <- oneway_test(Surv(survival, censor) ~ trtmnt, data = brain,
##                    ytrafo = function(data)
##                        trafo(data, surv_trafo = function(y)
##                            logrank_trafo(y, ties.method = "average",
##                                          type = "Prentice")),
##                    teststat = "quad")
## stopifnot(isequal(round(statistic(ota), 4), 4.9819)) # test statistic
## stopifnot(isequal(round(pvalue(ota), 4), 0.1731)) # p-value

##
## 'oneway_test.formula':
##

## Two-sided asymptotic, p. 305
ota <- oneway_test(Surv(survival, censor) ~ trtmnt, data = brain,
                   ytrafo = function(data)
                       trafo(data, surv_trafo = function(y)
                           logrank_trafo(y, ties.method = "average")),
                   teststat = "quad")
stopifnot(isequal(round(statistic(ota), 4), 5.0121)) # test statistic
stopifnot(isequal(round(pvalue(ota), 4), 0.1709)) # p-value
## <FIXME>  Why is this???  </FIXME>

## Clean-up
rm(brain, ota)
################################################################################


################################################################################
### 9.10    Linear-by-Linear Association Test
################################################################################

################################################################################
### 9.10.1  Example: The Space-Shutttle O-Ring Incidents, p. 311
load("SPACE.rda")

##
## 'oneway_test.formula':
##

## One-sided asymptotic, p. 311
otao <- oneway_test(temp ~ oring, data = space,
                    alternative = "less")
stopifnot(isequal(round(statistic(otao), 4), -2.6984)) # test statistic
stopifnot(isequal(round(pvalue(otao), 4), 0.0035)) # p-value

## Two-sided asymptotic, p. 311
ota <- oneway_test(temp ~ oring, data = space)
stopifnot(isequal(round(pvalue(ota), 4), 0.0070)) # p-value

## ## One-sided exact, p. 311
## oteo <- oneway_test(temp ~ oring, data = space,
##                     distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(oteo), 4), 0.0047)) # p-value

## ## Two-sided exact, p. 311
## ote <- oneway_test(temp ~ oring, data = space,
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(ote), 4), 0.0050)) # p-value
## stopifnot(isequal(round(dperm(ote, statistic(ote)), 6), 0.0006)) # point prob

## One-sided approximative, p. 311
otxo <- oneway_test(temp ~ oring, data = space,
                    distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(otxo), "conf.int")
stopifnot(pci[1] < 0.0047 & pci[2] > 0.0047) # p-value

## Two-sided approximative, p. 311
otx <- oneway_test(temp ~ oring, data = space,
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(otx), "conf.int")
stopifnot(pci[1] < 0.0050 & pci[2] > 0.0050) # p-value

## Clean-up
rm(space, otao, ota, #oteo, ote,
   otxo, otx, pci)
################################################################################


################################################################################
### 9.11    K-Sample Logrank Test
################################################################################

################################################################################
### 9.11.1  Example: Brain Tumor Data, p. 313
load("BRAIN.rda")

##
## 'surv_test.formula':
##

## Two-sided asymptotic, p. 314
sta <- surv_test(Surv(survival, censor) ~ trtmnt, data = brain,
                 ties.method = "average")
stopifnot(isequal(round(statistic(sta), 4), 5.0121)) # test statistic
stopifnot(isequal(round(pvalue(sta), 4), 0.1709)) # p-value
## <FIXME>  Why is this???  </FIXME>

## Clean-up
rm(brain, sta)
################################################################################


## ################################################################################
## ### 9.12    K-Sample Wilcoxon-Gehan Test
## ################################################################################

## ################################################################################
## ### 9.12.1  Example: Brain Tumor Data, p. 316
## load("BRAIN.rda")

## ##
## ## 'surv_test.formula':
## ##

## ## Two-sided asymptotic, p. 317
## sta <- surv_test(Surv(survival, censor) ~ trtmnt, data = brain,
##                  ties.method = "average", type = "Prentice")
## stopifnot(isequal(round(statistic(sta), 4), 4.9819)) # test statistic
## stopifnot(isequal(round(pvalue(sta), 4), 0.1731)) # p-value

## ## Clean-up
## rm(brain, sta)
## ################################################################################


################################################################################
### 9.13    Test for Trend with Censored Survival Data
################################################################################

################################################################################
### 9.13.1  Example: Censored Trend Test using Linear-by-Linear Option, p. 318
load("BRAIN.rda")

brain$trtmnt <- ordered(brain$trtmnt)

##
## 'oneway_test.formula':
##

## One-sided asymptotic, p. 320
otao <- oneway_test(Surv(survival, censor) ~ trtmnt, data = brain,
                    ytrafo = function(data)
                        trafo(data, surv_trafo = function(y)
                            logrank_trafo(y, ties.method = "average")),
                    alternative = "greater")
stopifnot(isequal(round(statistic(otao), 4), 1.7732)) # test statistic
stopifnot(isequal(round(pvalue(otao), 4), 0.0381)) # p-value
## <FIXME>  Why is this???  </FIXME>

## Two-sided asymptotic, p. 320
ota <- oneway_test(Surv(survival, censor) ~ trtmnt, data = brain,
                   ytrafo = function(data)
                       trafo(data, surv_trafo = function(y)
                           logrank_trafo(y, ties.method = "average")))
stopifnot(isequal(round(pvalue(ota), 4), 0.0762)) # p-value
## <FIXME>  Why is this???  </FIXME>

##
## 'surv_test.formula':
##

## Two-sided asymptotic, p. 320
sta <- surv_test(Surv(survival, censor) ~ trtmnt, data = brain,
                 ties.method = "average")
stopifnot(isequal(sqrt(statistic(sta)), unname(statistic(ota)))) # test statistic

## Clean-up
rm(brain, otao, ota, sta)
################################################################################


################################################################################
### 14      Two Independent Binomial Samples
################################################################################

################################################################################
### 14.4    Examples
################################################################################

################################################################################
### 14.4.1  Conditional Exact Tests, p. 429
load("CLNTRLT.rda")

dta <- table2df(clntrlt)

## Pearson's Chi-square Exact Test: Clinical Trial Data, p. 430

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 431
itao <- independence_test(Outcome ~ Drug, data = dta,
                          alternative = "less")
## Convert to Pearson's statistic
n <- length(itao@statistic@ytrans)
adj <- (n - 1) / n # Pearson
ts <- (statistic(itao, "linear") - expectation(itao)) / sqrt(variance(itao) * adj)
itao@statistic@teststatistic <- as.vector(ts)
stopifnot(isequal(round(statistic(itao)^2, 4), 3.8095)) # test statistic
stopifnot(isequal(round(pvalue(itao), 4), 0.0255)) # p-value

## Two-sided asymptotic, p. 431
ita <- independence_test(Outcome ~ Drug, data = dta)
## Convert to Pearson's statistic
ita@statistic@teststatistic <- as.vector(ts)
stopifnot(isequal(round(pvalue(ita), 4), 0.0510)) # p-value

## One-sided exact, p. 431
iteo <- independence_test(Outcome ~ Drug, data = dta,
                          distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(iteo), 4), 0.0704)) # p-value
stopifnot(isequal(round(dperm(iteo, statistic(iteo)), 4), 0.0650)) # point prob

## Two-sided exact, p. 431
ite <- independence_test(Outcome ~ Drug, data = dta,
                          distribution = "exact")
stopifnot(isequal(round(pvalue(ite), 4), 0.1409)) # p-value

## One-sided approximative, p. 431
itxo <- independence_test(Outcome ~ Drug, data = dta,
                          distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0704 & pci[2] > 0.0704) # p-value

## Two-sided approximative, p. 431
itx <- independence_test(Outcome ~ Drug, data = dta,
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.1409 & pci[2] > 0.1409) # p-value

##
## 'chisq_test.formula':
##

## Two-sided asymptotic, p. 431
cta <- chisq_test(Outcome ~ Drug, data = dta)
stopifnot(isequal(statistic(cta), statistic(ita)^2)) # test statistic

## Two-sided approximative, p. 431
ctx <- chisq_test(Outcome ~ Drug, data = dta,
                  distribution = approximate(B = 10000))
pci <- attr(pvalue(ctx), "conf.int")
stopifnot(pci[1] < 0.1409 & pci[2] > 0.1409) # p-value

##
## 'chisq_test.table':
##

## Two-sided asymptotic, p. 431
cta2 <- chisq_test(clntrlt)
stopifnot(isequal(statistic(cta2), statistic(ita)^2)) # test statistic

## Two-sided approximative, p. 431
ctx2 <- chisq_test(clntrlt,
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(ctx2), "conf.int")
stopifnot(pci[1] < 0.1409 & pci[2] > 0.1409) # p-value

## Clean-up
rm(clntrlt, dta, itao, n, adj, ts, ita, iteo, ite, itxo, itx, cta, ctx, cta2,
   ctx2, pci)
################################################################################


################################################################################
### 17      C Ordered Binomials
################################################################################

################################################################################
### 17.4    Trend Test with Equally Spaced Scores
################################################################################

################################################################################
### 17.4.1  Maternal Drinking Example, p. 538
load("KORN.rda")

dta <- table2df(korn)
dta$Alcohol.Consumption <- ordered(dta$Alcohol.Consumption)
scrs <- 0:4

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 541
itao <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
                          alternative = "less")
stopifnot(isequal(round(statistic(itao), 4), -1.3519)) # test statistic
stopifnot(isequal(round(pvalue(itao), 4), 0.0882)) # p-value

## Two-sided asymptotic, p. 541
ita <- independence_test(Alcohol.Consumption ~ Malformation, data = dta)
stopifnot(isequal(round(pvalue(ita), 4), 0.1764)) # p-value

## ## Works, but is too time consuming...
## ## One-sided exact, p. 541
## iteo <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
##                           distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(iteo), 4), 0.1046)) # p-value

## ## Works, but is too time consuming...
## ## Two-sided exact, p. 541
## ite <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.1790)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 4), 0.0280)) # point prob

## One-sided approximative, p. 541
itxo <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
                          distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.1046 & pci[2] > 0.1046) # p-value

## Two-sided approximative, p. 541
itx <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.1790 & pci[2] > 0.1790) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 541
itao2 <- independence_test(korn,
                           scores = list("Alcohol.Consumption" = scrs),
                           alternative = "less")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 541
ita2 <- independence_test(korn,
                          scores = list("Alcohol.Consumption" = scrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## ## Works, but is too time consuming...
## ## One-sided exact, p. 541
## iteo2 <- independence_test(korn,
##                            scores = list("Alcohol.Consumption" = scrs),
##                            distribution = "exact", alternative = "less")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Works, but is too time consuming...
## ## Two-sided exact, p. 541
## ite2 <- independence_test(korn,
##                           scores = list("Alcohol.Consumption" = scrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value
## stopifnot(isequal(dperm(ite2, statistic(ite2)), dperm(ite, statistic(ite)))) # point prob

## One-sided approximative, p. 541
itxo2 <- independence_test(korn,
                           scores = list("Alcohol.Consumption" = scrs),
                           distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.1046 & pci[2] > 0.1046) # p-value

## Two-sided approximative, p. 541
itx2 <- independence_test(korn,
                          scores = list("Alcohol.Consumption" = scrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.1790 & pci[2] > 0.1790) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 541
lta <- lbl_test(Alcohol.Consumption ~ Malformation, data = dta)
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 541
ltx <- lbl_test(Alcohol.Consumption ~ Malformation, data = dta,
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.1790 & pci[2] > 0.1790) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 541
lta2 <- lbl_test(korn)
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 541
ltx2 <- lbl_test(korn,
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.1790 & pci[2] > 0.1790) # p-value

## Clean-up
rm(korn, dta, scrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)
################################################################################

################################################################################
### 17.4.2  Example 1: Oral Contraceptives: Stratified Binomial Trend Test, p. 541
load("LANCET.rda")

dta <- table2df(lancet)
dta$Cigarettes.Smoked <- ordered(dta$Cigarettes.Smoked)
scrs <- 1:3

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 545
itao <- independence_test(Cigarettes.Smoked ~ Disease.Status | Age.Group,
                          data = dta,
                          alternative = "greater")
stopifnot(isequal(round(statistic(itao), 7), 4.3346371)) # test statistic
stopifnot(isequal(round(pvalue(itao), 7), 0.0000073)) # p-value

## Two-sided asymptotic, p. 545
ita <- independence_test(Cigarettes.Smoked ~ Disease.Status | Age.Group,
                         data = dta)
stopifnot(isequal(round(pvalue(ita), 7), 0.0000146)) # p-value

## ## One-sided exact, p. 545
## iteo <- independence_test(Cigarettes.Smoked ~ Disease.Status | Age.Group,
##                           data = dta,
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 7), 0.0000060)) # p-value

## ## Two-sided exact, p. 545
## ite <- independence_test(Cigarettes.Smoked ~ Disease.Status | Age.Group,
##                          data = dta,
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 7), 0.0000115)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 7), 0.0000045)) # point prob

## One-sided approximative, p. 545
itxo <- independence_test(Cigarettes.Smoked ~ Disease.Status | Age.Group,
                          data = dta,
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0000060 & pci[2] > 0.0000060) # p-value

## Two-sided approximative, p. 545
itx <- independence_test(Cigarettes.Smoked ~ Disease.Status | Age.Group,
                         data = dta,
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0000115 & pci[2] > 0.0000115) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 545
itao2 <- independence_test(lancet, scores = list(Cigarettes.Smoked = scrs),
                           alternative = "greater")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 545
ita2 <- independence_test(lancet, scores = list(Cigarettes.Smoked = scrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## ## One-sided exact, p. 545
## iteo2 <- independence_test(lancet, scores = list(Cigarettes.Smoked = scrs),
##                            distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Two-sided exact, p. 545
## ite2 <- independence_test(lancet, scores = list(Cigarettes.Smoked = scrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value
## stopifnot(isequal(dperm(ite2, statistic(ite2)), dperm(ite, statistic(ite)))) # point prob

## One-sided approximative, p. 545
itxo2 <- independence_test(lancet, scores = list(Cigarettes.Smoked = scrs),
                           distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.0000060 & pci[2] > 0.0000060) # p-value

## Two-sided approximative, p. 545
itx2 <- independence_test(lancet, scores = list(Cigarettes.Smoked = scrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.0000115 & pci[2] > 0.0000115) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 545
lta <- lbl_test(Cigarettes.Smoked ~ Disease.Status | Age.Group, data = dta)
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 545
ltx <- lbl_test(Cigarettes.Smoked ~ Disease.Status | Age.Group, data = dta,
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.0000115 & pci[2] > 0.0000115) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 545
lta2 <- lbl_test(lancet)
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 545
ltx2 <- lbl_test(lancet,
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.0000115 & pci[2] > 0.0000115) # p-value

## Clean-up
rm(lancet, dta, scrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)
################################################################################


################################################################################
### 17.5    Permutation Test with General Scores
################################################################################

################################################################################
### 17.5.1  Example: Maternal Alcohol Example: Unstratified Binomal Trend, p. 548
load("KORN.rda")
dta <- table2df(korn)
scrs <- c(0, 0.5, 1.5, 4, 7)

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 550
itao <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
                          scores = list("Alcohol.Consumption" = scrs),
                          alternative = "less")
stopifnot(isequal(round(statistic(itao), 3), -2.563)) # test statistic
stopifnot(isequal(round(pvalue(itao), 6), 0.005186)) # p-value

## Two-sided asymptotic, p. 550
ita <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
                         scores = list("Alcohol.Consumption" = scrs))
stopifnot(isequal(round(pvalue(ita), 5), 0.01037)) # p-value

## ## Works, but is too time consuming...
## ## One-sided exact, p. 550
## iteo <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
##                           scores = list("Alcohol.Consumption" = scrs),
##                           distribution = "exact", alternative = "less")
## stopifnot(isequal(round(pvalue(iteo), 5), 0.01678)) # p-value

## ## Works, but is too time consuming...
## ## Two-sided exact, p. 550
## ite <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
##                          scores = list("Alcohol.Consumption" = scrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 5), 0.01725)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 5), 0.00291)) # point prob

## One-sided approximative, p. 550
itxo <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
                          scores = list("Alcohol.Consumption" = scrs),
                          distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.01678 & pci[2] > 0.01678) # p-value

## Two-sided approximative, p. 550
itx <- independence_test(Alcohol.Consumption ~ Malformation, data = dta,
                         scores = list("Alcohol.Consumption" = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.01725 & pci[2] > 0.01725) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 550
itao2 <- independence_test(korn,
                           scores = list("Alcohol.Consumption" = scrs),
                           alternative = "less")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 550
ita2 <- independence_test(korn,
                          scores = list("Alcohol.Consumption" = scrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## ## Works, but is too time consuming...
## ## One-sided exact, p. 550
## iteo2 <- independence_test(korn,
##                            scores = list("Alcohol.Consumption" = scrs),
##                            distribution = "exact", alternative = "less")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Works, but is too time consuming...
## ## Two-sided exact, p. 550
## ite2 <- independence_test(korn,
##                           scores = list("Alcohol.Consumption" = scrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value
## stopifnot(isequal(dperm(ite2, statistic(ite2)), dperm(ite, statistic(ite)))) # point prob

## One-sided approximative, p. 550
itxo2 <- independence_test(korn,
                           scores = list("Alcohol.Consumption" = scrs),
                           distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.01678 & pci[2] > 0.01678) # p-value

## Two-sided approximative, p. 550
itx2 <- independence_test(korn,
                          scores = list("Alcohol.Consumption" = scrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.01725 & pci[2] > 0.01725) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 550
lta <- lbl_test(Alcohol.Consumption ~ Malformation, data = dta,
                scores = list("Alcohol.Consumption" = scrs))
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 550
ltx <- lbl_test(Alcohol.Consumption ~ Malformation, data = dta,
                scores = list("Alcohol.Consumption" = scrs),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.01725 & pci[2] > 0.01725) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 550
lta2 <- lbl_test(korn, scores = list("Alcohol.Consumption" = scrs))
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 550
ltx2 <- lbl_test(korn,scores = list("Alcohol.Consumption" = scrs),
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.01725 & pci[2] > 0.01725) # p-value

## Clean-up
rm(korn, dta, scrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)
################################################################################

################################################################################
### 17.5.2  Example 1: Animal Toxicology: Stratified Binomial Trend, p. 550
load("FDA1.rda")
scrs <- c(0, 1, 5, 50)

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 554
itao <- independence_test(dose ~ disease | stratum, data = fda1, weights = ~ freq,
                          scores = list(dose = scrs),
                          alternative = "greater")
stopifnot(isequal(round(statistic(itao), 3), 1.739)) # test statistic
stopifnot(isequal(round(pvalue(itao), 5), 0.04098)) # p-value

## Two-sided asymptotic, p. 554
ita <- independence_test(dose ~ disease | stratum, data = fda1, weights = ~ freq,
                         scores = list(dose = scrs))
stopifnot(isequal(round(pvalue(ita), 5), 0.08196)) # p-value

## ## One-sided exact, p. 554
## iteo <- independence_test(dose ~ disease | stratum, data = fda1, weights = ~ freq,
##                           scores = list(dose = scrs),
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 4), 0.0651)) # p-value

## ## Two-sided exact, p. 554
## ite <- independence_test(dose ~ disease | stratum, data = fda1, weights = ~ freq,
##                          scores = list(dose = scrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 5), 0.07685)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 5), 0.01087)) # point prob

## One-sided approximative, p. 554
itxo <- independence_test(dose ~ disease | stratum, data = fda1, weights = ~ freq,
                          scores = list(dose = scrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0651 & pci[2] > 0.0651) # p-value

## Two-sided approximative, p. 554
itx <- independence_test(dose ~ disease | stratum, data = fda1, weights = ~ freq,
                         scores = list(dose = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.07685 & pci[2] > 0.07685) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 554
lta <- lbl_test(xtabs(freq ~ dose + disease + stratum, data = fda1),
                scores = list(dose = scrs))
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 554
ltx <- lbl_test(xtabs(freq ~ dose + disease + stratum, data = fda1),
                scores = list(dose = scrs),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.07685 & pci[2] > 0.07685) # p-value

## Clean-up
rm(fda1, scrs, itao, ita, #iteo, ite,
   itxo, itx, lta, ltx, pci)
################################################################################

################################################################################
### 17.5.4  Example: Estimating Trend Parameters: Hiroshima Atomic Bomb Surviors, p. 556

### (I)   First strata
load("HIROSH.rda")

dtaI <- droplevels(subset(table2df(hirosh), Age.Group == "0-9"))
hiroshI <- as.table(hirosh[, , 1])
scrs <- c(0, 4.5, 30, 75)

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 560
itao <- independence_test(Radiation.Dose ~ Disease.Status, data = dtaI,
                          scores = list("Radiation.Dose" = scrs),
                          alternative = "greater")
stopifnot(isequal(round(statistic(itao), 3), 1.68)) # test statistic
stopifnot(isequal(round(pvalue(itao), 5), 0.04646)) # p-value

## Two-sided asymptotic, p. 560
ita <- independence_test(Radiation.Dose ~ Disease.Status, data = dtaI,
                         scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(round(pvalue(ita), 5), 0.09291)) # p-value

## One-sided exact, p. 560
iteo <- independence_test(Radiation.Dose ~ Disease.Status, data = dtaI,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = "exact", alternative = "greater")
stopifnot(isequal(round(pvalue(iteo), 5), 0.06529)) # p-value

## Two-sided exact, p. 560
ite <- independence_test(Radiation.Dose ~ Disease.Status, data = dtaI,
                         scores = list("Radiation.Dose" = scrs),
                         distribution = "exact")
stopifnot(isequal(round(pvalue(ite), 5), 0.06825)) # p-value
stopifnot(isequal(round(dperm(ite, statistic(ite)), 6), 0.002704)) # point prob

## One-sided approximative, p. 560
itxo <- independence_test(Radiation.Dose ~ Disease.Status, data = dtaI,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.06529 & pci[2] > 0.06529) # p-value

## Two-sided approximative, p. 560
itx <- independence_test(Radiation.Dose ~ Disease.Status, data = dtaI,
                         scores = list("Radiation.Dose" = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.06825 & pci[2] > 0.06825) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 560
itao2 <- independence_test(hiroshI,
                           scores = list("Radiation.Dose" = scrs),
                           alternative = "greater")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 560
ita2 <- independence_test(hiroshI,
                          scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## ## One-sided exact, p. 560
## iteo2 <- independence_test(hiroshI,
##                            scores = list("Radiation.Dose" = scrs),
##                            distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Two-sided exact, p. 560
## ite2 <- independence_test(hiroshI,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value
## stopifnot(isequal(dperm(ite2, statistic(ite2)), dperm(ite, statistic(ite)))) # point prob

## One-sided approximative, p. 560
itxo2 <- independence_test(hiroshI,
                           scores = list("Radiation.Dose" = scrs),
                           distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.06529 & pci[2] > 0.06529) # p-value

## Two-sided approximative, p. 560
itx2 <- independence_test(hiroshI,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.06825 & pci[2] > 0.06825) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 560
lta <- lbl_test(Radiation.Dose ~ Disease.Status, data = dtaI,
                scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 560
ltx <- lbl_test(Radiation.Dose ~ Disease.Status, data = dtaI,
                scores = list("Radiation.Dose" = scrs),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.06825 & pci[2] > 0.06825) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 560
lta2 <- lbl_test(hiroshI, scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 560
ltx2 <- lbl_test(hiroshI, scores = list("Radiation.Dose" = scrs),
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.06825 & pci[2] > 0.06825) # p-value

## Clean-up
rm(hirosh, dtaI, hiroshI, scrs, itao, ita, iteo, ite, itxo, itx, itao2, ita2,
   #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)

### (II)  Results for All survivors

### (a) Equally Spaced Scores: (0, 1, 2, 3)
load("HIROSH.rda")

dta <- table2df(hirosh)
scrs <- 0:3

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 566
itao <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                          scores = list("Radiation.Dose" = scrs),
                          alternative = "greater")
stopifnot(isequal(round(statistic(itao), 3), 3.391)) # test statistic
stopifnot(isequal(round(pvalue(itao), 7), 0.0003484)) # p-value

## Two-sided asymptotic, p. 566
ita <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                         scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(round(pvalue(ita), 7), 0.0006967)) # p-value

## ## One-sided exact, p. 566
## iteo <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 7), 0.0006466)) # p-value

## ## Two-sided exact, p. 566
## ite <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
##                          scores = list("Radiation.Dose" = scrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 7), 0.0008836)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 7), 0.0002471)) # point prob

## One-sided approximative, p. 566
itxo <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0006466 & pci[2] > 0.0006466) # p-value

## Two-sided approximative, p. 566
itx <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                         scores = list("Radiation.Dose" = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0008836 & pci[2] > 0.0008836) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 566
itao2 <- independence_test(hirosh,
                           scores = list("Radiation.Dose" = scrs),
                           alternative = "greater")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 566
ita2 <- independence_test(hirosh,
                          scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## ## One-sided exact, p. 566
## iteo2 <- independence_test(hirosh,
##                            scores = list("Radiation.Dose" = scrs),
##                            distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Two-sided exact, p. 566
## ite2 <- independence_test(hirosh,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value
## stopifnot(isequal(dperm(ite2, statistic(ite2)), dperm(ite, statistic(ite)))) # point prob

## One-sided approximative, p. 566d
itxo2 <- independence_test(hirosh,
                           scores = list("Radiation.Dose" = scrs),
                           distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.0006466 & pci[2] > 0.0006466) # p-value

## Two-sided approximative, p. 566
itx2 <- independence_test(hirosh,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.0008836 & pci[2] > 0.0008836) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 566
lta <- lbl_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta)
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 566
ltx <- lbl_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0008836 & pci[2] > 0.0008836) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 566
lta2 <- lbl_test(hirosh)
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 566
ltx2<- lbl_test(hirosh,
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.0008836 & pci[2] > 0.0008836) # p-value

## Clean-up
rm(hirosh, dta, scrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)

### (b) Mid-Range Scores: (0, 4.5, 30, 75)
load("HIROSH.rda")

dta <- table2df(hirosh)
scrs <- c(0, 4.5, 30, 75)

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 567
itao <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                          scores = list("Radiation.Dose" = scrs),
                          alternative = "greater")
stopifnot(isequal(round(statistic(itao), 4), 3.1630)) # test statistic
stopifnot(isequal(round(pvalue(itao), 4), 0.0008)) # p-value

## Two-sided asymptotic, p. 567
ita <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                         scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(round(pvalue(ita), 4), 0.0016)) # p-value

## ## One-sided exact, p. 567
## iteo <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 4), 0.0023)) # p-value

## ## Two-sided exact, p. 567
## ite <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
##                          scores = list("Radiation.Dose" = scrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.0023)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 4), 0.0001)) # point prob

## One-sided approximative, p. 567
itxo <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0023 & pci[2] > 0.0023) # p-value

## Two-sided approximative, p. 567
itx <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                         scores = list("Radiation.Dose" = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0023 & pci[2] > 0.0023) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 567
itao2 <- independence_test(hirosh,
                           scores = list("Radiation.Dose" = scrs),
                           alternative = "greater")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 567
ita2 <- independence_test(hirosh,
                          scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## ## One-sided exact, p. 567
## iteo2 <- independence_test(hirosh,
##                            scores = list("Radiation.Dose" = scrs),
##                            distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Two-sided exact, p. 567
## ite2 <- independence_test(hirosh,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value
## stopifnot(isequal(dperm(ite2, statistic(ite2)), dperm(ite, statistic(ite)))) # point prob

## One-sided approximative, p. 567
itxo2 <- independence_test(hirosh,
                           scores = list("Radiation.Dose" = scrs),
                           distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.0023 & pci[2] > 0.0023) # p-value

## Two-sided approximative, p. 567
itx2 <- independence_test(hirosh,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.0023 & pci[2] > 0.0023) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 567
lta <- lbl_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 567
ltx <- lbl_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dta,
                scores = list("Radiation.Dose" = scrs),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.0023 & pci[2] > 0.0023) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 567
lta2 <- lbl_test(hirosh, scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 567
ltx2 <- lbl_test(hirosh, scores = list("Radiation.Dose" = scrs),
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.0023 & pci[2] > 0.0023) # p-value

## Clean-up
rm(hirosh, dta, scrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)

### (III) Results Excluding People Not in the City during Explosion

### (a) Equally Spaced Scores: (0, 1, 2, 3)
load("HIROSH.rda")

dtaIII <- droplevels(subset(table2df(hirosh),
                            Radiation.Dose != "Not in City"))
hiroshIII <- as.table(hirosh[, -1, ])
scrs <- 1:3

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 569
itao <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                          scores = list("Radiation.Dose" = scrs),
                          alternative = "greater")
stopifnot(isequal(round(pvalue(itao), 4), 0.0089)) # p-value

## Two-sided asymptotic, p. 569
ita <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                         scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(round(pvalue(ita), 4), 0.0178)) # p-value

## ## One-sided exact, p. 569
## iteo <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 4), 0.0156)) # p-value

## ## Two-sided exact, p. 569
## ite <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
##                          scores = list("Radiation.Dose" = scrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.0208)) # p-value

## One-sided approximative, p. 569
itxo <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0156 & pci[2] > 0.0156) # p-value

## Two-sided approximative, p. 569
itx <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                         scores = list("Radiation.Dose" = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0208 & pci[2] > 0.0208) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 569
itao2 <- independence_test(hiroshIII,
                           scores = list("Radiation.Dose" = scrs),
                           alternative = "greater")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 569
ita2 <- independence_test(hiroshIII,
                          scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## ## One-sided exact, p. 569
## iteo2 <- independence_test(hiroshIII,
##                            scores = list("Radiation.Dose" = scrs),
##                            distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Two-sided exact, p. 569
## ite2 <- independence_test(hiroshIII,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value

## One-sided approximative, p. 569
itxo2 <- independence_test(hiroshIII,
                           scores = list("Radiation.Dose" = scrs),
                           distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.0156 & pci[2] > 0.0156) # p-value

## Two-sided approximative, p. 569
itx2 <- independence_test(hiroshIII,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.0208 & pci[2] > 0.0208) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 569
lta <- lbl_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII)
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 569
ltx <- lbl_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.0208 & pci[2] > 0.0208) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 569
lta2 <- lbl_test(hiroshIII)
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 569
ltx2 <- lbl_test(hiroshIII,
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.0208 & pci[2] > 0.0208) # p-value

## Clean-up
rm(hirosh, dtaIII, hiroshIII, scrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)

### (b) Mid-Range Scores: (0, 4.5, 30, 75)
load("HIROSH.rda")

dtaIII <- droplevels(subset(table2df(hirosh),
                            Radiation.Dose != "Not in City"))
hiroshIII <- as.table(hirosh[, -1, ])
scrs <- c(4.5, 30, 75)

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 569
itao <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                          scores = list("Radiation.Dose" = scrs),
                          alternative = "greater")
stopifnot(isequal(round(pvalue(itao), 4), 0.0099)) # p-value

## Two-sided asymptotic, p. 569
ita <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                         scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(round(pvalue(ita), 4), 0.0198)) # p-value

## ## One-sided exact, p. 569
## iteo <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 4), 0.0155)) # p-value

## ## Two-sided exact, p. 569
## ite <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
##                          scores = list("Radiation.Dose" = scrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.0200)) # p-value

## One-sided approximative, p. 569
itxo <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0155 & pci[2] > 0.0155) # p-value

## Two-sided approximative, p. 569
itx <- independence_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                         scores = list("Radiation.Dose" = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0200 & pci[2] > 0.0200) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 569
itao2 <- independence_test(hiroshIII,
                           scores = list("Radiation.Dose" = scrs),
                           alternative = "greater")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 569
ita2 <- independence_test(hiroshIII,
                          scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## ## One-sided exact, p. 569
## iteo2 <- independence_test(hiroshIII,
##                            scores = list("Radiation.Dose" = scrs),
##                            distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Two-sided exact, p. 569
## ite2 <- independence_test(hiroshIII,
##                           scores = list("Radiation.Dose" = scrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value

## One-sided approximative, p. 569
itxo2 <- independence_test(hiroshIII,
                           scores = list("Radiation.Dose" = scrs),
                           distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.0155 & pci[2] > 0.0155) # p-value

## Two-sided approximative, p. 569
itx2 <- independence_test(hiroshIII,
                          scores = list("Radiation.Dose" = scrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.0200 & pci[2] > 0.0200) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 569
lta <- lbl_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 569
ltx <- lbl_test(Radiation.Dose ~ Disease.Status | Age.Group, data = dtaIII,
                scores = list("Radiation.Dose" = scrs),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.0200 & pci[2] > 0.0200) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 569
lta2 <- lbl_test(hiroshIII, scores = list("Radiation.Dose" = scrs))
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 569
ltx2 <- lbl_test(hiroshIII, scores = list("Radiation.Dose" = scrs),
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.0200 & pci[2] > 0.0200) # p-value

## Clean-up
rm(hirosh, dtaIII, hiroshIII, scrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)
################################################################################

################################################################################
### 17.5.5  Example: Endometrial Cancer: Marginal Homogeneity Test, p.569
load("LANDIS.rda")

diag(landis) <- 0
scrs <- c(0, 0.2, 0.5125, 0.7)

##
## 'symmetry_test.table':
##

## One-sided asymptotic, p. 573
stao <- symmetry_test(landis, scores = list(response = scrs),
                      alternative = "greater")
stopifnot(isequal(round(statistic(stao), 7), 3.7345758)) # test statistic
stopifnot(isequal(round(pvalue(stao), 7), 0.0000940)) # p-value

## Two-sided asymptotic, p. 573
sta <- symmetry_test(landis, scores = list(response = scrs))
stopifnot(isequal(round(pvalue(sta), 7), 0.0001880)) # p-value

## ## One-sided exact, p. 573
## steo <- symmetry_test(landis, scores = list(response = scrs),
##                       distribution = exact(fact = 1e4), alternative = "greater")
## stopifnot(isequal(round(pvalue(steo), 7), 0.0000940)) # p-value

## ## Two-sided exact, p. 573
## ste <- symmetry_test(landis, scores = list(response = scrs),
##                      distribution = exact(fact = 1e4))
## stopifnot(isequal(round(pvalue(ste), 7), 0.0001880)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 7), 0.0000014)) # point prob

## One-sided approximative, p. 573
stxo <- symmetry_test(landis, scores = list(response = scrs),
                      distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0000940 & pci[2] > 0.0000940) # p-value

## Two-sided approximative, p. 573
stx <- symmetry_test(landis, scores = list(response = scrs),
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0001880 & pci[2] > 0.0001880) # p-value

##
## 'mh_test.table':
##

## Two-sided asymptotic, p. 573
mta <- mh_test(landis, scores = list(response = scrs))
stopifnot(isequal(statistic(mta), unname(statistic(sta)^2))) # test statistic

## Two-sided approximative, p. 573
mtx <- mh_test(landis, scores = list(response = scrs),
               distribution = approximate(B = 10000))
pci <- attr(pvalue(mtx), "conf.int")
stopifnot(pci[1] < 0.0001120 & pci[2] > 0.0001120) # p-value

## Clean-up
rm(landis, scrs, stao, sta, #steo, ste,
   stxo, stx, mta, mtx, pci)
################################################################################


################################################################################
### 18      Two Ordered Multinomials
################################################################################

################################################################################
### 18.4    Wilcoxon Rank Sum Test
################################################################################

################################################################################
### 18.4.1  Example: Childhood Chicken Pox, p. 618
load("VARI.rda")

dta <- table2df(vari)
dta$Adverse.Effect <- as.numeric(dta$Adverse.Effect)

##
## 'wilcox_test.formula':
##

## One-sided asymptotic, p. 620
wtao <- wilcox_test(Adverse.Effect ~ Treatment, data = dta,
                    alternative = "less")
stopifnot(isequal(round(statistic(wtao), 4), -1.6430)) # test statistic
stopifnot(isequal(round(pvalue(wtao), 4), 0.0502)) # p-value

## Two-sided asymptotic, p. 620
wta <- wilcox_test(Adverse.Effect ~ Treatment, data = dta)
stopifnot(isequal(round(pvalue(wta), 4), 0.1004)) # p-value

## One-sided exact, p. 620
wteo <- wilcox_test(Adverse.Effect ~ Treatment, data = dta,
                    distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(wteo), 4), 0.0743)) # p-value

## Two-sided exact, p. 620
wte <- wilcox_test(Adverse.Effect ~ Treatment, data = dta,
                   distribution = "exact")
stopifnot(isequal(round(pvalue(wte), 4), 0.1061)) # p-value
stopifnot(isequal(round(dperm(wte, statistic(wte)), 4), 0.0083)) # point prob

## One-sided approximative, p. 620
wtxo <- wilcox_test(Adverse.Effect ~ Treatment, data = dta,
                    distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(wtxo), "conf.int")
stopifnot(pci[1] < 0.0743 & pci[2] > 0.0743) # p-value

## Two-sided approximative, p. 620
wtx <- wilcox_test(Adverse.Effect ~ Treatment, data = dta,
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(wtx), "conf.int")
stopifnot(pci[1] < 0.1061 & pci[2] > 0.1061) # p-value

## Clean-up
rm(vari, dta, wtao, wta, wteo, wte, wtxo, wtx, pci)
################################################################################

################################################################################
### 18.4.2  Example: Bone Marrow Transplant with Stratified Data, p. 621
load("GVHD.rda")

dta <- table2df(gvhd)
dta$Severity.of.GVHD.Toxicity <- as.numeric(dta$Severity.of.GVHD.Toxicity)

##
## 'wilcox_test.formula':
##

## One-sided asymptotic, p. 623
wtao <- wilcox_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                    data = dta,
                    alternative = "greater")
stopifnot(isequal(round(statistic(wtao), 4), 1.5625)) # test statistic
stopifnot(isequal(round(pvalue(wtao), 4), 0.0591)) # p-value

## Two-sided asymptotic, p. 623
wta <- wilcox_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                   data = dta)
stopifnot(isequal(round(pvalue(wta), 4), 0.1182)) # p-value

## ## One-sided exact, p. 623
## wteo <- wilcox_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
##                     data = dta,
##                     distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(wteo), 4), 0.0639)) # p-value

## ## Two-sided exact, p. 623
## wte <- wilcox_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
##                    data = dta,
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(wte), 4), 0.1281)) # p-value
## stopifnot(isequal(round(dperm(wte, statistic(wte)), 4), 0.0084)) # point prob

## One-sided approximative, p. 623
wtxo <- wilcox_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                    data = dta,
                    distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(wtxo), "conf.int")
stopifnot(pci[1] < 0.0639 & pci[2] > 0.0639) # p-value

## Two-sided approximative, p. 623
wtx <- wilcox_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                   data = dta,
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(wtx), "conf.int")
stopifnot(pci[1] < 0.1281 & pci[2] > 0.1281) # p-value

## Clean-up
rm(gvhd, dta, wtao, wta, #wteo, wte,
   wtxo, wtx, pci)
################################################################################


################################################################################
### 18.5    Normal Scores Test
################################################################################

################################################################################
### 18.5.1  Example: Childhood Chicken Pox, p. 625
load("VARI.rda")

dta <- table2df(vari)
dta$Adverse.Effect <- as.numeric(dta$Adverse.Effect)

##
## 'normal_test.formula':
##

## One-sided asymptotic, p. 627
ntao <- normal_test(Adverse.Effect ~ Treatment, data = dta,
                    ties.method = "average",
                    alternative = "less")
stopifnot(isequal(round(statistic(ntao), 4), -1.6467)) # test statistic
stopifnot(isequal(round(pvalue(ntao), 4), 0.0498)) # p-value

## Two-sided asymptotic, p. 627
nta <- normal_test(Adverse.Effect ~ Treatment, data = dta,
                   ties.method = "average")
stopifnot(isequal(round(pvalue(nta), 4), 0.0996)) # p-value

## One-sided exact, p. 627
nteo <- normal_test(Adverse.Effect ~ Treatment, data = dta,
                    ties.method = "average",
                    distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(nteo), 4), 0.0577)) # p-value

## Two-sided exact, p. 627
nte <- normal_test(Adverse.Effect ~ Treatment, data = dta,
                   ties.method = "average",
                   distribution = "exact")
stopifnot(isequal(round(pvalue(nte), 4), 0.1024)) # p-value
stopifnot(isequal(round(dperm(nte, statistic(nte)), 4), 0.0083)) # point prob
## <FIXME>  Why is this???  </FIXME>

## One-sided approximative, p. 627
ntxo <- normal_test(Adverse.Effect ~ Treatment, data = dta,
                    ties.method = "average",
                    distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(ntxo), "conf.int")
stopifnot(pci[1] < 0.0577 & pci[2] > 0.0577) # p-value

## Two-sided approximative, p. 627
ntx <- normal_test(Adverse.Effect ~ Treatment, data = dta,
                   ties.method = "average",
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(ntx), "conf.int")
stopifnot(pci[1] < 0.1024 & pci[2] > 0.1024) # p-value

## Clean-up
rm(vari, dta, ntao, nta, nteo, nte, ntxo, ntx, pci)
################################################################################

################################################################################
### 18.5.2  Example: Bone Marrow Transplant with Stratified Data, p. 627
load("GVHD.rda")

dta <- table2df(gvhd)
dta$Severity.of.GVHD.Toxicity <- as.numeric(dta$Severity.of.GVHD.Toxicity)

##
## 'normal_test.formula':
##

## One-sided asymptotic, p. 629
ntao <- normal_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                    data = dta, ties.method = "average",
                    alternative = "greater")
stopifnot(isequal(round(statistic(ntao), 4), 1.6322)) # test statistic
stopifnot(isequal(round(pvalue(ntao), 4), 0.0513)) # p-value

## Two-sided asymptotic, p. 629
nta <- normal_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                   data = dta, ties.method = "average")
stopifnot(isequal(round(pvalue(nta), 4), 0.1026)) # p-value

## ## One-sided exact, p. 629
## nteo <- normal_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
##                     data = dta, ties.method = "average",
##                     distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(nteo), 4), 0.0514)) # p-value

## ## Two-sided exact, p. 629
## nte <- normal_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
##                    data = dta, ties.method = "average",
##                    distribution = "exact")
## stopifnot(isequal(round(pvalue(nte), 4), 0.1035)) # p-value
## stopifnot(isequal(round(dperm(nte, statistic(nte)), 4), 0.0009)) # point prob

## One-sided approximative, p. 629
ntxo <- normal_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                    data = dta, ties.method = "average",
                    distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(ntxo), "conf.int")
stopifnot(pci[1] < 0.0514 & pci[2] > 0.0514) # p-value

## Two-sided approximative, p. 629
ntx <- normal_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                   data = dta, ties.method = "average",
                   distribution = approximate(B = 10000))
pci <- attr(pvalue(ntx), "conf.int")
stopifnot(pci[1] < 0.1035 & pci[2] > 0.1035) # p-value

## Clean-up
rm(gvhd, dta, ntao, nta, #nteo, nte,
   ntxo, ntx, pci)
################################################################################


################################################################################
### 18.7    Permutation Test with General Scores
################################################################################

################################################################################
### 18.7.1  Example: Childhood Chicken Pox, p. 631
load("VARI.rda")

dta <- table2df(vari)
scrs <- c(0, 1, 1000, 1005)

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 633
itao <- independence_test(Adverse.Effect ~ Treatment, data = dta,
                          scores = list(Adverse.Effect = scrs),
                          alternative = "less")
stopifnot(isequal(round(statistic(itao), 4), -0.5822)) # test statistic
stopifnot(isequal(round(pvalue(itao), 4), 0.2802)) # p-value

## Two-sided asymptotic, p. 633
ita <- independence_test(Adverse.Effect ~ Treatment, data = dta,
                         scores = list(Adverse.Effect = scrs))
stopifnot(isequal(round(pvalue(ita), 4), 0.5604)) # p-value

## One-sided exact, p. 633
iteo <- independence_test(Adverse.Effect ~ Treatment, data = dta,
                          scores = list(Adverse.Effect = scrs),
                          distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(iteo), 4), 0.1538)) # p-value

## Two-sided exact, p. 633
ite <- independence_test(Adverse.Effect ~ Treatment, data = dta,
                         scores = list(Adverse.Effect = scrs),
                         distribution = "exact")
stopifnot(isequal(round(pvalue(ite), 4), 0.3560)) # p-value
stopifnot(isequal(round(dperm(ite, statistic(ite)), 4), 0.0083)) # point prob

## One-sided approximative, p. 633
itxo <- independence_test(Adverse.Effect ~ Treatment, data = dta,
                          scores = list(Adverse.Effect = scrs),
                          distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.1538 & pci[2] > 0.1538) # p-value

## Two-sided approximative, p. 633
itx <- independence_test(Adverse.Effect ~ Treatment, data = dta,
                         scores = list(Adverse.Effect = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.3560 & pci[2] > 0.3560) # p-value

## Clean-up
rm(vari, dta, scrs, itao, ita, iteo, ite, itxo, itx, pci)
################################################################################

################################################################################
### 18.5.2  Example: Bone Marrow Transplant with Stratified Data, p. 633
load("GVHD.rda")

dta <- table2df(gvhd)
scrs <- c(1, 4, 9, 16, 25)

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 636
itao <- independence_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                          data = dta,
                          scores = list(Severity.of.GVHD.Toxicity = scrs),
                          alternative = "greater")
stopifnot(isequal(round(statistic(itao), 4), 1.8335)) # test statistic
stopifnot(isequal(round(pvalue(itao), 4), 0.0334)) # p-value

## Two-sided asymptotic, p. 636
ita <- independence_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                         data = dta,
                         scores = list(Severity.of.GVHD.Toxicity = scrs))
stopifnot(isequal(round(pvalue(ita), 4), 0.0667)) # p-value

## ## One-sided exact, p. 636
## iteo <- independence_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
##                           data = dta,
##                           scores = list(Severity.of.GVHD.Toxicity = scrs),
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 4), 0.0348)) # p-value

## ## Two-sided exact, p. 636
## ite <- independence_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
##                          data = dta,
##                          scores = list(Severity.of.GVHD.Toxicity = scrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.0659)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 4), 0.0025)) # point prob

## One-sided approximative, p. 636
itxo <- independence_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                          data = dta,
                          scores = list(Severity.of.GVHD.Toxicity = scrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0348 & pci[2] > 0.0348) # p-value

## Two-sided approximative, p. 636
itx <- independence_test(Severity.of.GVHD.Toxicity ~ MHC.Status | Stratum,
                         data = dta,
                         scores = list(Severity.of.GVHD.Toxicity = scrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0659 & pci[2] > 0.0659) # p-value

## Clean-up
rm(gvhd, dta, scrs, itao, ita, #iteo, ite,
   itxo, itx, pci)
################################################################################


################################################################################
### 20      Unordered R x C Contingency Tables
################################################################################

################################################################################
### 20.5    Pearson's Chi-Square Test, p. 663
load("ORAL.rda")

##
## 'chisq_test.table':
##

## Two-sided asymptotic, p. 666
cta <- chisq_test(oral)
stopifnot(isequal(round(statistic(cta), 4), 22.0992)) # test statistic
stopifnot(isequal(round(pvalue(cta), 4), 0.1400)) # p-value

## ## Two-sided exact, p. 666
## cte <- chisq_test(oral,
##                   distribution = "exact")
## stopifnot(isequal(round(pvalue(cte), 4), 0.0269)) # p-value
## stopifnot(isequal(round(dperm(cte, statistic(cte)), 4), 0.0000)) # point prob

## Two-sided approximative, p. 666 (667)
ctx <- chisq_test(oral,
                  distribution = approximate(B = 10000))
pci <- attr(pvalue(ctx), "conf.int")
stopifnot(pci[1] < 0.0269 & pci[2] > 0.0269) # p-value

## Clean-up
rm(oral, cta, #cte,
   ctx, pci)
################################################################################


################################################################################
### 21      Single Ordered R x C Contingency Tables
################################################################################

################################################################################
### 21.5    The Kruskal-Wallis Test, p. 678
load("TUMOR.rda")

dta <- table2df(tumor)
dta$Regression <- as.numeric(dta$Regression)

##
## 'kruskal_test.formula':
##

## Two-sided asymptotic, p. 680
kta <- kruskal_test(Regression ~ CHEMO, data = dta)
stopifnot(isequal(round(statistic(kta), 4), 8.6824)) # test statistic
stopifnot(isequal(round(pvalue(kta), 4), 0.0695)) # p-value

## ## Two-sided exact, p. 680
## kte <- kruskal_test(Regression ~ CHEMO, data = dta,
##                     distribution = "exact")
## stopifnot(isequal(round(pvalue(kte), 4), 0.0390)) # p-value
## stopifnot(isequal(round(dperm(kte, statistic(kte)), 4), 0.0015)) # point prob

## Two-sided approximative, p. 680
ktx <- kruskal_test(Regression ~ CHEMO, data = dta,,
                    distribution = approximate(B = 10000))
pci <- attr(pvalue(ktx), "conf.int")
stopifnot(pci[1] < 0.0390 & pci[2] > 0.0390) # p-value

## Clean-up
rm(tumor, dta, kta, #kte,
   ktx, pci)
################################################################################


################################################################################
### 21.6    The Normal Scores Test, p. 682
load("TUMOR.rda")

dta <- table2df(tumor)
dta$Regression <- as.numeric(dta$Regression)

##
## 'independence_test.formula':
##

## Two-sided asymptotic, p. 684
ita <- independence_test(Regression ~ CHEMO, data = dta,
                         ytrafo = function(data)
                             trafo(data, numeric_trafo = function(y)
                                 normal_trafo(y, ties.method = "average")),
                         teststat = "quad")
stopifnot(isequal(round(statistic(ita), 4), 8.8528)) # test statistic
stopifnot(isequal(round(pvalue(ita), 4), 0.0649)) # p-value

## ## Two-sided exact, p. 684
## ite <- independence_test(Regression ~ CHEMO, data = dta,
##                          ytrafo = function(data)
##                              trafo(data, numeric_trafo = function(y)
##                                  normal_trafo(y, ties.method = "average")),
##                          teststat = "quad",
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.03890)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 4), 0.0015)) # point prob

## Two-sided approximative, p. 684
itx <- independence_test(Regression ~ CHEMO, data = dta,
                         ytrafo = function(data)
                             trafo(data, numeric_trafo = function(y)
                                 normal_trafo(y, ties.method = "average")),
                         teststat = "quad",
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0390 & pci[2] > 0.0390) # p-value

## Clean-up
rm(tumor, dta, ita, #ite,
   itx, pci)
################################################################################


################################################################################
### 21.7    The Savage Scores Test, p. 686
load("TUMOR.rda")

dta <- table2df(tumor)
dta$Regression <- as.numeric(dta$Regression)

##
## 'independence_test.formula':
##

## Two-sided asymptotic, p. 687
sta <- surv_test(Surv(Regression) ~ CHEMO, data = dta, ties.method = "average")
stopifnot(isequal(round(statistic(sta), 4), 9.3366)) # test statistic
stopifnot(isequal(round(pvalue(sta), 4), 0.0532)) # p-value

## ## Two-sided exact, p. 687
## ste <- surv_test(Surv(Regression) ~ CHEMO, data = dta, ties.method = "average",
##                  distribution = "exact")
## stopifnot(isequal(round(pvalue(ste), 4), 0.0281)) # p-value
## stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0015)) # point prob

## Two-sided approximative, p. 687 (688)
stx <- surv_test(Surv(Regression) ~ CHEMO, data = dta, ties.method = "average",
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0281 & pci[2] > 0.0281) # p-value

## Clean-up
rm(tumor, dta, sta, #ste,
   stx, pci)
################################################################################


################################################################################
### 21.8    One-Way ANOVA with Arbitrary Scores, p.689
load("TUMOR.rda")

dta <- table2df(tumor)
scrs <- c(0, 100, 150)

##
## 'independence_test.formula':
##

## Two-sided asymptotic, p. 691
ita <- independence_test(Regression ~ CHEMO, data = dta,
                         scores = list(Regression = scrs), teststat = "quad")
stopifnot(isequal(round(statistic(ita), 4), 8.5069)) # test statistic
stopifnot(isequal(round(pvalue(ita), 4), 0.0747)) # p-value

## ## Two-sided exact, p. 691
## ite <- independence_test(Regression ~ CHEMO, data = dta,
##                          scores = list(Regression = scrs), teststat = "quad",
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.0450)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 4), 0.0021)) # point prob

## Two-sided approximative, p. 691 (693)
itx <- independence_test(Regression ~ CHEMO, data = dta,
                         scores = list(Regression = scrs), teststat = "quad",
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0450 & pci[2] > 0.0450) # p-value

## Clean-up
rm(tumor, dta, scrs, ita, #ite,
   itx, pci)
################################################################################


################################################################################
### 22      Doubly Ordered R x C Contingency Tables
################################################################################

################################################################################
### 22.6    Linear-by-Linear Association Test, p.703

### (a) Equally Spaced Column Scores: (0, 1, 2, 3)
load("TOX.rda")

dta <- table2df(tox)
tscrs <- c(100, 200, 300, 400)
dscrs <- 1:4

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 705
itao <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                          scores = list("Drug.Toxicity" = tscrs,
                                        "Drug.Dose" = dscrs),
                          alternative = "greater")
stopifnot(isequal(round(statistic(itao), 4), 1.8066)) # test statistic
stopifnot(isequal(round(pvalue(itao), 4), 0.0354)) # p-value

## Two-sided asymptotic, p. 705
ita <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                         scores = list("Drug.Toxicity" = tscrs,
                                       "Drug.Dose" = dscrs))
stopifnot(isequal(round(pvalue(ita), 4), 0.0708)) # p-value

## ## One-sided exact, p. 705
## iteo <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
##                           scores = list("Drug.Toxicity" = tscrs,
##                                         "Drug.Dose" = dscrs),
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 4), 0.0442)) # p-value

## ## Two-sided exact, p. 705
## ite <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
##                          scores = list("Drug.Toxicity" = tscrs,
##                                        "Drug.Dose" = dscrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.0792)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 4), 0.0121)) # point prob

## One-sided approximative, p. 705
itxo <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                          scores = list("Drug.Toxicity" = tscrs,
                                        "Drug.Dose" = dscrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0442 & pci[2] > 0.0442) # p-value

## Two-sided approximative, p. 705
itx <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                         scores = list("Drug.Toxicity" = tscrs,
                                       "Drug.Dose" = dscrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0792 & pci[2] > 0.0792) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 705
itao2 <- independence_test(tox,
                           scores = list("Drug.Toxicity" = tscrs,
                                         "Drug.Dose" = dscrs),
                           alternative = "greater")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 705
ita2 <- independence_test(tox,
                           scores = list("Drug.Toxicity" = tscrs,
                                         "Drug.Dose" = dscrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## One-sided exact, p. 705
## iteo2 <- independence_test(tox,
##                            scores = list("Drug.Toxicity" = tscrs,
##                                          "Drug.Dose" = dscrs),
##                            distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Two-sided exact, p. 705
## ite2 <- independence_test(tox,
##                           scores = list("Drug.Toxicity" = tscrs,
##                                         "Drug.Dose" = dscrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value
## stopifnot(isequal(dperm(ite2, statistic(ite2)), dperm(ite, statistic(ite)))) # point prob

## One-sided approximative, p. 705
itxo2 <- independence_test(tox,
                           scores = list("Drug.Toxicity" = tscrs,
                                         "Drug.Dose" = dscrs),
                           distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.0442 & pci[2] > 0.0442) # p-value

## Two-sided approximative, p. 705
itx2 <- independence_test(tox,
                          scores = list("Drug.Toxicity" = tscrs,
                                        "Drug.Dose" = dscrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.0792 & pci[2] > 0.0792) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 705
lta <- lbl_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                scores = list("Drug.Toxicity" = tscrs,
                              "Drug.Dose" = dscrs))
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 705
ltx <- lbl_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                scores = list("Drug.Toxicity" = tscrs,
                              "Drug.Dose" = dscrs),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.0792 & pci[2] > 0.0792) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 705
lta2 <- lbl_test(tox, scores = list("Drug.Toxicity" = tscrs,
                                    "Drug.Dose" = dscrs))
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 705
ltx2 <- lbl_test(tox, scores = list("Drug.Toxicity" = tscrs,
                                    "Drug.Dose" = dscrs),
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.0792 & pci[2] > 0.0792) # p-value

## Clean-up
rm(tox, dta, tscrs, dscrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)

### (b) Unequally Spaced Column Scores: (1, 3, 9, 27)
load("TOX.rda")

dta <- table2df(tox)
tscrs <- c(1, 3, 9, 27)
dscrs <- 1:4

##
## 'independence_test.formula':
##

## One-sided asymptotic, p. 705
itao <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                          scores = list("Drug.Toxicity" = tscrs,
                                        "Drug.Dose" = dscrs),
                          alternative = "greater")
stopifnot(isequal(round(statistic(itao), 4), 1.7344)) # test statistic
stopifnot(isequal(round(pvalue(itao), 4), 0.0414)) # p-value

## Two-sided asymptotic, p. 705
ita <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                         scores = list("Drug.Toxicity" = tscrs,
                                       "Drug.Dose" = dscrs))
stopifnot(isequal(round(pvalue(ita), 4), 0.0828)) # p-value

## ## One-sided exact, p. 705
## iteo <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
##                           scores = list("Drug.Toxicity" = tscrs,
##                                         "Drug.Dose" = dscrs),
##                           distribution = "exact", alternative = "greater")
## stopifnot(isequal(round(pvalue(iteo), 4), 0.0501)) # p-value

## ## Two-sided exact, p. 705
## ite <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
##                          scores = list("Drug.Toxicity" = tscrs,
##                                        "Drug.Dose" = dscrs),
##                          distribution = "exact")
## stopifnot(isequal(round(pvalue(ite), 4), 0.0780)) # p-value
## stopifnot(isequal(round(dperm(ite, statistic(ite)), 4), 0.0052)) # point prob

## One-sided approximative, p. 705
itxo <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                          scores = list("Drug.Toxicity" = tscrs,
                                        "Drug.Dose" = dscrs),
                          distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo), "conf.int")
stopifnot(pci[1] < 0.0501 & pci[2] > 0.0501) # p-value

## Two-sided approximative, p. 705
itx <- independence_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                         scores = list("Drug.Toxicity" = tscrs,
                                       "Drug.Dose" = dscrs),
                         distribution = approximate(B = 10000))
pci <- attr(pvalue(itx), "conf.int")
stopifnot(pci[1] < 0.0780 & pci[2] > 0.0780) # p-value

##
## 'independence_test.table':
##

## One-sided asymptotic, p. 705
itao2 <- independence_test(tox,
                           scores = list("Drug.Toxicity" = tscrs,
                                         "Drug.Dose" = dscrs),
                           alternative = "greater")
stopifnot(isequal(statistic(itao2), statistic(itao))) # test statistic
stopifnot(isequal(pvalue(itao2), pvalue(itao))) # p-value

## Two-sided asymptotic, p. 705
ita2 <- independence_test(tox,
                           scores = list("Drug.Toxicity" = tscrs,
                                         "Drug.Dose" = dscrs))
stopifnot(isequal(pvalue(ita2), pvalue(ita))) # p-value

## One-sided exact, p. 705
## iteo2 <- independence_test(tox,
##                            scores = list("Drug.Toxicity" = tscrs,
##                                          "Drug.Dose" = dscrs),
##                            distribution = "exact", alternative = "greater")
## stopifnot(isequal(pvalue(iteo2), pvalue(iteo))) # p-value

## ## Two-sided exact, p. 705
## ite2 <- independence_test(tox,
##                           scores = list("Drug.Toxicity" = tscrs,
##                                         "Drug.Dose" = dscrs),
##                           distribution = "exact")
## stopifnot(isequal(pvalue(ite2), pvalue(ite))) # p-value
## stopifnot(isequal(dperm(ite2, statistic(ite2)), dperm(ite, statistic(ite)))) # point prob

## One-sided approximative, p. 705
itxo2 <- independence_test(tox,
                           scores = list("Drug.Toxicity" = tscrs,
                                         "Drug.Dose" = dscrs),
                           distribution = approximate(B = 10000), alternative = "greater")
pci <- attr(pvalue(itxo2), "conf.int")
stopifnot(pci[1] < 0.0501 & pci[2] > 0.0501) # p-value

## Two-sided approximative, p. 705
itx2 <- independence_test(tox,
                          scores = list("Drug.Toxicity" = tscrs,
                                        "Drug.Dose" = dscrs),
                          distribution = approximate(B = 10000))
pci <- attr(pvalue(itx2), "conf.int")
stopifnot(pci[1] < 0.0780 & pci[2] > 0.0780) # p-value

##
## 'lbl_test.formula':
##

## Two-sided asymptotic, p. 705
lta <- lbl_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                scores = list("Drug.Toxicity" = tscrs,
                              "Drug.Dose" = dscrs))
stopifnot(isequal(statistic(lta), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 705
ltx <- lbl_test(Drug.Toxicity ~ Drug.Dose, data = dta,
                scores = list("Drug.Toxicity" = tscrs,
                              "Drug.Dose" = dscrs),
                distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx), "conf.int")
stopifnot(pci[1] < 0.0780 & pci[2] > 0.0780) # p-value

##
## 'lbl_test.table':
##

## Two-sided asymptotic, p. 705
lta2 <- lbl_test(tox, scores = list("Drug.Toxicity" = tscrs,
                                    "Drug.Dose" = dscrs))
stopifnot(isequal(statistic(lta2), unname(statistic(ita)^2))) # test statistic

## Two-sided approximative, p. 705
ltx2 <- lbl_test(tox, scores = list("Drug.Toxicity" = tscrs,
                                    "Drug.Dose" = dscrs),
                 distribution = approximate(B = 10000))
pci <- attr(pvalue(ltx2), "conf.int")
stopifnot(pci[1] < 0.0780 & pci[2] > 0.0780) # p-value

## Clean-up
rm(tox, dta, tscrs, dscrs, itao, ita, #iteo, ite,
   itxo, itx, itao2, ita2, #iteo2, ite2,
   itxo2, itx2, lta, ltx, lta2, ltx2, pci)
################################################################################


################################################################################
### 23      Stratified R x C Contingency Tables
################################################################################

################################################################################
### 23.4    Example: Types of Death in the U.S. Military, p.720
load("ARMY.rda")

##
## 'cmh_test.table':
##

## Two-sided asymptotic, p. 723
cta <- cmh_test(army)
stopifnot(isequal(round(statistic(cta), 4), 25.1846)) # test statistic
stopifnot(isequal(round(pvalue(cta), 4), 0.0028)) # p-value

## Clean-up
rm(army, cta)
################################################################################

################################################################################
### 23.6    Unordered Stratified R x C Table, p. 724
load("JOB.rda")

##
## 'cmh_test.table':
##

## Two-sided asymptotic, p. 726
cta <- cmh_test(job)
stopifnot(isequal(round(statistic(cta), 4), 10.2001)) # test statistic
stopifnot(isequal(round(pvalue(cta), 4), 0.3345)) # p-value

## Clean-up
rm(job, cta)
################################################################################

################################################################################
### 23.7    Singly Ordered R x C Table, p. 726

### (a) Equally Spaced Column Scores: (1, 2, 3, 4)
load("JOB.rda")

scrs <- c(1, 2, 3, 4)

##
## 'cmh_test.table':
##

## Two-sided asymptotic, p. 723
cta <- cmh_test(job,
                scores = list("Job.Satisfaction" = 1:4))
stopifnot(isequal(round(statistic(cta), 4), 9.2259)) # test statistic
stopifnot(isequal(round(pvalue(cta), 4), 0.0264)) # p-value

## Clean-up
rm(army, cta)
################################################################################

# teststatistic, page 1018, is _wrong_ (but StatXact itself
# computes the correct result)
# (isequal(round(statistic(cta), 3), 9.226))
# Agresti, 2002, Table 7.12, page 297
stopifnot(isequal(round(statistic(cta), 4), 9.0342))

# asymptotical p-value, page 1018, is _wrong_ (but StatXact itself
# computes the correct result)
# (isequal(round(pvalue(cta), 4), 0.02643))
# Agresti, 2002, Table 7.12, page 297
stopifnot(isequal(round(pvalue(cta), 4), 0.0288))
################################################################################

################################################################################
### 23.8    Doubly Ordered R x C Table, p. 727

### (a) Equally Spaced Column Scores: (1, 2, 3, 4)
load("JOB.rda")

jscrs <- c(1, 3, 4, 5)
iscrs <- c(3, 10, 20, 35)

##
## 'independence_test.table
##

## One-sided asymptotic, p. 728
itao <- independence_test(job,
                          scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                                        Income..Dollars. = c(3, 10, 20, 35)))

##
## 'lbl_test.table
##

lta <- lbl_test(jobsatisfaction,
                scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                              Income..Dollars. = c(3, 10, 20, 35)))

# teststatistic, page 1020
stopifnot(isequal(round(sqrt(statistic(lta)), 3), 2.481))

# asymptotical p-value, page 1020
stopifnot(isequal(round(pvalue(lta), 5), 0.01309))

lta <- cmh_test(jobsatisfaction,
    scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                  Income = c(3, 10, 20, 35)))

# teststatistic, page 1020
stopifnot(isequal(round(sqrt(statistic(lta)), 3), 2.481))

# asymptotical p-value, page 1020
stopifnot(isequal(round(pvalue(lta), 5), 0.01309))


### --------------------------------------------------------- ###

# some additional checks (always add new tests at the end because of the RNG's)

lta <- surv_test(Surv(time, event) ~ treatment | gender, data = srv)
stopifnot(isequal(round(pvalue(lta), 4), 0.1224))

### example from Callaert (2003), AmStat 57, 214-217
exdata <- data.frame(time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
                     event = rep(TRUE, 15),
                     group = factor(c(rep(0, 7), rep(1, 8))))
p <- pvalue(surv_test(Surv(time, event) ~ group, data = exdata,
          distribution = exact()))
stopifnot(isequal(round(p, 4), 0.0505))
p <- pvalue(surv_test(Surv(time, event) ~ group, data = exdata,
          distribution = exact(), ties = "average"))
stopifnot(isequal(round(p, 4), 0.0468))





################################################################################


### <FIXME>  Rename to azt1, or ???
load("AIDS.rda")
### </FIXME>
diff <- with(AIDS[1:8, ], post - pre)
y <- as.vector(t(cbind(abs(diff) * (diff < 0), abs(diff) * (diff >= 0))))
x <- factor(rep(c("neg", "pos"), length(diff)), levels = c("pos", "neg"))
block <- gl(length(diff), 2)

## One-sided asymptotic, p. 126
stao <- symmetry_test(y ~ x | block,
                      alternative = "less")
stopifnot(isequal(round(statistic(stao), 4), -1.7072)) # test statistic
stopifnot(isequal(round(pvalue(stao), 4), 0.0439)) # p-value

## Two-sided asymptotic, p. 126
sta <- symmetry_test(y ~ x | block)
stopifnot(isequal(round(pvalue(sta), 4), 0.0878)) # p-value

## One-sided exact, p. 126
steo <- symmetry_test(y ~ x | block,
                      distribution = "exact", alternative = "less")
stopifnot(isequal(round(pvalue(steo), 4), 0.0011)) # p-value

## Two-sided exact, p. 126
ste <- symmetry_test(y ~ x | block,
                     distribution = "exact")
stopifnot(isequal(round(pvalue(ste), 4), 0.0021)) # p-value
stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0000)) # point prob

## One-sided approximate, p. 126
stxo <- symmetry_test(y ~ x | block,
                      distribution = approximate(B = 10000), alternative = "less")
pci <- attr(pvalue(stxo), "conf.int")
stopifnot(pci[1] < 0.0011 & pci[2] > 0.0011) # p-value

## Two-sided approximate, p. 126
stx <- symmetry_test(y ~ x | block,
                     distribution = approximate(B = 10000))
pci <- attr(pvalue(stx), "conf.int")
stopifnot(pci[1] < 0.0021 & pci[2] > 0.0021) # p-value
################################################################################







## ## Two-sided exact, p. 248
fte <- independence_test(y ~ x | block,
                         distribution = "exact")
stopifnot(isequal(round(pvalue(fte), 4), 0.6875)) # p-value
stopifnot(isequal(round(dperm(fte, statistic(fte)), 4), 0.4688)) # point prob


load("AIDS.rda")
### </FIXME>
diff <- with(AIDS[1:8, ], post - pre)
diff0 <- diff[abs(diff) > 0]
y <- as.vector(t(cbind(abs(diff0) * (diff0 < 0), abs(diff0) * (diff0 >= 0))))
x <- factor(rep(c("neg", "pos"), length(diff0)), levels = c("pos", "neg"))
block <- gl(length(diff0), 2)

fta <-
    friedman_test(y ~ x | block)

fta






#######
### HÄR!
#######























################################################################################
### Workaround for exact tests with blocks
## One-sided exact, p. 219
steo <- independence_test(Surv(response, censor) ~ trtmnt | gender,
                          data = sur_trt, distr = "exact", alternative = "less",
                          ytrafo = function(data)
                              trafo(data, surv_trafo = function(y)
                                  round(logrank_trafo(y, ties = "average"), 5) * 1e5))
stopifnot(isequal(round(pvalue(steo), 4), 0.0900)) # p-value

## Two-sided exact, p. 219
ste <- independence_test(Surv(response, censor) ~ trtmnt | gender,
                          data = sur_trt, distr = "exact",
                          ytrafo = function(data)
                              trafo(data, surv_trafo = function(y)
                                  round(logrank_trafo(y, ties = "average"), 5) * 1e5))
stopifnot(isequal(round(pvalue(ste), 4), 0.1600)) # p-value
stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0100)) # point prob
## <FIXME>  Why is this???  </FIXME>


## ## One-sided exact, p. 134
yy <- round(as.numeric(as.character(y)), 4) * 1e4

independence_test(yy ~ x | block, distr = "exact", alt = "greater")
stopifnot(isequal(round(pvalue(steo), 4), 0.0001)) # p-value

## Two-sided exact, p. 134
ste <- independencesymmetry_test(y ~ x | block, scores = list(y = doses),
                     distr = exact(fact = 1e5))
stopifnot(isequal(round(pvalue(ste), 4), 0.0001)) # p-value
stopifnot(isequal(round(dperm(ste, statistic(ste)), 4), 0.0000)) # point prob


y <- with(sur_trt, logrank_trafo(Surv(response, censor), ties = "average"))
yy <- round(y, 5) * 1e5

independence_test(yy ~ trtmnt | gender, data = sur_trt, distr = "exact")

yyy <- yy * 10000

independence_test(yyy ~ trtmnt | gender, data = sur_trt, distr = "exact")

independence_test(yy ~ trtmnt | gender, data = sur_trt,
                  distr = exact(algo = "shift", fact = 100000))

distr = exact(algo = "shift", fact = 1e7))

stopifnot(isequal(round(pvalue(sta), 4), 0.1224)) # p-value





mat <- as.table(matrix(c(21, 9, 2, 12), nrow = 2, byrow = T))
mat0 <- mat
#diag(mat0) <- 0 # no contribution from the diagonal elements
dta <- table2df_sym(mat0)
y <- dta$response
x <- factor(dta$groups)#, levels = c("Var1", "B"))
block <- factor(rep.int(seq_len(sum(mat0)), 2))

symmetry_test(y ~ x | block)
symmetry_test(y ~ x | block,  distr = "exact")
