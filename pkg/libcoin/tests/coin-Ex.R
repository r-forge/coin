pkgname <- "coin"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('coin')
library("libcoin")

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CWD")
### * CWD

flush(stderr()); flush(stdout())

### Name: CWD
### Title: Coarse Woody Debris
### Aliases: CWD
### Keywords: datasets

### ** Examples

## Zeileis and Hothorn (2013, pp. 942-944)
## Approximative (Monte Carlo) maximally selected statistics
CWD[1:6] <- 100 * CWD[1:6] # scaling (to avoid harmless warning)
mt <- lc("maxstat_test",sample2 + sample3 + sample4 +
                   sample6 + sample7 + sample8 ~ trend, data = CWD,
                   distribution = approximate(B = 100000))

## Absolute maximum of standardized statistics (t = 3.08)
statistic(mt)

## 5 % critical value (t_0.05 = 2.86)
(c <- qperm(mt, 0.95))

## Only 'sample8' exceeds the 5 % critical value
sts <- statistic(mt, "standardized")
idx <- which(sts > c, arr.ind = TRUE)
sts[unique(idx[, 1]), unique(idx[, 2]), drop = FALSE]



cleanEx()
nameEx("ContingencyTests")
### * ContingencyTests

flush(stderr()); flush(stdout())

### Name: ContingencyTests
### Title: Tests of Independence in Two- or Three-Way Contingency Tables
### Aliases: chisq_test chisq_test.formula chisq_test.table
###   chisq_test.IndependenceProblem cmh_test cmh_test.formula
###   cmh_test.table cmh_test.IndependenceProblem lbl_test lbl_test.formula
###   lbl_test.table lbl_test.IndependenceProblem
### Keywords: htest

### ** Examples

## Example data
## Davis (1986, p. 140)
davis <- matrix(
    c(3,  6,
      2, 19),
    nrow = 2, byrow = TRUE
)
davis <- as.table(davis)

## Asymptotic Pearson chi-squared test
lc("chisq_test",davis)

## Approximative (Monte Carlo) Pearson chi-squared test
ct <- lc("chisq_test",davis,
                 distribution = approximate(B = 10000))
pvalue(ct)          # standard p-value
midpvalue(ct)       # mid-p-value
pvalue_interval(ct) # p-value interval

## Exact Pearson chi-squared test (Davis, 1986)
## Note: disagrees with Fisher's exact test
ct <- lc("chisq_test",davis,
                 distribution = "exact")
pvalue(ct)          # standard p-value
midpvalue(ct)       # mid-p-value
pvalue_interval(ct) # p-value interval
fisher.test(davis)


## Laryngeal cancer data
## Agresti (2002, p. 107, Tab. 3.13)
cancer <- matrix(
    c(21, 2,
      15, 3),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        "Treatment" = c("Surgery", "Radiation"),
           "Cancer" = c("Controlled", "Not Controlled")
    )
)
cancer <- as.table(cancer)

## Exact Pearson chi-squared test (Agresti, 2002, p. 108, Tab. 3.14)
## Note: agrees with Fishers's exact test
(ct <- lc("chisq_test",cancer,
                  distribution = "exact"))
midpvalue(ct)       # mid-p-value
pvalue_interval(ct) # p-value interval
fisher.test(cancer)


## Homework conditions and teacher's rating
## Yates (1948, Tab. 1)
yates <- matrix(
    c(141, 67, 114, 79, 39,
      131, 66, 143, 72, 35,
       36, 14,  38, 28, 16),
    byrow = TRUE, ncol = 5,
    dimnames = list(
           "Rating" = c("A", "B", "C"),
        "Condition" = c("A", "B", "C", "D", "E")
    )
)
yates <- as.table(yates)

## Asymptotic Pearson chi-squared test (Yates, 1948, p. 176)
lc("chisq_test",yates)

## Asymptotic Pearson-Yates chi-squared test (Yates, 1948, pp. 180-181)
## Note: 'Rating' and 'Condition' as ordinal
(ct <- lc("chisq_test",yates,
                  alternative = "less",
                  scores = list("Rating" = c(-1, 0, 1),
                                "Condition" = c(2, 1, 0, -1, -2))))
statistic(ct)^2 # chi^2 = 2.332

## Asymptotic Pearson-Yates chi-squared test (Yates, 1948, p. 181)
## Note: 'Rating' as ordinal
lc("chisq_test",yates,
           scores = list("Rating" = c(-1, 0, 1))) # Q = 3.825


## Change in clinical condition and degree of infiltration
## Cochran (1954, Tab. 6)
cochran <- matrix(
    c(11,  7,
      27, 15,
      42, 16,
      53, 13,
      11,  1),
    byrow = TRUE, ncol = 2,
    dimnames = list(
              "Change" = c("Marked", "Moderate", "Slight",
                           "Stationary", "Worse"),
        "Infiltration" = c("0-7", "8-15")
    )
)
cochran <- as.table(cochran)

## Asymptotic Pearson chi-squared test (Cochran, 1954, p. 435)
lc("chisq_test",cochran) # X^2 = 6.88

## Asymptotic Cochran-Armitage test (Cochran, 1954, p. 436)
## Note: 'Change' as ordinal
(ct <- lc("chisq_test",cochran,
                  scores = list("Change" = c(3, 2, 1, 0, -1))))
statistic(ct)^2 # X^2 = 6.66


## Change in size of ulcer crater for two treatment groups
## Armitage (1955, Tab. 2)
armitage <- matrix(
    c( 6, 4, 10, 12,
      11, 8,  8,  5),
    byrow = TRUE, ncol = 4,
    dimnames = list(
        "Treatment" = c("A", "B"),
           "Crater" = c("Larger", "< 2/3 healed",
                        "=> 2/3 healed", "Healed")
    )
)
armitage <- as.table(armitage)

## Approximative (Monte Carlo) Pearson chi-squared test (Armitage, 1955, p. 379)
lc("chisq_test",armitage,
           distribution = approximate(B = 10000)) # chi^2 = 5.91

## Approximative (Monte Carlo) Cochran-Armitage test (Armitage, 1955, p. 379)
(ct <- lc("chisq_test",armitage,
                  distribution = approximate(B = 10000),
                  scores = list("Crater" = c(-1.5, -0.5, 0.5, 1.5))))
statistic(ct)^2 # chi_0^2 = 5.26


## Relationship between job satisfaction and income stratified by gender
## Agresti (2002, p. 288, Tab. 7.8)

## Asymptotic generalized Cochran-Mantel-Haenszel test (Agresti, p. 297)
lc("cmh_test",jobsatisfaction) # CMH = 10.2001

## Asymptotic generalized Cochran-Mantel-Haenszel test (Agresti, p. 297)
## Note: 'Job.Satisfaction' as ordinal
lc("cmh_test",jobsatisfaction,
         scores = list("Job.Satisfaction" = c(1, 3, 4, 5))) # L^2 = 9.0342

## Asymptotic linear-by-linear association test (Agresti, p. 297)
## Note: 'Job.Satisfaction' and 'Income' as ordinal
(lt <- lc("lbl_test",jobsatisfaction,
                scores = list("Job.Satisfaction" = c(1, 3, 4, 5),
                              "Income" = c(3, 10, 20, 35))))
statistic(lt)^2 # M^2 = 6.1563



cleanEx()
nameEx("CorrelationTests")
### * CorrelationTests

flush(stderr()); flush(stdout())

### Name: CorrelationTests
### Title: Correlation Tests
### Aliases: spearman_test spearman_test.formula
###   spearman_test.IndependenceProblem fisyat_test fisyat_test.formula
###   fisyat_test.IndependenceProblem quadrant_test quadrant_test.formula
###   quadrant_test.IndependenceProblem koziol_test koziol_test.formula
###   koziol_test.IndependenceProblem
### Keywords: htest

### ** Examples

## Asymptotic Spearman test
lc("spearman_test",CONT ~ INTG, data = USJudgeRatings)

## Asymptotic Fisher-Yates test
lc("fisyat_test",CONT ~ INTG, data = USJudgeRatings)

## Asymptotic quadrant test
lc("quadrant_test",CONT ~ INTG, data = USJudgeRatings)

## Asymptotic Koziol-Nemec test
lc("koziol_test",CONT ~ INTG, data = USJudgeRatings)



cleanEx()
nameEx("GTSG")
### * GTSG

flush(stderr()); flush(stdout())

### Name: GTSG
### Title: Gastrointestinal Tumor Study Group
### Aliases: GTSG
### Keywords: datasets

### ** Examples

## Plot Kaplan-Meier estimates
plot(survfit(Surv(time / (365.25 / 12), event) ~ group, data = GTSG),
     lty = 1:2, ylab = "% Survival", xlab = "Survival Time in Months")
legend("topright", lty = 1:2,
       c("Chemotherapy+Radiation", "Chemotherapy"), bty = "n")

## Asymptotic logrank test
lc("logrank_test",Surv(time, event) ~ group, data = GTSG)

## Asymptotic Prentice test
lc("logrank_test",Surv(time, event) ~ group, data = GTSG, type = "Prentice")

## Asymptotic test against Weibull-type alternatives (Moreau et al., 1992)
moreau_weight <- function(time, n.risk, n.event)
    1 + log(-log(cumprod(n.risk / (n.risk + n.event))))

lc("independence_test",Surv(time, event) ~ group, data = GTSG,
                  ytrafo = function(data)
                      trafo(data, surv_trafo = function(y)
                          logrank_trafo(y, weight = moreau_weight)))

## Asymptotic test against crossing-curve alternatives (Shen and Le, 2000)
shen_trafo <- function(x)
    ansari_trafo(logrank_trafo(x, type = "Prentice"))

lc("independence_test",Surv(time, event) ~ group, data = GTSG,
                  ytrafo = function(data)
                      trafo(data, surv_trafo = shen_trafo))



cleanEx()
nameEx("IndependenceTest")
### * IndependenceTest

flush(stderr()); flush(stdout())

### Name: IndependenceTest
### Title: General Independence Test
### Aliases: independence_test independence_test.formula
###   independence_test.table independence_test.IndependenceProblem
### Keywords: htest

### ** Examples

## One-sided exact van der Waerden (normal scores) test...
lc("independence_test",asat ~ group, data = asat,
                  ## exact null distribution
                  distribution = "exact",
                  ## one-sided test
                  alternative = "greater",
                  ## apply normal scores to asat$asat
                  ytrafo = function(data)
                      trafo(data, numeric_trafo = normal_trafo),
                  ## indicator matrix of 1st level of asat$group
                  xtrafo = function(data)
                      trafo(data, factor_trafo = function(x)
                          matrix(x == levels(x)[1], ncol = 1)))

## ...or more conveniently
lc("normal_test",asat ~ group, data = asat,
            ## exact null distribution
            distribution = "exact",
            ## one-sided test
            alternative = "greater")


## Receptor binding assay of benzodiazepines
## Johnson, Mercante and May (1993, Tab. 1)
benzos <- data.frame(
      cerebellum = c( 3.41,  3.50,  2.85,  4.43,
                      4.04,  7.40,  5.63, 12.86,
                      6.03,  6.08,  5.75,  8.09,  7.56),
       brainstem = c( 3.46,  2.73,  2.22,  3.16,
                      2.59,  4.18,  3.10,  4.49,
                      6.78,  7.54,  5.29,  4.57,  5.39),
          cortex = c(10.52,  7.52,  4.57,  5.48,
                      7.16, 12.00,  9.36,  9.35,
                     11.54, 11.05,  9.92, 13.59, 13.21),
    hypothalamus = c(19.51, 10.00,  8.27, 10.26,
                     11.43, 19.13, 14.03, 15.59,
                     24.87, 14.16, 22.68, 19.93, 29.32),
        striatum = c( 6.98,  5.07,  3.57,  5.34,
                      4.57,  8.82,  5.76, 11.72,
                      6.98,  7.54,  7.66,  9.69,  8.09),
     hippocampus = c(20.31, 13.20,  8.58, 11.42,
                     13.79, 23.71, 18.35, 38.52,
                     21.56, 18.66, 19.24, 27.39, 26.55),
       treatment = factor(rep(c("Lorazepam", "Alprazolam", "Saline"),
                          c(4, 4, 5)))
)

## Approximative (Monte Carlo) multivariate Kruskal-Wallis test
## Johnson, Mercante and May (1993, Tab. 2)
lc("independence_test",cerebellum + brainstem + cortex +
                  hypothalamus + striatum + hippocampus ~ treatment,
                  data = benzos,
                  teststat = "quadratic",
                  distribution = approximate(B = 10000),
                  ytrafo = function(data)
                      trafo(data, numeric_trafo = rank_trafo)) # Q = 16.129



cleanEx()
nameEx("LocationTests")
### * LocationTests

flush(stderr()); flush(stdout())

### Name: LocationTests
### Title: Two- and K-Sample Location Tests
### Aliases: oneway_test oneway_test.formula
###   oneway_test.IndependenceProblem wilcox_test wilcox_test.formula
###   wilcox_test.IndependenceProblem kruskal_test kruskal_test.formula
###   kruskal_test.IndependenceProblem normal_test normal_test.formula
###   normal_test.IndependenceProblem median_test median_test.formula
###   median_test.IndependenceProblem savage_test savage_test.formula
###   savage_test.IndependenceProblem
### Keywords: htest

### ** Examples
## Don't show: 
options(useFancyQuotes = FALSE)
## End(Don't show)
## Tritiated Water Diffusion Across Human Chorioamnion
## Hollander and Wolfe (1999, p. 110, Tab. 4.1)
diffusion <- data.frame(
    pd = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46,
           1.15, 0.88, 0.90, 0.74, 1.21),
    age = factor(rep(c("At term", "12-26 Weeks"), c(10, 5)))
)

## Exact Wilcoxon-Mann-Whitney test
## Hollander and Wolfe (1999, p. 111)
## (At term - 12-26 Weeks)
(wt <- lc("wilcox_test",pd ~ age, data = diffusion,
                   distribution = "exact", conf.int = TRUE))

## Extract observed Wilcoxon statistic
## Note: this is the sum of the ranks for age = "12-26 Weeks"
statistic(wt, "linear")

## Expectation, variance, two-sided pvalue and confidence interval
expectation(wt)
covariance(wt)
pvalue(wt)
confint(wt)

## For two samples, the Kruskal-Wallis test is equivalent to the W-M-W test
lc("kruskal_test",pd ~ age, data = diffusion,
             distribution = "exact")

## Asymptotic Fisher-Pitman test
lc("oneway_test",pd ~ age, data = diffusion)

## Approximative (Monte Carlo) Fisher-Pitman test
pvalue(lc("oneway_test",pd ~ age, data = diffusion,
                   distribution = approximate(B = 10000)))

## Exact Fisher-Pitman test
pvalue(ot <- lc("oneway_test",pd ~ age, data = diffusion,
                         distribution = "exact"))

## Plot density and distribution of the standardized test statistic
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:2, nrow = 2))
s <- support(ot)
d <- sapply(s, function(x) dperm(ot, x))
p <- sapply(s, function(x) pperm(ot, x))
plot(s, d, type = "S", xlab = "Test Statistic", ylab = "Density")
plot(s, p, type = "S", xlab = "Test Statistic", ylab = "Cum. Probability")
par(op) # reset


## Example data
ex <- data.frame(
    y = c(3, 4, 8, 9, 1, 2, 5, 6, 7),
    x = factor(rep(c("no", "yes"), c(4, 5)))
)

## Boxplots
boxplot(y ~ x, data = ex)

## Exact Brown-Mood median test with different mid-scores
(mt1 <- lc("median_test",y ~ x, data = ex, distribution = "exact"))
(mt2 <- lc("median_test",y ~ x, data = ex, distribution = "exact",
                    mid.score = "0.5"))
(mt3 <- lc("median_test",y ~ x, data = ex, distribution = "exact",
                    mid.score = "1")) # sign change!

## Plot density and distribution of the standardized test statistics
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:3, nrow = 3))
s1 <- support(mt1); d1 <- dperm(mt1, s1)
plot(s1, d1, type = "h", main = "Mid-score: 0",
     xlab = "Test Statistic", ylab = "Density")
s2 <- support(mt2); d2 <- dperm(mt2, s2)
plot(s2, d2, type = "h", main = "Mid-score: 0.5",
     xlab = "Test Statistic", ylab = "Density")
s3 <- support(mt3); d3 <- dperm(mt3, s3)
plot(s3, d3, type = "h", main = "Mid-score: 1",
     xlab = "Test Statistic", ylab = "Density")
par(op) # reset


## Length of YOY Gizzard Shad
## Hollander and Wolfe (1999, p. 200, Tab. 6.3)
yoy <- data.frame(
    length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
               42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
               38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
               31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
    site = gl(4, 10, labels = as.roman(1:4))
)

## Approximative (Monte Carlo) Kruskal-Wallis test
lc("kruskal_test",length ~ site, data = yoy,
             distribution = approximate(B = 10000))

## Approximative (Monte Carlo) Nemenyi-Damico-Wolfe-Dunn test (joint ranking)
## Hollander and Wolfe (1999, p. 244)
## (where Steel-Dwass results are given)
it <- lc("independence_test",length ~ site, data = yoy,
                        distribution = approximate(B = 50000),
                        ytrafo = function(data)
                            trafo(data, numeric_trafo = rank_trafo),
                        xtrafo = mcp_trafo(site = "Tukey"))

## Global p-value
pvalue(it)

## Sites (I = II) != (III = IV) at alpha = 0.01 (p. 244)
pvalue(it, method = "single-step") # subset pivotality is violated


## Asymptotic Jonckheere-Terpstra test for ordered groups
pieces <- data.frame(
    control = c(40, 35, 38, 43, 44, 41),
    rough = c(38, 40, 47, 44, 40, 42),
    accurate = c(48, 40, 45, 43, 46, 44)
)
pieces <- stack(pieces)
pieces$ind <- ordered(pieces$ind,
                      levels = c("control", "rough", "accurate"))

## Look at K: the second line just sums up.
ff <- function(x) {
    K <- multcomp::contrMat(table(x), "Tukey")[, x]
    as.vector(rep(1, nrow(K)) %*% K)
}

lc("independence_test",values ~ ind, data = pieces,
                  alternative = "greater",
                  ytrafo = function(data)
                      trafo(data, numeric_trafo = rank_trafo),
                  xtrafo = function(data)
                      trafo(data, ordered_trafo = ff))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("MarginalHomogeneityTests")
### * MarginalHomogeneityTests

flush(stderr()); flush(stdout())

### Name: MarginalHomogeneityTests
### Title: Marginal Homogeneity Tests
### Aliases: mh_test mh_test.formula mh_test.table mh_test.SymmetryProblem
### Keywords: htest

### ** Examples

## Performance of prime minister
## Agresti (2002, p. 409)
performance <- matrix(
    c(794, 150,
       86, 570),
    nrow = 2, byrow = TRUE,
    dimnames = list(
         "First" = c("Approve", "Disprove"),
        "Second" = c("Approve", "Disprove")
    )
)
performance <- as.table(performance)
diag(performance) <- 0 # speed-up: only off-diagonal elements contribute

## Asymptotic McNemar Test
lc("mh_test",performance)

## Exact McNemar Test
lc("mh_test",performance, distribution = "exact")


## Effectiveness of different media for the growth of diphtheria
## Cochran (1950, Tab. 2)
cases <- c(4, 2, 3, 1, 59)
n <- sum(cases)
cochran <- data.frame(
    diphtheria = factor(
        unlist(rep(list(c(1, 1, 1, 1),
                        c(1, 1, 0, 1),
                        c(0, 1, 1, 1),
                        c(0, 1, 0, 1),
                        c(0, 0, 0, 0)),
                   cases))
    ),
    media = factor(rep(LETTERS[1:4], n)),
    case =  factor(rep(seq_len(n), each = 4))
)

## Asymptotic Cochran Q test (Cochran, 1950, p. 260)
lc("mh_test",diphtheria ~ media | case, data = cochran) # Q = 8.05

## Approximative Cochran Q test
mt <- lc("mh_test",diphtheria ~ media | case, data = cochran,
              distribution = approximate(B = 10000))
pvalue(mt)          # standard p-value
midpvalue(mt)       # mid-p-value
pvalue_interval(mt) # p-value interval


## Opinions on Pre- and Extramarital Sex
## Agresti (2002, p. 421)
opinions <- c("Always wrong", "Almost always wrong",
              "Wrong only sometimes", "Not wrong at all")
PreExSex <- matrix(
    c(144, 33, 84, 126,
        2,  4, 14,  29,
        0,  2,  6,  25,
        0,  0,  1,   5),
    nrow = 4,
    dimnames = list(
          "PreMaritalSex" = opinions,
        "ExtraMaritalSex" = opinions
    )
)
PreExSex <- as.table(PreExSex)

## Asymptotic Stuart test
lc("mh_test",PreExSex)

## Asymptotic Stuart-Birch test
## Note: response as ordinal
lc("mh_test",PreExSex, scores = list(response = 1:length(opinions)))


## Vote intention
## Madansky (1963, pp. 107-108)
vote <- array(
    c(120, 1,  8, 2,   2,  1, 2, 1,  7,
        6, 2,  1, 1, 103,  5, 1, 4,  8,
       20, 3, 31, 1,   6, 30, 2, 1, 81),
    dim = c(3, 3, 3),
    dimnames = list(
          "July" = c("Republican", "Democratic", "Uncertain"),
        "August" = c("Republican", "Democratic", "Uncertain"),
          "June" = c("Republican", "Democratic", "Uncertain")
    )
)
vote <- as.table(vote)

## Asymptotic Madansky test (Q = 70.77)
lc("mh_test",vote)


## Cross-over study
## http://www.nesug.org/proceedings/nesug00/st/st9005.pdf
dysmenorrhea <- array(
    c(6, 2, 1,  4, 3, 0,  5, 2, 2,
      3, 1, 0, 13, 3, 0, 10, 1, 0,
      1, 2, 1,  8, 1, 1, 14, 2, 0),
    dim = c(3, 3, 3),
    dimnames =  list(
          "Placebo" = c("None", "Moderate", "Complete"),
        "High dose" = c("None", "Moderate", "Complete"),
         "Low dose" = c("None", "Moderate", "Complete")
    )
)
dysmenorrhea <- as.table(dysmenorrhea)

## Asymptotic Madansky-Birch test (Q = 53.76)
## Note: response as ordinal
lc("mh_test",dysmenorrhea, scores = list(response = 1:3))

## Asymptotic Madansky-Birch test (Q = 47.29)
## Note: response and measurement conditions as ordinal
lc("mh_test",dysmenorrhea, scores = list(response = 1:3,
                                    conditions = 3:1))



cleanEx()
nameEx("MaximallySelectedStatisticsTests")
### * MaximallySelectedStatisticsTests

flush(stderr()); flush(stdout())

### Name: MaximallySelectedStatisticsTests
### Title: Generalized Maximally Selected Statistics
### Aliases: maxstat_test maxstat_test.formula maxstat_test.table
###   maxstat_test.IndependenceProblem
### Keywords: htest

### ** Examples

## Don't show: 
options(useFancyQuotes = FALSE)
## End(Don't show)
## Tree pipit data (Mueller and Hothorn, 2004)
## Asymptotic maximally selected statistics
lc("maxstat_test",counts ~ coverstorey, data = treepipit)

## Asymptotic maximally selected statistics
## Note: all covariates simultaneously
mt <- lc("maxstat_test",counts ~ ., data = treepipit)
mt@estimates$estimate


## Malignant arrythmias data (Hothorn and Lausen, 2003, Sec. 7.2)
## Asymptotic maximally selected statistics
lc("maxstat_test",Surv(time, event) ~  EF, data = hohnloser,
             ytrafo = function(data)
                 trafo(data, surv_trafo = function(y)
                     logrank_trafo(y, ties.method = "Hothorn-Lausen")))


## Breast cancer data (Hothorn and Lausen, 2003, Sec. 7.3)
## Asymptotic maximally selected statistics
data("sphase", package = "TH.data")
lc("maxstat_test",Surv(RFS, event) ~  SPF, data = sphase,
             ytrafo = function(data)
                 trafo(data, surv_trafo = function(y)
                     logrank_trafo(y, ties.method = "Hothorn-Lausen")))


## Job satisfaction data (Agresti, 2002, p. 288, Tab. 7.8)
## Asymptotic maximally selected statistics
lc("maxstat_test",jobsatisfaction)

## Asymptotic maximally selected statistics
## Note: 'Job.Satisfaction' and 'Income' as ordinal
lc("maxstat_test",jobsatisfaction,
             scores = list("Job.Satisfaction" = 1:4,
                           "Income" = 1:4))



cleanEx()
nameEx("NullDistribution")
### * NullDistribution

flush(stderr()); flush(stdout())

### Name: NullDistribution
### Title: Specification of the Reference Distribution
### Aliases: asymptotic approximate exact
### Keywords: htest

### ** Examples

## Approximative (Monte Carlo) Cochran-Mantel-Haenszel test

## Serial operation
set.seed(123)
lc("cmh_test",disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(B = 100000))

## Not run: 
##D ## Multicore with 8 processes (not for MS Windows)
##D set.seed(123, kind = "L'Ecuyer-CMRG")
##D lc("cmh_test",disease ~ smoking | gender, data = alzheimer,
##D          distribution = approximate(B = 100000,
##D                                     parallel = "multicore", ncpus = 8))
##D 
##D ## Automatic PSOCK cluster with 4 processes
##D set.seed(123, kind = "L'Ecuyer-CMRG")
##D lc("cmh_test",disease ~ smoking | gender, data = alzheimer,
##D          distribution = approximate(B = 100000,
##D                                     parallel = "snow", ncpus = 4))
##D 
##D ## Registered FORK cluster with 12 processes (not for MS Windows)
##D fork12 <- parallel::makeCluster(12, "FORK") # set-up cluster
##D parallel::setDefaultCluster(fork12) # register default cluster
##D set.seed(123, kind = "L'Ecuyer-CMRG")
##D lc("cmh_test",disease ~ smoking | gender, data = alzheimer,
##D          distribution = approximate(B = 100000,
##D                                     parallel = "snow"))
##D parallel::stopCluster(fork12) # clean-up
##D 
##D ## User-specified PSOCK cluster with 8 processes
##D psock8 <- parallel::makeCluster(8, "PSOCK") # set-up cluster
##D set.seed(123, kind = "L'Ecuyer-CMRG")
##D lc("cmh_test",disease ~ smoking | gender, data = alzheimer,
##D          distribution = approximate(B = 100000,
##D                                     parallel = "snow", cl = psock8))
##D parallel::stopCluster(psock8) # clean-up
## End(Not run)



cleanEx()
nameEx("PermutationDistribution-methods")
### * PermutationDistribution-methods

flush(stderr()); flush(stdout())

### Name: PermutationDistribution-methods
### Title: Computation of the Permutation Distribution
### Aliases: dperm dperm-methods dperm,AsymptNullDistribution-method
###   dperm,IndependenceTest-method dperm,NullDistribution-method pperm
###   pperm-methods pperm,AsymptNullDistribution-method
###   pperm,IndependenceTest-method pperm,NullDistribution-method qperm
###   qperm-methods qperm,AsymptNullDistribution-method
###   qperm,IndependenceTest-method qperm,NullDistribution-method rperm
###   rperm-methods rperm,IndependenceTest-method
###   rperm,NullDistribution-method support support-methods
###   support,IndependenceTest-method support,NullDistribution-method
### Keywords: methods htest distribution

### ** Examples

## Two-sample problem
dta <- data.frame(
    y = rnorm(20),
    x = gl(2, 10)
)

## Exact Ansari-Bradley test
at <- lc("ansari_test",y ~ x, data = dta, distribution = "exact")

## Support of the exact distribution of the Ansari-Bradley statistic
supp <- support(at)

## Density of the exact distribution of the Ansari-Bradley statistic
dens <- dperm(at, supp)

## Plotting the density
plot(supp, dens, type = "s")

## 95 % quantile
qperm(at, 0.95)

## One-sided p-value
pperm(at, statistic(at))

## Random number generation
rperm(at, 5)



cleanEx()
nameEx("ScaleTests")
### * ScaleTests

flush(stderr()); flush(stdout())

### Name: ScaleTests
### Title: Two- and K-Sample Scale Tests
### Aliases: taha_test taha_test.formula taha_test.IndependenceProblem
###   klotz_test klotz_test.formula klotz_test.IndependenceProblem
###   mood_test mood_test.formula mood_test.IndependenceProblem ansari_test
###   ansari_test.formula ansari_test.IndependenceProblem fligner_test
###   fligner_test.formula fligner_test.IndependenceProblem conover_test
###   conover_test.formula conover_test.IndependenceProblem
### Keywords: htest

### ** Examples

## Serum Iron Determination Using Hyland Control Sera
## Hollander and Wolfe (1999, p. 147, Tab 5.1)
sid <- data.frame(
    serum = c(111, 107, 100, 99, 102, 106, 109, 108, 104, 99,
              101, 96, 97, 102, 107, 113, 116, 113, 110, 98,
              107, 108, 106, 98, 105, 103, 110, 105, 104,
              100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99),
    method = gl(2, 20, labels = c("Ramsay", "Jung-Parekh"))
)

## Asymptotic Ansari-Bradley test
lc("ansari_test",serum ~ method, data = sid)

## Exact Ansari-Bradley test
pvalue(lc("ansari_test",serum ~ method, data = sid,
                   distribution = "exact"))


## Platelet Counts of Newborn Infants
## Hollander and Wolfe (1999, p. 171, Tab. 5.4)
platelet <- data.frame(
    counts = c(120, 124, 215, 90, 67, 95, 190, 180, 135, 399,
               12, 20, 112, 32, 60, 40),
    treatment = factor(rep(c("Prednisone", "Control"), c(10, 6)))
)

## Approximative (Monte Carlo) Lepage test
## Hollander and Wolfe (1999, p. 172)
lepage_trafo <- function(y)
    cbind("Location" = rank_trafo(y), "Scale" = ansari_trafo(y))

lc("independence_test",counts ~ treatment, data = platelet,
                  distribution = approximate(B = 10000),
                  ytrafo = function(data)
                      trafo(data, numeric_trafo = lepage_trafo),
                  teststat = "quadratic")

## Why was the null hypothesis rejected?
## Note: maximum statistic instead of quadratic form
ltm <- lc("independence_test",counts ~ treatment, data = platelet,
                         distribution = approximate(B = 10000),
                         ytrafo = function(data)
                             trafo(data, numeric_trafo = lepage_trafo))

## Step-down adjustment suggests a difference in location
pvalue(ltm, method = "step-down")

## The same results are obtained from the simple Sidak-Holm procedure since the
## correlation between Wilcoxon and Ansari-Bradley test statistics is zero
cov2cor(covariance(ltm))
pvalue(ltm, method = "step-down", distribution = "marginal", type = "Sidak")



cleanEx()
nameEx("SurvivalTests")
### * SurvivalTests

flush(stderr()); flush(stdout())

### Name: SurvivalTests
### Title: Two- and K-Sample Tests for Censored Data
### Aliases: surv_test logrank_test logrank_test.formula
###   logrank_test.IndependenceProblem
### Keywords: htest survival

### ** Examples

## Example data (Callaert, 2003, Tab.1)
callaert <- data.frame(
    time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
    group = factor(rep(0:1, c(7, 8)))
)

## Logrank scores using mid-ranks (Callaert, 2003, Tab.2)
with(callaert,
     logrank_trafo(Surv(time)))

## Classical asymptotic logrank test (p = 0.0523)
survdiff(Surv(time) ~ group, data = callaert)

## Exact logrank test using mid-ranks (p = 0.0505)
lc("logrank_test",Surv(time) ~ group, data = callaert, distribution = "exact")

## Exact logrank test using average-scores (p = 0.0468)
lc("logrank_test",Surv(time) ~ group, data = callaert, distribution = "exact",
             ties.method = "average-scores")


## Lung cancer data (StatXact 9 manual, p. 213, Tab. 7.19)
lungcancer <- data.frame(
    time = c(257, 476, 355, 1779, 355,
             191, 563, 242, 285, 16, 16, 16, 257, 16),
    event = c(0, 0, 1, 1, 0,
              1, 1, 1, 1, 1, 1, 1, 1, 1),
    group = factor(rep(1:2, c(5, 9)),
                   labels = c("newdrug", "control"))
)

## Logrank scores using average-scores (StatXact 9 manual, p. 214)
with(lungcancer,
     logrank_trafo(Surv(time, event), ties.method = "average-scores"))

## Exact logrank test using average-scores (StatXact 9 manual, p. 215)
lc("logrank_test",Surv(time, event) ~ group, data = lungcancer,
             distribution = "exact", ties.method = "average-scores")

## Exact Prentice test using average-scores (StatXact 9 manual, p. 222)
lc("logrank_test",Surv(time, event) ~ group, data = lungcancer,
             distribution = "exact", ties.method = "average-scores",
             type = "Prentice")


## Approximative (Monte Carlo) versatile test (Lee, 1996)
rho.gamma <- expand.grid(rho = seq(0, 2, 1), gamma = seq(0, 2, 1))
lee_trafo <- function(y)
    logrank_trafo(y, ties.method = "average-scores",
                  type = "Fleming-Harrington",
                  rho = rho.gamma["rho"], gamma = rho.gamma["gamma"])

it <- lc("independence_test",Surv(time, event) ~ group, data = lungcancer,
                        distribution = approximate(B = 10000),
                        ytrafo = function(data)
                            trafo(data, surv_trafo = lee_trafo))
pvalue(it, method = "step-down")



cleanEx()
nameEx("SymmetryTest")
### * SymmetryTest

flush(stderr()); flush(stdout())

### Name: SymmetryTest
### Title: General Symmetry Test
### Aliases: symmetry_test symmetry_test.formula symmetry_test.table
###   symmetry_test.SymmetryProblem
### Keywords: htest

### ** Examples

## One-sided exact Fisher-Pitman test for paired observations
y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
dta <- data.frame(
    y = c(y1, y2),
    x = gl(2, length(y1)),
    block = factor(rep(seq_along(y1), 2))
)

lc("symmetry_test",y ~ x | block, data = dta,
              distribution = "exact", alternative = "greater")

## Alternatively: transform data and set 'paired = TRUE'
delta <- y1 - y2
y <- as.vector(rbind(abs(delta) * (delta >= 0), abs(delta) * (delta < 0)))
x <- factor(rep(0:1, length(delta)), labels = c("pos", "neg"))
block <- gl(length(delta), 2)

lc("symmetry_test",y ~ x | block,
              distribution = "exact", alternative = "greater",
              paired = TRUE)


### Example data
### Gerig (1969, p. 1597)
gerig <- data.frame(
    y1 = c( 0.547, 1.811, 2.561,
            1.706, 2.509, 1.414,
           -0.288, 2.524, 3.310,
            1.417, 0.703, 0.961,
            0.878, 0.094, 1.682,
           -0.680, 2.077, 3.181,
            0.056, 0.542, 2.983,
            0.711, 0.269, 1.662,
           -1.335, 1.545, 2.920,
            1.635, 0.200, 2.065),
    y2 = c(-0.575, 1.840, 2.399,
            1.252, 1.574, 3.059,
           -0.310, 1.553, 0.560,
            0.932, 1.390, 3.083,
            0.819, 0.045, 3.348,
            0.497, 1.747, 1.355,
           -0.285, 0.760, 2.332,
            0.089, 1.076, 0.960,
           -0.349, 1.471, 4.121,
            0.845, 1.480, 3.391),
    x = factor(rep(1:3, 10)),
    b = factor(rep(1:10, each = 3))
)

### Asymptotic multivariate Friedman test
### Gerig (1969, p. 1599)
lc("symmetry_test",y1 + y2 ~ x | b, data = gerig, teststat = "quadratic",
              ytrafo = function(data)
                  trafo(data, numeric_trafo = rank_trafo,
                        block = gerig$b)) # L_n = 17.238

### Asymptotic multivariate Page test
(st <- lc("symmetry_test",y1 + y2 ~ x | b, data = gerig,
                     ytrafo = function(data)
                         trafo(data, numeric_trafo = rank_trafo,
                               block = gerig$b),
                     scores = list(x = 1:3)))
pvalue(st, method = "step-down")



cleanEx()
nameEx("SymmetryTests")
### * SymmetryTests

flush(stderr()); flush(stdout())

### Name: SymmetryTests
### Title: Symmetry Tests
### Aliases: sign_test sign_test.formula sign_test.SymmetryProblem
###   wilcoxsign_test wilcoxsign_test.formula
###   wilcoxsign_test.SymmetryProblem friedman_test friedman_test.formula
###   friedman_test.SymmetryProblem quade_test quade_test.formula
###   quade_test.SymmetryProblem
### Keywords: htest

### ** Examples

## Example data from ?wilcox.test
y1 <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y2 <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)

## One-sided exact sign test
(st <- lc("sign_test",y1 ~ y2, distribution = "exact",
                 alternative = "greater"))
midpvalue(st) # mid-p-value

## One-sided exact Wilcoxon signed-rank test
(wt <- lc("wilcoxsign_test",y1 ~ y2, distribution = "exact",
                       alternative = "greater"))
statistic(wt, type = "linear")
midpvalue(wt) # mid-p-value

## Comparison with R's wilcox.test() function
wilcox.test(y1, y2, paired = TRUE, alternative = "greater")


## Data with explicit group and block information
dta <- data.frame(y = c(y1, y2), x = gl(2, length(y1)),
                  block = factor(rep(seq_along(y1), 2)))

## For two samples, the sign test is equivalent to the Friedman test...
lc("sign_test",y ~ x | block, data = dta, distribution = "exact")
lc("friedman_test",y ~ x | block, data = dta, distribution = "exact")

## ...and the signed-rank test is equivalent to the Quade test
lc("wilcoxsign_test",y ~ x | block, data = dta, distribution = "exact")
lc("quade_test",y ~ x | block, data = dta, distribution = "exact")


## Comparison of three methods ("round out", "narrow angle", and "wide angle")
## for rounding first base.
## Hollander and Wolfe (1999, p. 274, Tab. 7.1)
rounding <- data.frame(
    times = c(5.40, 5.50, 5.55,
              5.85, 5.70, 5.75,
              5.20, 5.60, 5.50,
              5.55, 5.50, 5.40,
              5.90, 5.85, 5.70,
              5.45, 5.55, 5.60,
              5.40, 5.40, 5.35,
              5.45, 5.50, 5.35,
              5.25, 5.15, 5.00,
              5.85, 5.80, 5.70,
              5.25, 5.20, 5.10,
              5.65, 5.55, 5.45,
              5.60, 5.35, 5.45,
              5.05, 5.00, 4.95,
              5.50, 5.50, 5.40,
              5.45, 5.55, 5.50,
              5.55, 5.55, 5.35,
              5.45, 5.50, 5.55,
              5.50, 5.45, 5.25,
              5.65, 5.60, 5.40,
              5.70, 5.65, 5.55,
              6.30, 6.30, 6.25),
    methods = factor(rep(1:3, 22),
                     labels = c("Round Out", "Narrow Angle", "Wide Angle")),
    block = gl(22, 3)
)

## Asymptotic Friedman test
lc("friedman_test",times ~ methods | block, data = rounding)

## Parallel coordinates plot
with(rounding, {
    matplot(t(matrix(times, ncol = 3, byrow = TRUE)),
            type = "l", lty = 1, col = 1, ylab = "Time", xlim = c(0.5, 3.5),
            axes = FALSE)
    axis(1, at = 1:3, labels = levels(methods))
    axis(2)
})

## Where do the differences come from?
## Wilcoxon-Nemenyi-McDonald-Thompson test (Hollander and Wolfe, 1999, p. 295)
## Note: all pairwise comparisons
(st <- lc("symmetry_test",times ~ methods | block, data = rounding,
                     ytrafo = function(data)
                         trafo(data, numeric_trafo = rank_trafo,
                               block = rounding$block),
                     xtrafo = mcp_trafo(methods = "Tukey")))

## Simultaneous test of all pairwise comparisons
## Wide Angle vs. Round Out differ (Hollander and Wolfe, 1999, p. 296)
pvalue(st, method = "single-step") # subset pivotality is violated


## Strength Index of Cotton
## Hollander and Wolfe (1999, p. 286, Tab. 7.5)
cotton <- data.frame(
    strength = c(7.46, 7.17, 7.76, 8.14, 7.63,
                 7.68, 7.57, 7.73, 8.15, 8.00,
                 7.21, 7.80, 7.74, 7.87, 7.93),
    potash = ordered(rep(c(144, 108, 72, 54, 36), 3),
                     levels = c(144, 108, 72, 54, 36)),
    block = gl(3, 5)
)

## One-sided asymptotic Page test
lc("friedman_test",strength ~ potash | block, data = cotton, alternative = "greater")

## One-sided approximative (Monte Carlo) Page test
lc("friedman_test",strength ~ potash | block, data = cotton, alternative = "greater",
              distribution = approximate(B = 10000))


## Data from Quade (1979, p. 683)
dta <- data.frame(
    y = c(52, 45, 38,
          63, 79, 50,
          45, 57, 39,
          53, 51, 43,
          47, 50, 56,
          62, 72, 49,
          49, 52, 40),
     x = factor(rep(LETTERS[1:3], 7)),
     b = factor(rep(1:7, each = 3))
)

## Approximative (Monte Carlo) Friedman test
## Quade (1979, p. 683)
lc("friedman_test",y ~ x | b, data = dta,
              distribution = approximate(B = 10000)) # chi^2 = 6.000

## Approximative (Monte Carlo) Quade test
## Quade (1979, p. 683)
(qt <- lc("quade_test",y ~ x | b, data = dta,
                  distribution = approximate(B = 10000))) # W = 8.157

## Comparison with R's quade.test() function
quade.test(y ~ x | b, data = dta)

## quade.test() uses an F-statistic
b <- nlevels(qt@statistic@block)
A <- sum(qt@statistic@y^2)
B <- sum(statistic(qt, "linear")^2) / b
(b - 1) * B / (A - B) # F = 8.3765



cleanEx()
nameEx("Transformations")
### * Transformations

flush(stderr()); flush(stdout())

### Name: Transformations
### Title: Functions for Data Transformation
### Aliases: id_trafo rank_trafo normal_trafo median_trafo savage_trafo
###   consal_trafo koziol_trafo klotz_trafo mood_trafo ansari_trafo
###   fligner_trafo logrank_trafo logrank_weight maxstat_trafo
###   fmaxstat_trafo ofmaxstat_trafo f_trafo of_trafo trafo mcp_trafo
### Keywords: manip

### ** Examples

## Dummy matrix, two-sample problem (only one column)
f_trafo(gl(2, 3))

## Dummy matrix, K-sample problem (K columns)
x <- gl(3, 2)
f_trafo(x)

## Score matrix
ox <- as.ordered(x)
of_trafo(ox)
of_trafo(ox, scores = c(1, 3:4))
of_trafo(ox, scores = list(s1 = 1:3, s2 = c(1, 3:4)))

## Normal scores
y <- runif(6)
normal_trafo(y)

## All together now
trafo(data.frame(x = x, ox = ox, y = y), numeric_trafo = normal_trafo)

## The same, but allows for fine-tuning
trafo(data.frame(x = x, ox = ox, y = y), var_trafo = list(y = normal_trafo))

## Transformations for maximally selected statistics
maxstat_trafo(y)
fmaxstat_trafo(x)
ofmaxstat_trafo(ox)

## Apply transformation blockwise (as in the Friedman test)
trafo(data.frame(y = 1:20), numeric_trafo = rank_trafo, block = gl(4, 5))

## Multiple comparisons
dta <- data.frame(x)
mcp_trafo(x = "Tukey")(dta)

## The same, but useful when specific contrasts are desired
K <- rbind("2 - 1" = c(-1,  1, 0),
           "3 - 1" = c(-1,  0, 1),
           "3 - 2" = c( 0, -1, 1))
mcp_trafo(x = K)(dta)



cleanEx()
nameEx("alpha")
### * alpha

flush(stderr()); flush(stdout())

### Name: alpha
### Title: Genetic Components of Alcoholism
### Aliases: alpha
### Keywords: datasets

### ** Examples

## Boxplots
boxplot(elevel ~ alength, data = alpha)

## Asymptotic Kruskal-Wallis test
lc("kruskal_test",elevel ~ alength, data = alpha)



cleanEx()
nameEx("alzheimer")
### * alzheimer

flush(stderr()); flush(stdout())

### Name: alzheimer
### Title: Smoking and Alzheimer's Disease
### Aliases: alzheimer
### Keywords: datasets

### ** Examples

## Spineplots
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:2, ncol = 2))
spineplot(disease ~ smoking, data = alzheimer,
          subset = gender == "Male", main = "Male")
spineplot(disease ~ smoking, data = alzheimer,
          subset = gender == "Female", main = "Female")
par(op) # reset

## Asymptotic Cochran-Mantel-Haenszel test
lc("cmh_test",disease ~ smoking | gender, data = alzheimer)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("asat")
### * asat

flush(stderr()); flush(stdout())

### Name: asat
### Title: Toxicological Study on Female Wistar Rats
### Aliases: asat
### Keywords: datasets

### ** Examples

## Proof-of-safety based on ratio of medians (Pflueger and Hothorn, 2002)
## One-sided exact Wilcoxon-Mann-Whitney test
wt <- lc("wilcox_test",I(log(asat)) ~ group, data = asat,
                  distribution = "exact", alternative = "less",
                  conf.int = TRUE)

## One-sided confidence set
## Note: Safety cannot be concluded since the effect of the compound
##       exceeds 20 % of the control median
exp(confint(wt)$conf.int)



cleanEx()
nameEx("coin-package")
### * coin-package

flush(stderr()); flush(stdout())

### Name: coin-package
### Title: General Information on the 'coin' Package
### Aliases: coin-package coin
### Keywords: package

### ** Examples

## Not run: 
##D ## Generate doxygen documentation if you are interested in the internals:
##D ## Download source package into a temporary directory
##D tmpdir <- tempdir()
##D tgz <- download.packages("coin", destdir = tmpdir, type = "source")[2]
##D ## Extract contents
##D untar(tgz, exdir = tmpdir)
##D ## Run doxygen (assuming it is installed)
##D wd <- setwd(file.path(tmpdir, "coin"))
##D system("doxygen inst/doxygen.cfg")
##D setwd(wd)
##D ## Have fun!
##D browseURL(file.path(tmpdir, "coin", "inst",
##D                     "documentation", "html", "index.html"))
## End(Not run)



cleanEx()
nameEx("expectation-methods")
### * expectation-methods

flush(stderr()); flush(stdout())

### Name: expectation-methods
### Title: Extraction of the Expectation, Variance and Covariance of the
###   Linear Statistic
### Aliases: expectation expectation-methods
###   expectation,IndependenceLinearStatistic-method
###   expectation,IndependenceTest-method variance variance-methods
###   variance,CovarianceMatrix-method
###   variance,IndependenceLinearStatistic-method
###   variance,IndependenceTest-method variance,Variance-method covariance
###   covariance-methods covariance,CovarianceMatrix-method
###   covariance,IndependenceLinearStatistic-method
###   covariance,IndependenceTest-method
### Keywords: methods

### ** Examples

## Example data
dta <- data.frame(
    y = gl(3, 2),
    x = sample(gl(3, 2))
)

## Asymptotic Cochran-Mantel-Haenszel Test
ct <- lc("cmh_test",y ~ x, data = dta)

## The linear statistic, i.e., the contingency table...
(l <- statistic(ct, type = "linear"))

## ...and its expectation...
(El <- expectation(ct))

## ...and covariance
(Vl <- covariance(ct))

## The standardized contingency table...
(l - El) / sqrt(variance(ct))

## ...is identical to the standardized linear statistic
statistic(ct, type = "standardized")



cleanEx()
nameEx("glioma")
### * glioma

flush(stderr()); flush(stdout())

### Name: glioma
### Title: Malignant Glioma Pilot Study
### Aliases: glioma
### Keywords: datasets

### ** Examples

## Grade III glioma
g3 <- subset(glioma, histology == "Grade3")

## Plot Kaplan-Meier estimates
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:2, ncol = 2))
plot(survfit(Surv(time, event) ~ group, data = g3),
     main = "Grade III Glioma", lty = 2:1,
     ylab = "Probability", xlab = "Survival Time in Month",
     xlim = c(-2, 72))
legend("bottomleft", lty = 2:1, c("Control", "Treated"), bty = "n")

## Exact logrank test
lc("logrank_test",Surv(time, event) ~ group, data = g3,
             distribution = "exact")


## Grade IV glioma
gbm <- subset(glioma, histology == "GBM")

## Plot Kaplan-Meier estimates
plot(survfit(Surv(time, event) ~ group, data = gbm),
     main = "Grade IV Glioma", lty = 2:1,
     ylab = "Probability", xlab = "Survival Time in Month",
     xlim = c(-2, 72))
legend("topright", lty = 2:1, c("Control", "Treated"), bty = "n")
par(op) # reset

## Exact logrank test
lc("logrank_test",Surv(time, event) ~ group, data = gbm,
             distribution = "exact")


## Stratified approximative (Monte Carlo) logrank test
lc("logrank_test",Surv(time, event) ~ group | histology, data = glioma,
             distribution = approximate(B = 10000))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("hohnloser")
### * hohnloser

flush(stderr()); flush(stdout())

### Name: hohnloser
### Title: Left Ventricular Ejection Fraction
### Aliases: hohnloser
### Keywords: datasets

### ** Examples

## Asymptotic maximally selected logrank statistics
lc("maxstat_test",Surv(time, event) ~ EF, data = hohnloser)



cleanEx()
nameEx("jobsatisfaction")
### * jobsatisfaction

flush(stderr()); flush(stdout())

### Name: jobsatisfaction
### Title: Income and Job Satisfaction
### Aliases: jobsatisfaction
### Keywords: datasets

### ** Examples

## Approximative (Monte Carlo) linear-by-linear association test
lc("lbl_test",jobsatisfaction, distribution = approximate(B = 10000))



cleanEx()
nameEx("malformations")
### * malformations

flush(stderr()); flush(stdout())

### Name: malformations
### Title: Maternal Drinking and Congenital Sex Organ Malformation
### Aliases: malformations
### Keywords: datasets

### ** Examples

## Graubard and Korn (1987, Tab. 3)

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: midpoint scores (p < 0.05)
midpoints <- c(0, 0.5, 1.5, 4.0, 7.0)
lc("chisq_test",malformation ~ consumption, data = malformations,
           distribution = approximate(B = 1000), alternative = "greater",
           scores = list(consumption = midpoints))

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: midrank scores (p > 0.05)
midranks <- c(8557.5, 24375.5, 32013.0, 32473.0, 32555.5)
lc("chisq_test",malformation ~ consumption, data = malformations,
           distribution = approximate(B = 1000), alternative = "greater",
           scores = list(consumption = midranks))

## One-sided approximative (Monte Carlo) Cochran-Armitage test
## Note: equally spaced scores (p > 0.05)
lc("chisq_test",malformation ~ consumption, data = malformations,
           distribution = approximate(B = 1000), alternative = "greater")



cleanEx()
nameEx("mercuryfish")
### * mercuryfish

flush(stderr()); flush(stdout())

### Name: mercuryfish
### Title: Chromosomal Effects of Mercury-Contaminated Fish Consumption
### Aliases: mercuryfish
### Keywords: datasets

### ** Examples

## Coherence criterion
coherence <- function(data) {
    x <- as.matrix(data)
    matrix(apply(x, 1, function(y)
        sum(colSums(t(x) < y) == ncol(x)) -
            sum(colSums(t(x) > y) == ncol(x))), ncol = 1)
}

## Asymptotic POSET test
poset <- lc("independence_test",mercury + abnormal + ccells ~ group,
                           data = mercuryfish, ytrafo = coherence)

## Linear statistic (T in the notation of Rosenbaum, 1994)
statistic(poset, type = "linear")

## Expectation
expectation(poset)

## Variance
## Note: typo in Rosenbaum (1994, p. 371, Sec. 2, last paragraph)
variance(poset)

## Standardized statistic
statistic(poset)

## P-value
pvalue(poset)

## Exact POSET test
lc("independence_test",mercury + abnormal + ccells ~ group,
                  data = mercuryfish, ytrafo = coherence,
                  distribution = "exact")

## Asymptotic multivariate test
mvtest <- lc("independence_test",mercury + abnormal + ccells ~ group,
                            data = mercuryfish)

## Global p-value
pvalue(mvtest)

## Single-step adjusted p-values
pvalue(mvtest, method = "single-step")

## Step-down adjusted p-values
pvalue(mvtest, method = "step-down")



cleanEx()
nameEx("neuropathy")
### * neuropathy

flush(stderr()); flush(stdout())

### Name: neuropathy
### Title: Acute Painful Diabetic Neuropathy
### Aliases: neuropathy
### Keywords: datasets

### ** Examples

## Conover and Salsburg (1988, Tab. 2)

## One-sided approximative Fisher-Pitman test
lc("oneway_test",pain ~ group, data = neuropathy,
            alternative = "less",
            distribution = approximate(B = 10000))

## One-sided approximative Wilcoxon-Mann-Whitney test
lc("wilcox_test",pain ~ group, data = neuropathy,
            alternative = "less",
            distribution = approximate(B = 10000))

## One-sided approximative Conover-Salsburg test
lc("oneway_test",pain ~ group, data = neuropathy,
            alternative = "less",
            distribution = approximate(B = 10000),
            ytrafo = function(data)
                trafo(data, numeric_trafo = consal_trafo))

## One-sided approximative maximum test for a range of 'a' values
it <- lc("independence_test",pain ~ group, data = neuropathy,
                        alternative = "less",
                        distribution = approximate(B = 10000),
                        ytrafo = function(data)
                            trafo(data, numeric_trafo = function(y)
                                consal_trafo(y, a = 2:7)))
pvalue(it, method = "single-step")



cleanEx()
nameEx("ocarcinoma")
### * ocarcinoma

flush(stderr()); flush(stdout())

### Name: ocarcinoma
### Title: Ovarian Carcinoma
### Aliases: ocarcinoma
### Keywords: datasets

### ** Examples

## Exact logrank test
lt <- lc("logrank_test",Surv(time, event) ~ stadium, data = ocarcinoma,
                   distribution = "exact")

## Test statistic
statistic(lt)

## P-value
pvalue(lt)



cleanEx()
nameEx("photocar")
### * photocar

flush(stderr()); flush(stdout())

### Name: photocar
### Title: Multiple Dosing Photococarcinogenicity Experiment
### Aliases: photocar
### Keywords: datasets

### ** Examples

## Plotting data
op <- par(no.readonly = TRUE) # save current settings
layout(matrix(1:3, ncol = 3))
with(photocar, {
    plot(survfit(Surv(time, event) ~ group),
         lty =  1:3, xmax = 50, main = "Survival Time")
    legend("bottomleft", lty = 1:3, levels(group), bty = "n")
    plot(survfit(Surv(dmin, tumor) ~ group),
         lty = 1:3, xmax = 50, main = "Time to First Tumor")
    legend("bottomleft", lty = 1:3, levels(group), bty = "n")
    boxplot(ntumor ~ group, main = "Number of Tumors")
})
par(op) # reset

## Approximative multivariate (all three responses) test
it <- lc("independence_test",Surv(time, event) + Surv(dmin, tumor) + ntumor ~ group,
                        data = photocar,
                        distribution = approximate(B = 10000))

## Global p-value
pvalue(it)

## Why was the global null hypothesis rejected?
statistic(it, type = "standardized")
pvalue(it, method = "single-step")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("pvalue-methods")
### * pvalue-methods

flush(stderr()); flush(stdout())

### Name: pvalue-methods
### Title: Computation of the p-Value, Mid-p-Value and p-Value Interval
### Aliases: pvalue pvalue-methods pvalue,IndependenceTest-method
###   pvalue,MaxTypeIndependenceTest-method pvalue,NullDistribution-method
###   midpvalue midpvalue-methods midpvalue,IndependenceTest-method
###   midpvalue,NullDistribution-method pvalue_interval
###   pvalue_interval-methods pvalue_interval,IndependenceTest-method
###   pvalue_interval,NullDistribution-method
### Keywords: methods htest

### ** Examples

## Two-sample problem
dta <- data.frame(
    y = rnorm(20),
    x = gl(2, 10)
)

## Exact Ansari-Bradley test
(at <- lc("ansari_test",y ~ x, data = dta, distribution = "exact"))
pvalue(at)
midpvalue(at)
pvalue_interval(at)


## Bivariate two-sample problem
dta2 <- data.frame(
    y1 = rnorm(20) + rep(0:1, each = 10),
    y2 = rnorm(20),
    x = gl(2, 10)
)

## Approximative (Monte Carlo) bivariate Fisher-Pitman test
(it <- lc("independence_test",y1 + y2 ~ x, data = dta2,
                         distribution = approximate(B = 10000)))

## Global p-value
pvalue(it)

## Joint distribution single-step p-values
pvalue(it, method = "single-step")

## Joint distribution step-down p-values
pvalue(it, method = "step-down")

## Sidak step-down p-values
pvalue(it, method = "step-down", distribution = "marginal", type = "Sidak")

## Unadjusted p-values
pvalue(it, method = "unadjusted")


## Length of YOY Gizzard Shad (Hollander and Wolfe, 1999, p. 200, Tab. 6.3)
yoy <- data.frame(
    length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44,
               42, 60, 32, 42, 45, 58, 27, 51, 42, 52,
               38, 33, 26, 25, 28, 28, 26, 27, 27, 27,
               31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
    site = gl(4, 10, labels = as.roman(1:4))
)

## Approximative (Monte Carlo) Fisher-Pitman test with contrasts
## Note: all pairwise comparisons
(it <- lc("independence_test",length ~ site, data = yoy,
                         distribution = approximate(B = 10000),
                         xtrafo = mcp_trafo(site = "Tukey")))

## Joint distribution step-down p-values
pvalue(it, method = "step-down") # subset pivotality is violated



cleanEx()
nameEx("rotarod")
### * rotarod

flush(stderr()); flush(stdout())

### Name: rotarod
### Title: Rotating Rats
### Aliases: rotarod
### Keywords: datasets

### ** Examples

## One-sided exact Wilcoxon-Mann-Whitney test (p = 0.0186)
lc("wilcox_test",time ~ group, data = rotarod, distribution = "exact",
            alternative = "greater")

## Two-sided exact Wilcoxon-Mann-Whitney test (p = 0.0373)
lc("wilcox_test",time ~ group, data = rotarod, distribution = "exact")

## Two-sided asymptotic Wilcoxon-Mann-Whitney test (p = 0.0147)
lc("wilcox_test",time ~ group, data = rotarod)



cleanEx()
nameEx("statistic-methods")
### * statistic-methods

flush(stderr()); flush(stdout())

### Name: statistic-methods
### Title: Extraction of the Test Statistic and Linear Statistic
### Aliases: statistic statistic-methods
###   statistic,IndependenceLinearStatistic-method
###   statistic,IndependenceTest-method
###   statistic,IndependenceTestStatistic-method
### Keywords: methods

### ** Examples

## Example data
dta <- data.frame(
    y = gl(4, 5),
    x = gl(5, 4)
)

## Asymptotic Cochran-Mantel-Haenszel Test
ct <- lc("cmh_test",y ~ x, data = dta)

## Test statistic
statistic(ct)

## The unstandardized linear statistic...
statistic(ct, type = "linear")

## ...is identical to the contingency table
xtabs(~ x + y, data = dta)

## Illustrating departures from the null hypothesis of independence using the
## standardized linear statistic
statistic(ct, type = "standardized")



cleanEx()
nameEx("treepipit")
### * treepipit

flush(stderr()); flush(stdout())

### Name: treepipit
### Title: Tree Pipits in Franconian Oak Forests
### Aliases: treepipit
### Keywords: datasets

### ** Examples

## Asymptotic maximally selected statistics
lc("maxstat_test",counts ~ age + coverstorey + coverregen + meanregen +
                      coniferous + deadtree + cbpiles + ivytree,
             data = treepipit)



cleanEx()
nameEx("vision")
### * vision

flush(stderr()); flush(stdout())

### Name: vision
### Title: Unaided Distance Vision
### Aliases: vision
### Keywords: datasets

### ** Examples

## Asymptotic Stuart(-Maxwell) test (Q = 11.96)
diag(vision) <- 0 # speed-up
lc("mh_test",vision)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
