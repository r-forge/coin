
### Regression tests for the r x c x K problem, i.e.,
### testing the independence of a factor
### `y' and a factor factor `x' (possibly blocked)

set.seed(290875)
library(coin)
isequal <- coin:::isequal

### generate data: 2 x 2 x K
dat <- data.frame(x = gl(2, 50), y = gl(2, 50)[sample(1:100)], 
                  block = gl(10, 10)[sample(1:100)])[sample(1:100, 75),]

### Pearsons Chisq Test, asymptotic distribution
ptwo <- chisq.test(table(dat$x, dat$y), correct = FALSE)$p.value

stopifnot(isequal(pvalue(chisq_test(y ~ x, data = dat)), ptwo))
stopifnot(isequal(pvalue(chisq_test(table(dat$y, dat$x))), ptwo))

### Cochran-Mantel-Haenzel Test, asymptotic distribution
ptwo <- drop(mantelhaen.test(table(dat$x, dat$y, dat$block), 
                             correct = FALSE)$p.value)

stopifnot(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
stopifnot(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))


### generate data: r x c x K
dat <- data.frame(x = gl(4, 25), y = gl(4, 25)[sample(1:100)], 
                  block = gl(2, 50)[sample(1:100)])

### Cochran-Mantel-Haenzel Test, asymptotic distribution
### _is wrong_ in R < 2.1.0!!!
ptwo <- drop(mantelhaen.test(table(dat$y, dat$x, dat$block), 
                             correct = FALSE)$p.value)

(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))

### generate data: r x c x K
dat <- data.frame(x = gl(4, 25), y = gl(5, 20)[sample(1:100)], 
                  block = gl(2, 50)[sample(1:100)])

### Cochran-Mantel-Haenzel Test, asymptotic distribution
### _is wrong_!!!
ptwo <- drop(mantelhaen.test(table(dat$y, dat$x, dat$block),
                             correct = FALSE)$p.value)

(isequal(pvalue(cmh_test(y ~ x | block, data = dat)), ptwo))
(isequal(pvalue(cmh_test(table(dat$y, dat$x, dat$block))), ptwo))

### Job Satisfaction, Agresti, 2002, Table 7.8, page 288
data(jobsatisfaction, package = "coin")

# both unordered, results in Agresti, 2002, Table 7.12, 3rd row.
stopifnot(isequal(round(pvalue(cmh_test(jobsatisfaction)), 4), 0.3345))

# linear-by-linear association tests
df <- coin:::table2df(jobsatisfaction)
df$Income <- ordered(df$Income)
cmh_test(Job.Satisfaction ~ Income | Gender, data = df)

df$Job.Satisfaction <- ordered(df$Job.Satisfaction)

attr(df$Income, "scores") <- c(3, 10, 20, 35)
attr(df$Job.Satisfaction, "scores") <- c(1, 3, 4, 5)
cmh_test(Job.Satisfaction ~ Income | Gender, data = df)

cmh_test(Job.Satisfaction ~ Income | Gender, data = df, 
         scores = list(Income = c(2500, 10000, 20000, 30000)))
lbl_test(Job.Satisfaction ~ Income | Gender, data = df, 
         scors = list(Income = c(2500, 10000, 20000, 30000)))

### StatXact 6 manual, page 
stopifnot(isequal(round(pvalue(lbl_test(jobsatisfaction, 
                      scores = list(Income = c(3, 10, 20, 35), 
                                    Job.Satisfaction = c(1, 3, 4, 5)))), 
             5), 0.01309))

### congenital sex organ malformation, StatXact 6 manual, page 793
csom <- as.table(matrix(c(17066, 48, 14464, 38, 788, 5, 126, 1, 37, 1), nrow = 2,
    dimnames = list(Malformation = c("Absent", "Present"),
                    Alcohol = c("0", "<1", "1-2", "3-5", ">=6"))))

p <- pvalue(lbl_test(csom))
### asymptotic p-value Cochran-Armitage test (page 796)
stopifnot(isequal(round(p, 4), 0.1764))
stopifnot(isequal(round(pvalue(lbl_test(csom)), 4),
                  round(prop.trend.test(csom[2,], colSums(csom))$p.value, 4)))

### asympototic p-value permutation test with scores (page 807)
p <- pvalue(lbl_test(csom, scores = list(Alcohol = c(0, 0.5, 1.5, 4, 7))))
stopifnot(isequal(round(p, 4), 0.0104))
stopifnot(isequal(round(pvalue(lbl_test(csom, 
                               scores = list(Alcohol = c(0, 0.5, 1.5, 4, 7)))), 4),
                  round(prop.trend.test(csom[2,], colSums(csom),
                        score = c(0, 0.5, 1.5, 4, 7))$p.value, 4)))

### case-control study on oral contraceptives and smoking / myocardial
### infact, StatXact 6 manual, page 797
load("oral_contraceptives.rda")

### standardized test statistic, page 799
stopifnot(isequal(round(sqrt(statistic(lbl_test(oral_contraceptives))), 3), 
          4.335))

### Dose-Response data, StatXact 6 manual, 984
dr <- as.table(matrix(c(100, 18, 50, 50, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1), nrow = 4,
    dimnames = list(dose = paste((1:4)*100, "mg"),
                    tox = c("mild", "moderate", "severe", "death"))))

### with default scores, page 993
stopifnot(isequal(round(pvalue(lbl_test(dr)), 4), 0.0708))

### with alternative scores, page 997
stopifnot(isequal(round(pvalue(lbl_test(dr, scores = list(tox = c(1, 3, 9, 27)))), 4), 
                  0.0828))

### army: with blocks, page 1014
load("army.rda")
stopifnot(isequal(round(pvalue(cmh_test(army)), 6), 0.002774))

