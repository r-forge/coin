
library("coin")

### example
"tab" <-
structure(as.integer(c(15, 10, 14, 11, 10, 15, 11, 14)), .Dim = as.integer(c(2, 
4)), .Dimnames = structure(list(c("1", "2"), c("1", "2", "3", 
"4")), .Names = c("", "")), class = "table")

### severe differences between chisq_test and chisq.test when simulating
### p-values
chisq.test(tab, correct = FALSE)
set.seed(290895)
chisq.test(tab, correct = FALSE, simulate = TRUE, B = 19999)
pvalue(ca <- chisq_test(tab, distribution = approximate(B = 19999)))

### chisq.test Innereien
B <- 19999
x <- tab
n <- sum(x)
nr <- nrow(x)
nc <- ncol(x)
sr <- rowSums(x)
sc <- colSums(x)
E <- outer(sr, sc, "*")/n
dimnames(E) <- dimnames(x)
set.seed(290875)
tmp <- .C("chisqsim", as.integer(nr), as.integer(nc), 
          as.integer(sr), as.integer(sc), as.integer(n), 
          as.integer(B), as.double(E), integer(nr * nc), 
          double(n + 1), integer(nc), results = double(B), 
          PACKAGE = "stats")
s <- round(tmp$results, 10)
stat <- round(chisq.test(tab, correct = FALSE)$statistic, 10)
### same result as in coin!
(1 + sum(s >= stat)) / (length(s) + 1)
pvalue(ca)


### in debug(chisq.test) ausprobieren
#(1 + sum(tmp$results >= STATISTIC))/(B + 1) ### = 0.4
#(1 + sum(round(tmp$results, 10) >= round(STATISTIC, 10)))/(B + 1) ### = 0.5

### Problem: Gleitkommadarstellung!
### d.h. chisq.test berechnet P(S > s) und nicht P(S >= s)
#mean(tmp$results == STATISTIC) ### = 0
#mean(round(tmp$results, 10) == round(STATISTIC, 10)) ### = 0.09
