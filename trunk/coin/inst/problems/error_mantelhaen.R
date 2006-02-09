
library("coin")
set.seed(290875)

### a 2 x 2 x K problem
"mydf" <-
structure(list(x = structure(as.integer(c(1, 1, 1, 2, 1, 2, 1,
2, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 2,
2, 2, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1,
2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 2, 2,
2, 2, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2,
2, 1, 2, 2, 2, 1, 2, 2, 2)), .Label = c("1", "2"), class = "factor"),
    y = structure(as.integer(c(1, 1, 1, 2, 2, 1, 2, 1, 2, 2,
    2, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 2,
    1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2,
    1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1,
    2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2,
    2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 1, 2, 1, 2)), .Label = c("1",
    "2"), class = "factor"), z = structure(as.integer(c(1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 10)), .Label = c("1", "2", "3", "4", "5", "6",
    "7", "8", "9", "10"), class = "factor")), .Names = c("x",
"y", "z"), row.names = c("1", "2", "3", "4", "5", "6", "7", "8",
"9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
"20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",
"31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41",
"42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52",
"53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63",
"64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74",
"75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85",
"86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96",
"97", "98", "99", "100"), class = "data.frame")

### mantelhaen.test (exact) vs. coin (Monte-Carlo)
### alt = "less": OK
mantelhaen.test(mydf$x, mydf$y, mydf$z, correct = FALSE, 
                alt = "l", exact = TRUE)$p.value
pvalue(b <- independence_test(y ~ x | z, data = mydf, 
       alt = "l", dist = approximate(B = 9999)))

### alt = "greater": OK
mantelhaen.test(mydf$x, mydf$y, mydf$z, correct = FALSE, 
                alt = "g", exact = TRUE)$p.value
pvalue(b <- independence_test(y ~ x | z, data = mydf, 
       alt = "g", dist = approximate(B = 9999)))

### alt = "two": oups!
mantelhaen.test(mydf$x, mydf$y, mydf$z, correct = FALSE, 
                alt = "t", exact = TRUE)$p.value
pvalue(b <- independence_test(y ~ x | z, data = mydf, 
       alt = "t", dist = approximate(B = 9999)))

expect <- expectation(b)

### jetzt die Innereien von mantelhaen.test
x <- xtabs(~ y + x + z, data = mydf)
m <- apply(x, c(2, 3), sum)[1, ]
n <- apply(x, c(2, 3), sum)[2, ]
t <- apply(x, c(1, 3), sum)[1, ]
s <- sum(x[1, 1, ])
lo <- sum(pmax(0, t - n))
hi <- sum(pmin(m, t))
support <- lo:hi
K <- dim(x)[3]

### das ist die Dichte, ja?
dc <- .C("d2x2xk", as.integer(K), as.double(m), as.double(n), 
         as.double(t), d = double(hi - lo + 1), PACKAGE = "stats")$d

### und so wird der zweiseitige p-Wert ausgerechnet (warum?)
relErr <- 1 + 10^(-7)
d <- dc
sum(d[d <= d[s - lo + 1] * relErr])

### so wuerde ich das machen (weil ich den E-Wert kenne)
sum(d[abs(support - expect) >= (s - expect)])
### schaut deutlich besser aus


