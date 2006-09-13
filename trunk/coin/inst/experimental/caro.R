
library("coin")
set.seed(290875)

### bivariates Problem y = N(mu, 1) mit mu = 0 (x < 0.5) oder 1 (x >= 0.5)
x <- runif(100)
y <- (x < 0.5) + rnorm(length(x), sd = 1.5)

### das Uebliche
maxstat_test(y ~ x)

### OK, mal nur drei cutpoints
xtmp <- ordered(cut(x, c(0, 0.25, 0.5, 0.75, 1)) )

### fuer jeden cutpoint "Verrauschen"
y1 <- y + rnorm(length(x), sd = 0.25)
y2 <- y + rnorm(length(x), sd = 0.25)
y3 <- y + rnorm(length(x), sd = 0.25)

### Teststatistik ausrechnen, ACHTUNG: p-Wert ist sinnlos!
mt <- maxstat_test(y1 + y2 + y3 ~ xtmp)
mt 

### nur die Diagonale ist wichtig!
T <- diag(statistic(mt, "standardized"))

### nur die entsprechenden Eintraege der Kovarianzmatrix:
tmp <- rep(0, 3)
index <- as.logical(c(1, tmp, 1, tmp, 1))

### Korrelationsmatrix ausrechnen
cr <- cov2cor(covariance(mt)[index, index])

### neue Teststatistik
Tmax <- max(abs(T))

### und neuer P-Wert
1 - pmvnorm(rep(-Tmax, 3), rep(Tmax, 3), corr = cr)

### Monte Carlo, etwas komplizierter
### resampling, 10000 mal, generiert lineare Statistiken
object <- mt@statistic
B <- 10000
pls <- plsraw <- .Call("R_MonteCarloIndependenceTest", object@xtrans,
                  object@ytrans, as.integer(object@block), as.integer(B), 
                  PACKAGE = "coin")

### standardisieren
dcov <- sqrt(variance(object))
expect <- expectation(object) 
pls <- (pls - expect) / dcov  

### und jetzt nur die interessanten auswaehlen
pls <- pls[index, ]

### p-Wert ausrechnen
mean(colSums(abs(round(pls, 10)) >= round(Tmax, 10)) > 0)
 
