
tfamily <- function(name, parm, p, q, d, dd, ddd, dd2d, lp2PI, PI2lp) {

    if (missing(ddd)) ddd <- NA
    if (missing(dd2d)) dd2d <- function(z) dd(z) / d(z)
    if (missing(lp2PI)) lp2PI <- plogis
    if (missing(PI2lp)) PI2lp <- qlogis
    ret <- list(name = name,
                parm = parm,
                p = p,
                q = q,
                d = d,
                dd = dd,
                ddd = ddd,
                dd2d = dd2d,
                lp2PI = lp2PI,
                PI2lp = PI2lp)
    class(ret) <- "tfamily"
    ret
}

Wilcoxon <- function()
    tfamily(name = "proportional odds model", parm = "log-odds ratio",
            p = plogis,
            q = qlogis,
            d = dlogis,
            dd = function(z) {
                p <- plogis(z)
                p * (1 - p)^2 - p^2 * (1 - p)
            },
            ddd = function(z) {
                ex <- exp(z)
                ifelse(is.finite(z), (ex - 4 * ex^2 + ex^3) / (1 + ex)^4, 0.0)
            },
            dd2d = function(z) {
                ex <- exp(z)
                (1 - ex) / (1 + ex)
            },
            lp2PI = function(x) {
               OR <- exp(x)
               ret <- OR * (OR - 1 - x)/(OR - 1)^2
               ret[abs(x) < .Machine$double.eps] <- 0.5
               return(ret)
            },
            PI2lp = function(p) {
               f <- function(x, PI)
                   x + (exp(-x) * (PI + exp(2 * x) * (PI - 1) + exp(x)* (1 - 2 * PI)))
               ret <- sapply(p, function(p) 
                   uniroot(f, PI = p, interval = 50 * c(-1, 1))$root)
               return(ret)
            }
    )

Lehmann <- function()
    tfamily(name = "Lehmann alternative", parm = "log-reverse time hazard ratio",
            p = function(z) exp(-exp(-z)),
            q = function(p) -log(-log(p)),
            d = function(z) ifelse(is.finite(z), exp(- z - exp(-z)), 0.0),
            dd = function(z) {
               ex <- exp(-z)
               ifelse(is.finite(z), exp(-ex - z) * (ex - 1.0), 0.0)
            },
            ddd = function(z) {
               ex <- exp(-z)
               ifelse(is.finite(z), exp(-z - ex) * (ex - 1)^2 - exp(-ex - 2 * z), 0.0)
            },
            dd2d = function(z)
               exp(-z) - 1,
            lp2PI = plogis,
            PI2lp = qlogis
    )

Savage <- function()
    tfamily(name = "proportional hazard model", parm = "log-hazard ratio",
            p = function(z) -expm1(-exp(z)),
            q = function(p) log(-log1p(- p)),
            d = function(z) ifelse(is.finite(z), exp(z - exp(z)), 0.0),
            dd = function(z) {
                ex <- exp(z)
                ifelse(is.finite(z), (ex - ex^2) / exp(ex), 0.0)
            },
            ddd = function(z) {
                ex <- exp(z)
                ifelse(is.finite(z), (ex - 3*ex^2 + ex^3) / exp(ex), 0.0)
            },
            dd2d = function(z)
               1 - exp(z),
            lp2PI = plogis,
            PI2lp = qlogis
    )

vdWaerden <- function()
    tfamily(name = "vdWaeren", parm = "generalised Cohen's d",
            p = pnorm,
            q = qnorm,
            d = dnorm,
            dd = function(z) ifelse(is.finite(z), -dnorm(x = z) * z, 0.0), 
            ddd = function(z) ifelse(is.finite(z), dnorm(x = z) * (z^2 - 1), 0.0),
            dd2d = function(z) -z,
            lp2PI = function(x) pnorm(x / sqrt(2)),
            PI2lp = function(p) qnorm(p) * sqrt(2)
    )

Cauchy <- function()
    tfamily(name = "Cauchy", parm = "good question",
            p = pcauchy,
            q = qcauchy,
            d = dcauchy,
            dd = function(z) 
                ifelse(is.finite(z), - 2 * z / (pi * (z^2 + 1)^2), 0.0),
            ddd = function(z) 
                ifelse(is.finite(z), 8 * z^2 / (pi * (z^2 + 1)^3) - 2 / (pi * (z^2 + 1)^2), 0.0)
    )
