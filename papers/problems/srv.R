
library("coin")
data("ovarian")

FH <- function(x, p = 1, q = 1 - p) 
    x$surv^p * (1 - x$surv)^q

TW <- function(x, f = sqrt) f(x$n.risk)

PPP <- function(n) cumprod(n / n + 1)

lg <- function(x, w = function(x) 1, ...) {

    KM <- survfit(x)
    nrisk <- KM$n.risk
    nevent <- KM$n.event
    weights <- w(KM, ...)

    ci <- weights - cumsum(weights * nevent / nrisk)
    Ci <- - cumsum(weights * nevent / nrisk)
    
    ties <-  -diff(c(nrisk, 0))
    ci <- rep(ci, ties)
    Ci <- rep(Ci, ties)

    rt <- rank(x[,1], ties.method = "min")
    event <- x[,2]
    event * ci[rt] + (1 - event) * Ci[rt]
}

yt <- function(data) trafo(data, surv_trafo = lg)
ytPP <- function(data) trafo(data, surv_trafo = function(x) 
    lg(x, w = TW, f = PPP))
ytFH <- function(data) trafo(data, surv_trafo = function(x) 
    lg(x, w = FH, p = 1))



ovarian$rx <- as.factor(ovarian$rx)

survdiff(Surv(futime, fustat) ~ rx, data = ovarian)
survdiff(Surv(futime, fustat) ~ rx, data = ovarian, rho = 1)

independence_test(Surv(futime, fustat) ~ rx, data = ovarian, 
    ytrafo = ytFH)

surv_test(Surv(futime, fustat) ~ rx, data = ovarian) 

independence_test(Surv(futime, fustat) ~ rx, data = ovarian, 
    ytrafo = ytPP)

lung$inst <- as.factor(lung$inst)
lung$pat.karno <- as.factor(lung$pat.karno)
a <- survdiff(Surv(time, status) ~ pat.karno, 
    data=subset(lung))

b <- independence_test(Surv(time, status) ~ pat.karno, 
    data = lung, ytrafo = yt, teststat = "quad")

a$obs - a$exp
statistic(b, "linear")
a$var
covariance(b)

tmp <- data.frame(time = rexp(1000), event = rbinom(1000, 1, 0.5), 
                  group = gl(2, 500))

a <- survdiff(Surv(time, event) ~ group, data = tmp)
b <- independence_test(Surv(time, event) ~ group,
    data = tmp, ytrafo = yt, teststat = "quad")
