twos_exactp_wh <- function(W, S, iavail = NULL, javail = NULL, STAT, ct = 0) { 

    stat <- W %*% S 
    if (stat > STAT) return(NA)

    if (length(iavail) == 0) return(1)
    np <- 1
    i <- max(iavail)

    for (j in javail) {
        Wij <- W
        Wij[,i] <- W[,j]
        Wij[,j] <- W[,i]
        a <- twos_exactp_wh(Wij, S, iavail[iavail != i], 
                            c(i, javail[javail != j]), 
                            STAT, ct + 1)
        if (is.na(a)) break;
        javail <- javail[javail != j]
        np <- np + a
    }
    return(np)
}

twos_exactp <- function(object) {

    if (!extends(class(object), "ScalarIndependenceTestStatistic"))
            stop("Argument ", sQuote("object"), " is not of class ",
                  sQuote("ScalarIndependenceTestStatistic"))

    ### 2 groups as  variable
    groups <- ncol(object@xtrans) == 1 && all(object@xtrans[,1] %in% c(0, 1))

    values <- ncol(object@ytrans) == 1 && is.numeric(object@y[[1]])

    if (!(groups && values))
        stop(sQuote("object"), " does not represent a two-sample statistic")

    alternative <- object@alternative
    if (alternative != "less")
        stop("currently only alternative = less implemented")

    statistic <- drop(object@teststatistic * sqrt(object@covariance) +
                      object@expectation)
    xtrans <- t(object@xtrans)
    ytrans <- object@ytrans

    m <- sum(xtrans)
    n <- ncol(xtrans) - m
    pX <- xtrans[ ,rev(order(xtrans)), drop = FALSE]
    pY <- ytrans[order(ytrans), ,drop = FALSE]
    iavail <- 1:m
    javail <- (m+1):(m+n)
    px <- twos_exact_wh(pX, pY, iavail, javail, statistic, 0)
    (prod(m:1)*prod(n:1) * px) / prod((m+n):1)
}

