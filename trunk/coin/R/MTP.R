

### single step maxT multiple testing procedure
singlestep <- function(object, ...) {

    if (object@statistic@alternative == "two.sided") {
        ts <- abs(statistic(object, "standardized"))
    } else {
        ts <- statistic(object, "standardized")
    }
    ### iterate over unique(ts) only and remap
    n <- length(ts)
    idx <- c(which(ts[-1L] != ts[-n]), n)
    uts <- ts[idx]
    ret <- 1 - sapply(uts, pperm, object = object, ...)
    ### HW: Why not 1 - pperm(object, uts) ?
    ret <- matrix(rep.int(ret, diff(c(0L, idx))),
                  nrow = nrow(ts), ncol = ncol(ts))
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}

### algorithm 2.8 (Free Step-Down Resampling Method) in
### Westfall & Young (1993), page 66 _using standardized 
### statistics instead of p-values_!
### <FIXME>
sdmaxT <- function(pls, ts) {

    ### order of original statistics
    rts <- order(ts)

    ### algorithm 2.8 (Free Step-Down Resampling Method) in 
    ### Westfall & Young (1993), page 66 _using standardized
    ### statistics instead of p-values_!
    q <- pls[,rts, drop = FALSE]

    if (ncol(q) > 1) {
        for (j in 2:ncol(q))
            q[,j] <- pmax(q[,j], q[,j-1])
    }
    ret <- rowMeans(GE(t(q), ts[rts]))
    ## enforce monotonicity
    for (i in (length(ret) - 1):1)
        ret[i] <- max(ret[i], ret[i + 1])

    ret <- matrix(ret[rank(ts)], nrow = nrow(ts), ncol = ncol(ts))
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}
    

### stepdown maxT multiple testing procedure
stepdown <- function(object, ...) {

    if (!(extends(class(object), "MaxTypeIndependenceTest") &&
          extends(class(object@distribution), "ApproxNullDistribution")))
        stop(sQuote("object"), " is not of class ",
             sQuote("MaxTypeIndependenceTest"), 
             " or distribution was not approximated via Monte-Carlo") 
 
    ### raw simulation results, scores have been handled already
    pls <- support(object, raw = TRUE)

    ### standardize
    dcov <- sqrt(variance(object))
    expect <- expectation(object) 
    pls <- t((pls - expect) / dcov)
    ts <- statistic(object, "standardized")

    if (object@statistic@alternative == "two.sided") {
        pls <- abs(pls)
        ts <- abs(ts)
    } 
    if (object@statistic@alternative == "less") {
        pls <- -pls
        ts <- -ts
    }

    sdmaxT(pls, ts)
}

### Discrete permutation method (Westfall & Wolfinger, 1997, AmStat 51, 3-8)
discrete <- function(object, method = c("bonferroni", "sidak"), ...) {

   ### <FIXME> this should be possible when the _exact_ marginal
   ### distributions are available
   ### </FIXME>

   if (!(extends(class(object), "MaxTypeIndependenceTest") &&
         extends(class(object@distribution), "ApproxNullDistribution")))
       stop(sQuote("object"), " is not of class ",
            sQuote("MaxTypeIndependenceTest"),
            " or distribution was not approximated via Monte-Carlo")

   method <- match.arg(method)

   alternative <- object@statistic@alternative

   ### raw simulation results, scores have been handled already
   pls <- support(object, raw = TRUE)

   ### standardize
   dcov <- sqrt(variance(object))
   expect <- expectation(object)
   pls <- t((pls - expect) / dcov)
   ts <- statistic(object, "standardized")

   if (object@statistic@alternative == "two.sided") {
       pls <- abs(pls)
       ts <- abs(ts)
   } 
   if (object@statistic@alternative == "less") {
       pls <- -pls
       ts <- -ts
    }

   ### reorder
   rts <- order(ts)
   pls <- pls[, rts, drop = FALSE]

   ### unadjusted p-values
   pvals <- rowMeans(GE(t(pls), ts[rts]))

   foo <- function(x, t) mean(GE(x, t))

   p <- vector(mode = "list", length = ncol(pls))
   for (i in 1:ncol(pls)) {
       ux <- unique(pls[, i])
       p[[i]] <- sapply(ux, foo, x = pls[, i])
   } 
   
   ### discrete adjustment
   if (method == "bonferroni") {
       adjp <- rep.int(0, length(ts))
       for (i in 1:length(pvals)) {
           for (q in 1:length(p)) {
               x <- p[[q]][p[[q]] <= pvals[i]] # Below Eq. 2
               if (length(x) > 0)
                   adjp[i] <- adjp[i] + max(x) # Eq. 4, Bonferroni
           }
       }
       adjp <- pmin(adjp, 1)
   } else {
       adjp <- rep.int(1, length(ts))
       for (i in 1:length(pvals)) {
           for (q in 1:length(p)) {
               x <- p[[q]][p[[q]] <= pvals[i]] # Below Eq. 2
               if (length(x) > 0)
                   adjp[i] <- adjp[i] * (1 - max(x)) # Eq. 2, Sidak
           }
       }
       adjp <- 1 - pmin(adjp, 1)                     
   }   
   ### enforce monotonicity
   for (i in (length(adjp) - 1):1)
       adjp[i] <- max(adjp[i], adjp[i + 1])
   
   ret <- matrix(adjp[rank(ts)], nrow = nrow(ts), ncol = ncol(ts))
   rownames(ret) <- rownames(ts)
   colnames(ret) <- colnames(ts)
   ret
}

#####################################
## Westfall (1997) method in coin ###
#####################################
## Basic code for npmcp() taken from pqfunctions.R in package 
## multcomp

### cf. mcp(x = "Tukey") in multcomp
mcp_trafo <- function(...) {

  args <- list(...)
  stopifnot(length(args) == 1)

  ret <- function(data) {

      x <- data[[names(args)]]
      stopifnot(is.factor(x))
      C <- args[[1]]
      if (is.character(C)) {
          C <- contrMat(table(x), "Tukey")
      } else {
          stopifnot(is.matrix(C))
          stopifnot(ncol(C) == nlevels(x))
      }
      ret <- trafo(data, 
                   factor_trafo = function(x) {   
                       model.matrix(~ x - 1) %*% t(C)
                   })
      attr(ret, "contrast") <- C
      ret
  }
  ret
}

### compute p-values under subset pivotality
npmcp <- function(object) {

    ### extract from object
    y <- object@statistic@y[[1]]
    x <- object@statistic@x[[1]]
    ytrafo <- object@statistic@ytrafo
    alternative <- object@statistic@alternative
    distribution <- object@call$distribution
    stand_tstat <- statistic(object, type = "standardized")
    tstat <- switch(alternative,
                    "less" = stand_tstat,
                    "greater" = -stand_tstat,
                    "two.sided" = -abs(stand_tstat))

    # get contrast matrix from xtrans
    C <- attr(object@statistic@xtrans, "contrast")
    stopifnot(inherits(C, "matrix"))
  
    # order test statistics, most "extreme" one comes first
    Corder <- C[order(tstat), , drop = FALSE]
  
    # compute allowed subsets of hypotheses
    # returns list consisting of lists (one for each rejection step of H0)
    ms <- multcomp:::maxsets(Corder)

    p <- sapply(ms, function(sub) { # for every list of allowed subsets
        max(sapply(sub, function(s) { # for every subset
            Ctmp <- Corder[s, , drop = FALSE] # current allowed subset
            # x levels in current subset
            xlev <- apply(Ctmp, MARGIN = 2, function(col) any(col != 0))
      
            dattmp <- subset(data.frame(y = y, x = x),
                x %in% names(xlev)[xlev]) # relevant data subset
            pvalue(
                independence_test(y ~ x,
                          data = dattmp,
                          xtrafo = mcp_trafo(x = Ctmp),
                          ytrafo = ytrafo, 
                          distribution = distribution,
                          alternative = alternative)
           )
        }))
    })

    for (i in 2:length(p))
        p[i] <- max(p[i-1], p[i]) # forces pvalue monotonicity

    ret <- matrix(p[rank(tstat)])
    attr(ret, "dimnames") <- attr(tstat, "dimnames")
    return(ret)
}

### unadjusted p-values
unadjusted <- function(object, ...) {    
    
    if (extends(class(object@distribution), "AsymptNullDistribution")) {                
        ts <- statistic(object, "standardized")    
        ret <- switch(object@statistic@alternative,
                      "less"      = pnorm(ts),
                      "greater"   = 1 - pnorm(ts),
                      "two.sided" = 2 * pmin(pnorm(ts), 1 - pnorm(ts)))
    } else {
        ### raw simulation results, scores have been handled already
        pls <- support(object, raw = TRUE)

        ## standardize
        dcov <- sqrt(variance(object))
        expect <- expectation(object) 
        pls <- t((pls - expect) / dcov)
        ts <- statistic(object, "standardized")

        if (object@statistic@alternative == "two.sided") {
            pls <- abs(pls)
            ts <- abs(ts)
        } 
        if (object@statistic@alternative == "less") {
            pls <- -pls
            ts <- -ts
        }
    
        ## unadjusted p-values
        ret <- matrix(rowMeans(GE(t(pls), as.vector(ts))),
                      nrow = nrow(ts), ncol = ncol(ts))
        rownames(ret) <- rownames(ts)
        colnames(ret) <- colnames(ts)
    }
    
    ret
}
