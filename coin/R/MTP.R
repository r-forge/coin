

### single step maxT multiple testing procedure
singlestep <- function(object, ...) {

    if (object@statistic@alternative == "two.sided") {
        ts <- abs(statistic(object, "standardized"))
    } else {
        ts <- statistic(object, "standardized")
    }
    ret <- 1 - sapply(ts, pperm, object = object, ...)  
    ret <- matrix(ret, nrow = nrow(ts), ncol = ncol(ts))
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
    if (object@statistic@alternative == "two.sided") {
        pls <- lapply(pls, function(x)
            (abs(x - expect) / dcov)
        )
        ts <- drop(abs(statistic(object, "standardized")))
    } else {
        pls <- lapply(pls, function(x)
            (x - expect) / dcov
        )
        ts <- drop(statistic(object, "standardized"))
    }

    ### order of original statistics
    rts <- order(ts)
    if (object@statistic@alternative == "less") rts <- rev(rts)

    ### algorithm 2.8 (Free Step-Down Resampling Method) in
    ### Westfall & Young (1993), page 66 _using standardized 
    ### statistics instead of p-values_!
    q <- matrix(unlist(pls), nrow = length(pls), 
                byrow = TRUE)[,rts]

    if (object@statistic@alternative == "less") {
        for (j in 2:ncol(q))
            q[,j] <- pmin(q[,j], q[,j-1])
        ret <- matrix(rowMeans(t(q) < ts[rts])[rev(rank(ts))], 
                      nrow = nrow(ts), ncol = ncol(ts))
    } else {
        for (j in 2:ncol(q))
            q[,j] <- pmax(q[,j], q[,j-1])
        ret <- matrix(rowMeans(t(q) >= ts[rts])[rank(ts)], 
                      nrow = nrow(ts), ncol = ncol(ts))
    }
    
    rownames(ret) <- rownames(ts)
    colnames(ret) <- colnames(ts)
    ret
}

### <FIXME>: implement one-sided cases </FIXME>

### Bonferroni permutation method (Westfall & Wolfinger, 1997, AmStat 51, 3-8)
dbonf <- function(object, ...) {

   ### <FIXME> this should be possible when the _exact_ marginal
   ### distributions are available
   ### </FIXME>

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

    if (object@statistic@alternative == "two.sided") {
        pls <- lapply(pls, function(x)
            (abs(x - expect) / dcov)
        )
        ### order of original statistics
        ts <- abs(statistic(object, "standardized"))
    } else {
        pls <- lapply(pls, function(x)
            (x - expect) / dcov
        )
        ts <- statistic(object, "standardized")
    }
   pls <- matrix(unlist(pls), nrow = length(pls), byrow = TRUE)

   ### Bonferroni adjustment (Westfall & Wolfinger, 1997)
   adjp <- rep(0, length(ts))
   for (i in 1:length(ts)) {
       for (q in 1:ncol(pls)) {
           x <- pls[pls[,q] >= ts[i],q]
           if (length(x) > 0)
               adjp[i] <- adjp[i] + mean(pls[,q] >= min(x))
       }
   }
   ret <- matrix(pmin(adjp, 1), nrow = nrow(ts), ncol = ncol(ts))
   rownames(ret) <- rownames(ts)
   colnames(ret) <- colnames(ts)
   ret
}

