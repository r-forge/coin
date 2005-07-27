
### stepdown maxT multiple testing procedure
stepdown <- function(object, aggregate = FALSE, MARGIN = 2) {

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
   pls <- lapply(pls, function(x)
       (abs(x - expect) / dcov)
   )

   ### order of original statistics
   ts <- abs(statistic(object, "standardized"))
   if (aggregate) {
       ts <- matrix(apply(ts, MARGIN, max), nrow = 1)
       pls <- lapply(pls, function(x) apply(matrix(x, ncol = length(ts)),
                                            MARGIN, max))
   }
   rts <- order(ts)

   ### algorithm 2.8 (Free Step-Down Resampling Method) in
   ### Westfall & Young (1993), page 66 _using standardized 
   ### statistics instead of p-values_!
   q <- matrix(unlist(pls), nrow = length(pls), 
               byrow = TRUE)[,rts]

   for (j in 2:ncol(q))
       q[,j] <- pmax(q[,j], q[,j-1])
    
   ret <- matrix(rowMeans(t(q) >= ts[rts])[rank(ts)], 
                 nrow = nrow(ts), ncol = ncol(ts))
   rownames(ret) <- rownames(ts)
   colnames(ret) <- colnames(ts)
   ret
}

