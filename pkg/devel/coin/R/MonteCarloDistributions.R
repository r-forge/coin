split_index <- function(n, by) {
    if (n < by)
        by <- n
    lengths(lapply(seq_len(by), function(i) seq.int(i, n, by)),
            use.names = FALSE)
}

MonteCarlo <- function(x, y, block, weights, nresample, standardise = FALSE, parallel, ncpus, cl) {

    montecarlo <- function(nresample) {
        ret <- LinStatExpCov(X = x, Y = y, weights = as.integer(weights),
                              block = as.factor(block), nresample = nresample,
                              standardise = as.integer(standardise))
        if (standardise)
            ret[c("Variance", "PermutedLinearStatistic",
                   "StandardisedPermutedLinearStatistic")]
        ret
    }

    if (parallel == "no")
        montecarlo(nresample)
    else {
        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
            ## advance stream in master process upon exit
            on.exit(assign(".Random.seed",
                           value = nextRNGStream(
                               get(".Random.seed", envir = .GlobalEnv,
                                   inherits = FALSE)),
                           envir = .GlobalEnv))

        if (parallel == "multicore") {
            if (.Platform$OS.type == "windows")
                stop(sQuote(paste0("parallel = ", dQuote("multicore"))),
                     " is not available for MS Windows")
            if (as.integer(ncpus) < 2L)
                warning("parallel operation requires at least two processes")
            ret <- mclapply(split_index(nresample, ncpus),
                            FUN = montecarlo, mc.cores = ncpus)
        } else {
            if (is.null(cl)) {
                ## has a default cluster been registered?
                ## see parallel:::defaultCluster
                ## <FIXME> R-3.5.0 has 'getDefaultCluster()' </FIXME>
                cl <- get("default",
                          envir = get(".reg", envir = getNamespace("parallel"),
                                      inherits = FALSE),
                          inherits = FALSE)
                if (is.null(cl)) {
                    ## no default cluster, so setup a PSOCK cluster
                    cl <- makePSOCKcluster(ncpus)
                    on.exit(stopCluster(cl), add = TRUE) # clean-up
                }
            }
            if (RNGkind()[1L] == "L'Ecuyer-CMRG")
                ## distribute streams (using master process) for reproducibility
                clusterSetRNGStream(cl)
            ncpus <- as.integer(length(cl))
            if (ncpus < 2L)
                warning("parallel operation requires at least two processes")
            ret <- clusterApply(cl, x = split_index(nresample, ncpus),
                                fun = montecarlo)
        }
        ret[[1]]$PermutedLinearStatistic <-
            do.call("cbind", lapply(ret, function(i)
                i$PermutedLinearStatistic))
        if (standardise)
            ret[[1]]$StandardisedPermutedLinearStatistic <-
                do.call("cbind", lapply(ret, function(i)
                    i$StandardisedPermutedLinearStatistic))
        ret[[1]]
    }
}
