
### identity transformation
id_trafo <- function(x) x

### Ansari-Bradley
ansari_trafo <- function(x) {
    r <- rank(x)
    pmin(r, length(x) - r + 1)
}

### Fligner
fligner_trafo <- function(x) 
    qnorm((1 + rank(abs(x))/(length(x) + 1))/2)

### Normal Scores (van der Waerden)
normal_trafo <- function(x)
    qnorm(rank(x)/(length(x) + 1))

### Median Scores
median_trafo <- function(x)
    ifelse(rank(x) <= (length(x) + 1)/2, 0, 1)

### Conover & Salsburg (1988)
consal_trafo <- function(x)
    (rank(x)/(length(x) + 1))^4

### maximally selected (rank, chi^2, whatsoever) statistics
maxstat_trafo <- function(x, minprob = 0.1, maxprob = 0.9) {
    qx <- quantile(x, prob = c(minprob, maxprob), type = 1)
    ux <- unique(x)
    cutpoints <- ux[ux > qx[1] & ux <= qx[2]] 
    cm <- matrix(unlist(sapply(cutpoints, function(cut)
        as.numeric(x <= cut))), nrow = length(x))
    cm
}

#### Logrank 
logrank_trafo <- function(x) {
    time <- x[,1]
    event <- x[,2]
    n <- length(time)
    ot <- order(time)
    rt <- rank(time, ties.method = "max")
    fact <- event/(n - rt + 1)
    event - cumsum(fact[ot])[rt]
}

### factor handling
f_trafo <- function(x) {
    mm <- model.matrix(~ x - 1)
    colnames(mm) <- levels(x)
    mm <- mm[,colSums(mm) > 0,drop = FALSE]
    ### the two-sample situations
    if (ncol(mm) == 2) mm <- mm[,-2,drop = FALSE]
    return(mm)
}

### transformation function
trafo <- function(data, numeric_trafo = id_trafo, factor_trafo = f_trafo, 
                 surv_trafo = logrank_trafo) {
    tr <- lapply(data, function(x) {
        if (is.factor(x))
            return(factor_trafo(x))
        if (class(x) == "Surv")
            return(surv_trafo(x))
        if (is.numeric(x))
            return(numeric_trafo(x))
        stop("data class ", class(x), " is not supported")
    })
    RET <- c()
    assignvar <- c()
    for (i in 1:length(tr)) {
        RET <- cbind(RET, tr[[i]])
        p <- ifelse(is.matrix(tr[[i]]), ncol(tr[[i]]), 1)
        assignvar <- c(assignvar, rep(i, p))
    }
    attr(RET, "assign") <- assignvar
    return(RET)
}
