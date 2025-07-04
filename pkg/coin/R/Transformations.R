### compute average scores, see Hajek, Sidak, Sen (page 131ff)
.average_scores <-
function(s, x)
    ave(s, factor(x))

### rank with NAs kept in place
.rank <-
function(x, ties.method = "average") {
    if (ties.method == "mid-ranks")
        ties.method <- "average"
    rank(x, na.last = "keep", ties.method = ties.method)
}

### identity transformation
id_trafo <- function(x) x

### rank transformation
rank_trafo <- function(x, ties.method = c("mid-ranks", "random")) {
    ties.method <- match.arg(ties.method)
    .rank(x, ties.method = ties.method)
}

## Klotz (1962)
klotz_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    switch(ties.method,
        "mid-ranks" = {
            qnorm(rank_trafo(x) / (sum(!is.na(x)) + 1))^2
        },
        "average-scores" = {
            s <- qnorm(rank_trafo(x, ties.method = "random") /
                         (sum(!is.na(x)) + 1))^2
            .average_scores(s, x)
        }
    )
}

## Mood
mood_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    switch(ties.method,
        "mid-ranks" = {
            (rank_trafo(x) - (sum(!is.na(x)) + 1) / 2)^2
        },
        "average-scores" = {
            s <- (rank_trafo(x, ties.method = "random") -
                    (sum(!is.na(x)) + 1) / 2)^2
            .average_scores(s, x)
        }
    )
}

### Ansari-Bradley
ansari_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    switch(ties.method,
        "mid-ranks" = {
            r <- rank_trafo(x)
            pmin.int(r, sum(!is.na(x)) - r + 1)
        },
        "average-scores" = {
            r <- rank_trafo(x, ties.method = "random")
            s <- pmin.int(r, sum(!is.na(x)) - r + 1)
            .average_scores(s, x)
        }
    )
}

### Fligner
fligner_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    switch(ties.method,
        "mid-ranks" = {
            qnorm((1 + rank_trafo(abs(x)) / (sum(!is.na(x)) + 1)) / 2)
        },
        "average-scores" = {
            s <- qnorm((1 + rank_trafo(abs(x), ties.method = "random") /
                          (sum(!is.na(x)) + 1)) / 2)
            .average_scores(s, x)
        }
    )
}

### Normal Scores (van der Waerden)
normal_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    switch(ties.method,
        "mid-ranks" = {
            qnorm(rank_trafo(x) / (sum(!is.na(x)) + 1))
        },
        "average-scores" = {
            s <- qnorm(rank_trafo(x, ties.method = "random") /
                         (sum(!is.na(x)) + 1))
            .average_scores(s, x)
        }
    )
}

### Median Scores
median_trafo <- function(x, mid.score = c("0", "0.5", "1")) {
    ## "0.5" => symmetric median scores (Randles & Wolfe, 1979, pp. 264--266)
    x <- as.numeric(x)
    mid.score <- match.arg(mid.score)
    md <- median(x, na.rm = TRUE)
    scores <- as.numeric(x > md)
    if (mid.score != "0")
        scores[x %EQ% md] <- as.numeric(mid.score)
    scores
}

### Savage scores
savage_trafo <- function(x, ties.method = c("mid-ranks", "average-scores")) {
    ties.method <- match.arg(ties.method)
    switch(ties.method,
        "mid-ranks" = {
            s <- 1 / (sum(!is.na(x)) - .rank(x, ties.method = "min") + 1)
            cumsum(s[order(x)])[.rank(x, ties.method = "max")] - 1
        },
        "average-scores" = {
            o <- order(x)
            s <- 1 / (sum(!is.na(x)) - .rank(x, ties.method = "first") + 1)
            .average_scores(cumsum(s[o])[order(o)], x) - 1
        }
    )
}

### Conover & Salsburg (1988)
consal_trafo <- function(x, ties.method = c("mid-ranks", "average-scores"),
                         a = 5) {
    ties.method <- match.arg(ties.method)

    cs <- function(a) {
        switch(ties.method,
            "mid-ranks" = {
                (rank_trafo(x) / (sum(!is.na(x)) + 1))^(a - 1)
            },
            "average-scores" = {
                 s <- (rank_trafo(x, ties.method = "random") /
                         (sum(!is.na(x)) + 1))^(a - 1)
                 .average_scores(s, x)
            }
        )
    }

    if (length(a) == 1) cs(a)
    else vapply(setNames(a, paste("a =", a)), cs, as.double(x))
}

## Koziol-Nemec (1979, p. 46, eq. 2.6)
koziol_trafo <- function(x, ties.method = c("mid-ranks", "average-scores"),
                         j = 1) {
    ties.method <- match.arg(ties.method)
    switch(ties.method,
        "mid-ranks" = {
            sqrt(2) * cospi(j * rank_trafo(x) / (sum(!is.na(x)) + 1))
        },
        "average-scores" = {
            s <- sqrt(2) * cospi(j * rank_trafo(x, ties.method = "random") /
                                   (sum(!is.na(x)) + 1))
            .average_scores(s, x)
        }
    )
}

### maximally selected (rank, chi^2, whatsoever) statistics
### numeric/ordered factor
maxstat_trafo <- ofmaxstat_trafo <-
    function(x, minprob = 0.1, maxprob = 1 - minprob)
{
    ORDERED <- is.ordered(x)
    if (ORDERED) {
        x0 <- factor(x) # drop unused levels
        lev <- levels(x0)
        x <- as.numeric(x0)
    }
    qx <- quantile(x, probs = c(minprob, maxprob), na.rm = TRUE, names = FALSE,
                   type = 1)
    if (diff(qx) < eps)
        return(NULL)
    cp <- sort(unique(x))
    cp <- cp[-length(cp)]
    cp <- if (mean(x <= qx[2], na.rm = TRUE) <= maxprob)
              cp[cp >= qx[1] & cp <= qx[2]]
          else
              cp[cp >= qx[1] & cp < qx[2]]
    cm <- .Call(R_maxstattrafo, x = as.double(x), cutpoints = as.double(cp))
    dimnames(cm) <- list(
        seq_along(x),
        if (ORDERED) {
            lapply(seq_along(cp), function(i) {
                idx <- lev %in% x0[as.logical(cm[, i])]
                paste0("{",
                       paste0(lev[idx], collapse = ", "),
                       "} vs. {",
                       paste0(lev[!idx], collapse = ", "),
                       "}")
            })
        } else
            paste0("x <= ", round(cp, 3))
    )
    cm
}

### compute index matrix of all 2^(nlevel - 1) possible splits
### code translated from package 'tree'
.fsplits <-
    function(nlevel)
{
    mi <- 2.0^(nlevel - 1L) - 1.0
    index <- matrix(0L, nrow = mi, ncol = nlevel)
    index[, 1L] <- 1L
    for (i in seq_len(mi)) {
        ii <- i - 1L
        for (j in seq_len(nlevel)[-1L]) {
            index[i, j] <- ii %% 2L
            ii <- ii %/% 2L
        }
    }
    storage.mode(index) <- "logical"
    index
}

### set up transformation g(x) for all possible binary splits in an unordered x
fmaxstat_trafo <-
    function(x, minprob = 0.1, maxprob = 1 - minprob)
{
    x <- factor(x) # drop unused levels
    lev <- levels(x)
    cp <- .fsplits(length(lev))
    n_cp <- nrow(cp)
    cm <- matrix(0.0, nrow = length(x), ncol = n_cp)
    nm <- vector(mode = "character", length = n_cp)
    for (i in seq_len(n_cp)) {
        idx <- cp[i, ]
        cm[, i] <- x %in% lev[idx]
        nm[i] <- paste0("{",
                        paste0(lev[idx], collapse = ", "),
                        "} vs. {",
                        paste0(lev[!idx], collapse = ", "),
                        "}")
    }
    cm[is.na(x), ] <- NA
    dimnames(cm) <- list(seq_along(x), nm)
    mn <- colMeans(cm, na.rm = TRUE)
    cm[, (mn %GE% minprob) & (mn %LE% maxprob), drop = FALSE]
}

### weighted logrank scores; with three different methods of handling ties
logrank_trafo <-
    function(x, ties.method = c("mid-ranks", "Hothorn-Lausen", "average-scores"),
             weight = logrank_weight, ...)
{
    ties.method <- match.arg(ties.method)

    if (!(is.Surv(x) && isTRUE(attr(x, "type") == "right")))
        stop(sQuote(deparse(substitute(x))), " is not of class ",
             dQuote("Surv"), " representing right-censored data")

    cc <- complete.cases(x)
    time <- x[cc, 1]
    event <- x[cc, 2]

    n <- length(time)

    if (ties.method == "average-scores") {
        noise <- runif(n, max = min(diff(sort(unique(time)))) / 2)
        time0 <- time
        time <- time - event * noise # break tied events at random
    }

    r <- rank(time, ties.method = if (ties.method != "Hothorn-Lausen") "min"
                                  else "max")
    o <- order(time, event)
    or <- r[o]

    ## number at risk, number of ties and events at the ordered unique times
    n_risk <- n - unique(or) + 1L
    n_ties <- if (ties.method != "Hothorn-Lausen") -diff(c(n_risk, 0L))
              else -diff(c(n - unique(rank(time, ties.method = "min")[o]) + 1L, 0L))
    n_event <- vapply(split(event[o], or), sum, NA_real_, USE.NAMES = FALSE)

    ## index: expands ties and returns in original order
    idx <- rep.int(seq_along(n_ties), n_ties)[r] # => unique(or)[idx] == r

    ## weights
    w <- weight(sort(unique(time)), n_risk, n_event, ...)

    ## weighted logrank scores
    nw <- NCOL(w)
    if (nw == 1L) {
        scores <- rep.int(NA_real_, length(cc))
        scores[cc] <-
            if (ties.method != "average-scores")
                cumsum(w * n_event / n_risk)[idx] - event * w[idx]
            else # average over events only
                .average_scores(
                    cumsum(w * n_event / n_risk)[idx] - event * w[idx],
                    time0 + (1 - event) * noise
                )
    } else {
        scores <- matrix(NA_real_, nrow = length(cc), ncol = nw,
                         dimnames = list(NULL, colnames(w)))
        scores[cc, ] <-
            vapply(seq_len(nw), function(i) {
                if (ties.method != "average-scores")
                    cumsum(w[, i] * n_event / n_risk)[idx] - event * w[idx, i]
                else # average over events only
                    .average_scores(
                        cumsum(w[, i] * n_event / n_risk)[idx] - event * w[idx, i],
                        time0 + (1 - event) * noise
                    )
            }, time)
    }
    scores
}

### some popular logrank weights
logrank_weight <-
function(time, n.risk, n.event,
         type = c("logrank", "Gehan-Breslow", "Tarone-Ware", "Peto-Peto",
                  "Prentice", "Prentice-Marek", "Andersen-Borgan-Gill-Keiding",
                  "Fleming-Harrington", "Gaugler-Kim-Liao", "Self"),
         rho = NULL, gamma = NULL)
{
    type <- match.arg(type)

    ## weight functions
    w <- function(rho, gamma) {
        switch(type,
            "logrank" = { # Mantel (1966), Peto and Peto (1972), Cox (1972)
                rep.int(1L, length(time))
            },
            "Gehan-Breslow" = { # Gehan (1965), Breslow (1970)
                n.risk
            },
            "Tarone-Ware" = { # Tarone and Ware (1977)
                n.risk^rho
            },
            "Peto-Peto" = { # Peto and Peto (1972), Leton and Zuluaga (2001)
                S <- cumprod(1 - n.event / n.risk) # S(t), Kaplan-Meier
                c(1, S[-length(S)]) # S(t-)
            },
            "Prentice" = { # Prentice (1978), Leton and Zuluaga (2001)
                cumprod(n.risk / (n.risk + n.event)) # S(t)
            },
            "Prentice-Marek" = { # Prentice and Marek (1979)
                cumprod(1 - n.event / (n.risk + 1)) # S(t)
            },
            "Andersen-Borgan-Gill-Keiding" = { # Andersen et al (1982)
                S <- cumprod(1 - n.event / (n.risk + 1)) # S(t)
                c(1, S[-length(S)]) * n.risk / (n.risk + 1) # S(t-), pred.
            },
            "Fleming-Harrington" = { # Fleming and Harrington (1991)
                S <- cumprod(1 - n.event / n.risk) # S(t), Kaplan-Meier
                S <- c(1, S[-length(S)]) # S(t-)
                S^rho * (1 - S)^gamma
            },
            "Gaugler-Kim-Liao" = { # Gaugler et al (2007)
                S <- cumprod(1 - n.event / (n.risk + 1)) # S(t)
                S^rho * (1 - S)^gamma
            },
            "Self" = { # Self (1991)
                ## NOTE: this allows for arbitrary follow-up times
                v <- (time - diff(c(0, time)) / 2) / max(time[n.event > 0])
                v^rho * (1 - v)^gamma
            }
        )
    }

    ## set defaults and eliminate 'rho' and 'gamma' when redundant
    if (type == "Tarone-Ware") {
        if (is.null(rho)) rho <- 0.5
        gamma <- NULL
    } else if (type %in% c("Fleming-Harrington", "Gaugler-Kim-Liao", "Self")) {
        if (is.null(rho)) rho <- 0
        if (is.null(gamma)) gamma <- 0
    } else rho <- gamma <- NULL

    ## find rho-gamma combinations, recycle if necessary, and re-assign
    rho_gamma <- suppressWarnings(cbind(rho, gamma)) # no warning on recycling
    if (!is.null(rho)) rho <- rho_gamma[, 1]
    if (!is.null(gamma)) gamma <- rho_gamma[, 2]

    ## weights
    if (length(rho) < 2 && length(gamma) < 2) w(rho, gamma)
    else setColnames(vapply(seq_len(nrow(rho_gamma)),
                            function(i) w(rho[i], gamma[i]), time),
                     ## compute names
                     paste0("rho = ", rho,
                            if (!is.null(gamma)) ", gamma = ", gamma))
}

### factor handling
f_trafo <- function(x) {
    mf <- model.frame(~ x, na.action = na.pass, drop.unused.levels = TRUE)
    if (nlevels(mf$x) == 1)
        stop("can't deal with factors containing only one level")
    ## construct design matrix _without_ intercept
    mm <- model.matrix(~ x - 1, data = mf)
    colnames(mm) <- levels(mf$x)
    ## the two-sample situations
    if (ncol(mm) == 2)
        mm <- mm[, -2, drop = FALSE]
    mm
}

### ordered factors
of_trafo <- function(x, scores = NULL) {
    if (!is.ordered(x))
        warning(sQuote(deparse(substitute(x))), " is not an ordered factor")
    nl <- nlevels(x)
    if (is.null(scores)) {
        scores <- if (!is.null(s <- attr(x, "scores")))
                      s
                  else
                      seq_len(nl)
    }
    ## two-sample case: scores must be normalized for exact p-values
    if (nl == 2L) {
        mn <- min(scores)
        scores <- (scores - mn) / (max(scores) - mn)
    }
    if (!is.list(scores))
        scores <- list(scores)
    if (all(lengths(scores) == nl))
        setRownames(do.call(cbind, scores)[x, , drop = FALSE], seq_along(x))
    else
        stop(sQuote("scores"), " does not match the number of levels")
}

### Zheng (2008)
.ordered_scores <- function(r, s) {
    n <- length(s)
    if (r == 1)      # 'x' has three levels => simplest case!
        matrix(s, nrow = 1, ncol = n)
    else if (n == 1) # => all order-preserving binary partitions
        matrix(s, nrow = r, ncol = 1)
    else {           # 'x' has four or more levels
        s1 <- .ordered_scores(r - 1, s)
        s2 <- .ordered_scores(r, s[-1])
        cbind(rbind(s[1], s1), s2)
    }
}

zheng_trafo <- function(x, increment = 0.1) {
    if (!is.ordered(x))
        warning(sQuote(deparse(substitute(x))), " is not an ordered factor")
    if (increment <= 0 || increment > 1)
        stop(sQuote("increment"),
             " must be greater than 0, but not greater than 1")
    r <- nlevels(x) - 2
    if(r == 0)
        stop(sQuote(deparse(substitute(x))), " has less than three levels")

    ## compute scores
    scores <- rbind(0, .ordered_scores(r, seq.int(0, 1, increment)), 1)

    ## compute colnames
    cn <- format(scores, digits = max(1, min(.ndecimals(increment), 4)),
                 scientific = FALSE)
    cn <- vapply(seq_len(ncol(cn)), function(i)
                     paste0(if (is_ytrafo()) "eta" else "gamma",
                            " = (", paste0(cn[, i], collapse = ", "), ")"),
                 NA_character_)

    setDimnames(scores[x, , drop = FALSE], list(seq_along(x), cn))
}

### transformation function
trafo <- function(data, numeric_trafo = id_trafo, factor_trafo = f_trafo,
                  ordered_trafo = of_trafo, surv_trafo = logrank_trafo,
                  var_trafo = NULL, block = NULL) {

    if (!(is.data.frame(data) || is.list(data)))
        stop(sQuote("data"), " is not of class ",
             dQuote("data.frame"), " or ", dQuote("list"))

### <FIXME> This two-pass procedure for 'block' is *very* expensive
###         for large datasets
    if (!is.null(block)) {
        if (!is.factor(block) || length(block) != nrow(data))
            stop(sQuote("block"), " is not a factor with ",
                 nrow(data), " elements")

        ## need to check dimension of matrix returned by
        ## user supplied functions
        ret <- trafo(data, numeric_trafo, factor_trafo, ordered_trafo, surv_trafo)

        ## apply trafo to each block separately
        for (lev in levels(block)) {
            ret[block == lev, ] <-
                trafo(data[block == lev, , drop = FALSE],
                      numeric_trafo, factor_trafo, ordered_trafo, surv_trafo)
        }
        return(ret)
    }
### </FIXME>

    if (!is.null(var_trafo)) {
        if (!is.list(var_trafo)) stop(sQuote("var_trafo"), " is not a list")
        if (!all(names(var_trafo) %in% names(data)))
            stop("variable(s) ",
                 names(var_trafo)[!(names(var_trafo) %in% names(data))],
                 " not found in ", sQuote("var_trafo"))
    }

    ## compute transformations for each variable
    tr <- vector(mode = "list", length = length(data))
    names(tr) <- names(data)
    for (nm in names(data)) {
        x <- data[[nm]]
        if (nm %in% names(var_trafo))
            tr[[nm]] <- as.matrix(var_trafo[[nm]](x))
        else if (is.ordered(x))
            tr[[nm]] <- as.matrix(ordered_trafo(x))
        else if (is.factor(x) || is.logical(x))
            tr[[nm]] <- as.matrix(factor_trafo(x))
        else if (is.Surv(x))
            tr[[nm]] <- as.matrix(surv_trafo(x))
        else if (is.numeric(x))
            tr[[nm]] <- as.matrix(numeric_trafo(x))
        else {
            if (idx <- inherits(x, "AsIs", TRUE))
                oldClass(x) <- oldClass(x)[-idx]
            stop("data class ", paste(dQuote(class(x)), collapse = ", "),
                 " is not supported")
        }
    }

    ## set up a matrix of transformations
    ## when more than one factor is in play, factor names
    ## _and_ colnames of the corresponding rows are combined by '.'
### <FIXME> Yet another *very* expensive operation for large datasets, due to
###         building up the 'ret' object dynamically
    ret <- c()
    assignvar <- c()
    cn <- c()
    for (i in 1:length(tr)) {
        if (nrow(tr[[i]]) != nrow(data))
            stop("transformation of variable ", names(tr)[i],
                 " are not of length / nrow", nrow(data))
        ret <- cbind(ret, tr[[i]])
        if (is.null(colnames(tr[[i]]))) {
            cn <- c(cn, rep.int("", ncol(tr[[i]])))
        } else {
            cn <- c(cn, paste0(if (length(tr) > 1) "." else "", colnames(tr[[i]])))
        }
        assignvar <- c(assignvar, rep.int(i, ncol(tr[[i]])))
    }
### </FIXME>
    attr(ret, "assign") <- assignvar
    if (length(tr) > 1) {
        colnames(ret) <- paste0(rep.int(names(tr), tabulate(assignvar)), cn)
    } else {
        colnames(ret) <- cn
    }
    ret
}

### multiple comparisons, cf. mcp(x = "Tukey") in multcomp
mcp_trafo <- function(...) {

    stopifnot(...length() == 1)

    ret <- function(data) {
        x <- data[[...names()]]
        stopifnot(is.factor(x))
        C <- ..1
        if (is.character(C)) {
            C <- contrMat(table(x), C)
        } else {
            stopifnot(is.matrix(C))
            stopifnot(ncol(C) == nlevels(x))
            if (is.null(colnames(C)))
                colnames(C) <- levels(x)
            attr(C, "type") <- "User-defined"
            class(C) <- c("contrMat", "matrix")
        }

        ret <- trafo(data, factor_trafo = function(x)
            tcrossprod(model.matrix(~ x - 1, data = model.frame(~ x, na.action = na.pass)), C))
        attr(ret, "contrast") <- C
        ret
    }
    ret
}
