
library("coin")
library("multcomp")

### synchronized permutation scheme (Fortunato Pesarin et al.)
itp2mcp <- function(object, K) {

    ### no weights or blocks
    stopifnot(all(object@weights == 1))
    stopifnot(nlevels(object@block) == 1)

    ### one-way layout (may be relaxed)
    stopifnot(ncol(object@x) == 1)
    grp <- object@x[[1]]
    stopifnot(is.factor(grp))

    ### contrast matrix -> 2 group comparisons only
    stopifnot(is.matrix(K))
    stopifnot(nlevels(grp) == ncol(K))

    ### balanced case only
    stopifnot(length(unique(table(grp))) == 1)

    ### compare two groups at a time
    stopifnot(all(K %in% c(-1, 0, 1)))
    stopifnot(all(rowSums(abs(K)) == 2))

    ### restructure response (synchronization)
    y <- apply(K, 1, function(k) {
        Y1 <- object@y[grp %in% levels(grp)[k == 1], , drop = FALSE]
        Y2 <- object@y[grp %in% levels(grp)[k == -1], , drop = FALSE]
        ### FIXME: second level is baseline -> is alternative correct?
        Y <- rbind(Y1, Y2)
        contr <- paste(levels(grp)[k == 1], levels(grp)[k == -1], sep = "-")
        colnames(Y) <- paste(colnames(Y), contr, sep = "_")
        return(Y)
    })
    y <- do.call("cbind", y)

    ### set up synchronized structure
    x <- gl(2, table(grp)[1])
    ret <- new("IndependenceProblem", x = data.frame(grp = x),
                                      y = y)
    return(ret)
}


### example
### Length of YOY Gizzard Shad from Kokosing Lake, Ohio,
### sampled in Summer 1984, Hollander & Wolfe (1999), Table 6.3, page 200
YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44, 
                             42, 60, 32, 42, 45, 58, 27, 51, 42, 52, 
                             38, 33, 26, 25, 28, 28, 26, 27, 27, 27, 
                             31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                  site = factor(c(rep("I", 10), rep("II", 10),
                                  rep("III", 10), rep("IV", 10))))

### all-pair comparison
sYOY <- new("IndependenceProblem", x = YOY[, "site", drop = FALSE], 
                                   y = YOY[, "length", drop = FALSE])
K <- contrMat(table(YOY$site), "Tukey")
it <- independence_test(itp2mcp(sYOY, K))

### asymptotic p-values
pvalue(it, method = "single-step")

it <- independence_test(itp2mcp(sYOY, K), distribution = approximate(B = 100000))

### asymptotic p-values
pvalue(it, method = "single-step")

