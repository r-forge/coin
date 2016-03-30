
set.seed(29)

dyn.load("table.so")

i2w_1 <- function(indx1, indx2)
    .Call("R_int_table", indx1, nlevels(indx1) + 1L, 
                         indx2, nlevels(indx2) + 1L)

i2w_w <- function(indx1, indx2, weights)
    .Call("R_int_w_table", indx1, nlevels(indx1) + 1L, 
                           indx2, nlevels(indx2) + 1L, weights)

i2w_s <- function(indx1, indx2, subset)
    .Call("R_int_s_table", indx1, nlevels(indx1) + 1L, 
                           indx2, nlevels(indx2) + 1L, subset - 1L)

i2w_s_w <- function(indx1, indx2, subset, weights)
    .Call("R_int_s_w_table", indx1, nlevels(indx1) + 1L, 
                             indx2, nlevels(indx2) + 1L, subset - 1L, weights)

i2w <- function(indx1, indx2, weights, subset, perm = FALSE) { 
    if (!perm) {
        if (missing(subset) & missing(weights))
            return(i2w_1(indx1, indx2))
        if (missing(subset) & !missing(weights))
            return(i2w_w(indx1, indx2, weights))
        if (!missing(subset) & missing(weights))
            return(i2w_s(indx1, indx2, subset))
        return(i2w_s_w(indx1, indx2, subset, weights))
    }
    if (missing(subset))
        subset <- 1:length(indx1)
    if (missing(weights))
        weights <- rep(1, length(indx1))
    indx1 <- rep(indx1[subset], weights[subset])
    indx2 <- rep(indx2[subset], weights[subset])
    i2w_1(sample(indx1), indx2)
}

fun <- function() {
n <- 10000
i1 <- factor(sample(1:7, size = n, replace = TRUE))
i2 <- factor(sample(1:5, size = n, replace = TRUE))
w <- as.integer(floor(runif(n, max = 7)))
sum(w)
s <- sample(1:n, size = n/10, replace = FALSE)

(X <- i2w(i1, i2))
(Xp <- i2w(i1, i2, perm = TRUE))
o1 <- all.equal(rowSums(X), rowSums(Xp)) & all.equal(colSums(X), colSums(Xp))

(X <- i2w(i1, i2, subset = s))
(Xp <- i2w(i1, i2, subset = s, perm = TRUE))
o2 <- all.equal(rowSums(X), rowSums(Xp)) & all.equal(colSums(X), colSums(Xp))

(X <- i2w(i1, i2, weights = w))
(Xp <- i2w(i1, i2, weights = w, perm = TRUE))
o3 <- all.equal(rowSums(X), rowSums(Xp)) & all.equal(colSums(X), colSums(Xp))

(X <- i2w(i1, i2, weights = w, subset = s))
(Xp <- i2w(i1, i2, weights = w, subset = s, perm = TRUE))
o4 <- all.equal(rowSums(X), rowSums(Xp)) & all.equal(colSums(X), colSums(Xp))

o1 & o2 & o3 & o4
}

Rprof("a")
all(replicate(100, fun()))
Rprof(NULL)

