
library("lehmann")
library("lattice")

set.seed(29)

N <- 23
x <- gl(2, N)[-(1:ceiling(N/2))]
y <- rlogis(length(x), location = c(0, 0)[x])

f <- formals(lehmann:::trafo.test.factor)
args <- expand.grid(i = 1:50, 
                    link = eval(f$link),
                    alternative = eval(f$alternative),
                    inference = eval(f$inference),
                    mu = c(0, 1), 
                    nbins = c(0, ceiling(N / 5)),
                    B = c(0, 10000, Inf), 
                    est = NA, cl = NA, cr = NA, stringsAsFactors = FALSE)
args <- subset(args, !(alternative != "two.sided" & inference == "LRatio"))
args <- subset(args, !(B > 0 & inference != "PermScore"))
args <- subset(args, !(mu > 0 & link != "logit"))
args <- subset(args, !(!is.finite(B) & link != "logit" & inference != "PermScore"))

for (i in 1:nrow(args)) {
    y <- rlogis(length(x), location = c(0, 0)[x])
    print(args[i,])
    tr <- try(trafo.test(y = y, x = x, link = args$link[i],
               alternative = args$alternative[i],
               inference = args$inference[i],
               mu = args$mu[i], B = args$B[i]))
    print(tr)
    if (!inherits(tr, "try-error")) {
        args$est[i] <- tr$estimate
        args$cl[i] <- tr$conf.int[1]
        args$cr[i] <- tr$conf.int[2]
    }
}

bwplot(est ~ inference + link | mu, data = args)
