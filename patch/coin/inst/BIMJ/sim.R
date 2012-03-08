

######################## TORSTEN #############################
library("coin")

assocAS <- function(x, scores = 1:ncol(x), ...) {

    if (!is.table(x))
        stop(sQuote("x"), " is not a table.")
    if (length(scores) != ncol(x))
        stop("Length of ", sQuote(scores), " does not match ", sQuote("ncol(x)"))
    xtab <- as.vector(t(x))

    tmpdat <- data.frame(groups = factor(rep(paste("G", 1:nrow(x), sep = ""), rowSums(x))),
                         scores1 = rep(rep(scores, nrow(x)), xtab),
                         scores2 = rep(rep(c(scores[1], scores[2], scores[2]), nrow(x)), xtab),
                         scores3 = rep(rep(c(scores[1], scores[1], scores[2]), nrow(x)), xtab))

    independence_test(scores1 + scores2 + scores3 ~ groups, data = tmpdat, ...)
}

set.seed(17057711)

SIMG=function (sims,R,S,p,f0,f1,f2)
{
#sims=10000
#R=200; S=200; 
#p=0.2; f0=0.1; f1=0.1; f2=0.3
phi=f2*(p**2)+(2*f1)*((1-p)*p)+f0*((1-p)**2)
p0=(f0*(1-p)**2)/phi
p1=(2*f1*(1-p)*p)/phi
p2=(f2*p**2)/phi
q0=((1-f0)*(1-p)**2)/(1-phi)
q1=(2*(1-f1)*(1-p)*p)/(1-phi)
q2=((1-f2)*p**2)/(1-phi)
Probs=list(c(p2,p1,p0), c(q2,q1,q0)) # ansteigend
ntotal=c(R,S)

assocASsim <- function(n = sims, ngroup = ntotal, 
                       probs =Probs,
                       scores = 1:3, ...) {

    g1 <- rmultinom(n, ngroup[1], probs[[1]])
    g2 <- rmultinom(n, ngroup[2], probs[[2]])

    padj <- matrix(0, nrow = n, ncol = length(scores))
    praw <- matrix(0, nrow = n, ncol = length(scores))

    for (i in 1:n) {
    
        x <- as.table(rbind(g1[,i], g2[,i]))
        it <- assocAS(x, scores = scores, ...)

        ### adjusted p-value
        padj[i,] <- pvalue(it, method = "single-step")
        
        ### raw p-value
        praw[i,] <-  pnorm(statistic(it, "standardized"))
    }
    list(padj = padj, praw = praw)
}

tw <- assocASsim(alternative = "less")

### minimum p-value erg=table(apply(tw$padj, 1, which.min))
welchmin=apply(tw$padj, 1, which.min)
nwelc=factor(welchmin, levels=c(1,2,3))
kleinst=apply(tw$padj, 1, min)
gg=kleinst<0.05
fgg=factor(gg, levels=c(TRUE, FALSE))
oma=table(nwelc,fgg)

est1=oma[1,1]/sims
est2=oma[2,1]/sims
est3=oma[3,1]/sims
estMax=length(which(tw$padj[,1]<0.05 | tw$padj[,2]<0.05 |tw$padj[,3]<0.05))/sims
estA=length(which(tw$praw[,1]<0.05))/sims
estR=length(which(tw$praw[,2]<0.05))/sims
estD=length(which(tw$praw[,3]<0.05))/sims
ergeb=c("Max"=estMax,"MaxA"=est1,"MaxR"=est2,"MaxD"=est3,"Add"=estA, "Rec"=estR, "Dom"=estD,
 "p"=p,"f0"=f0,"f1"=f1,"f2"=f2,"cas"=R, "con"=S )
ergeb
}

SIMG(sims=10000,R=200,S=200,p=0.5,f0=0.1,f1=0.1,f2=0.1)
SIMG(sims=10000,R=200,S=200,p=0.5,f0=0.1,f1=0.187,f2=0.187)
SIMG(sims=10000,R=200,S=200,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=10000,R=200,S=200,p=0.5,f0=0.1,f1=0.1,  f2=0.175)





SIMG(sims=1000,R=180,S=220,p=0.5,f0=0.1,f1=0.187,f2=0.187)
SIMG(sims=1000,R=150,S=250,p=0.5,f0=0.1,f1=0.187,f2=0.187)
SIMG(sims=1000,R=100,S=300,p=0.5,f0=0.1,f1=0.187,f2=0.187)
SIMG(sims=1000,R=180,S=220,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=1000,R=150,S=250,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=1000,R=100,S=300,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=1000,R=180,S=220,p=0.5,f0=0.1,f1=0.1,  f2=0.175)
SIMG(sims=1000,R=150,S=250,p=0.5,f0=0.1,f1=0.1,  f2=0.175)
SIMG(sims=1000,R=100,S=300,p=0.5,f0=0.1,f1=0.1,  f2=0.175)


SIMG(sims=1000,R=220,S=180,p=0.5,f0=0.1,f1=0.187,f2=0.187)
SIMG(sims=1000,R=250,S=150,p=0.5,f0=0.1,f1=0.187,f2=0.187)
SIMG(sims=1000,R=300,S=100,p=0.5,f0=0.1,f1=0.187,f2=0.187)
SIMG(sims=1000,R=220,S=180,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=1000,R=250,S=150,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=1000,R=300,S=100,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=1000,R=220,S=180,p=0.5,f0=0.1,f1=0.1,  f2=0.175)
SIMG(sims=1000,R=250,S=150,p=0.5,f0=0.1,f1=0.1,  f2=0.175)
SIMG(sims=1000,R=300,S=100,p=0.5,f0=0.1,f1=0.1,  f2=0.175)



SIMG(sims=10000,R=200,S=200,p=0.5,f0=0.1,f1=0.1,  f2=0.1)
SIMG(sims=10000,R=180,S=220,p=0.5,f0=0.1,f1=0.1,  f2=0.1)
SIMG(sims=10000,R=150,S=250,p=0.5,f0=0.1,f1=0.1,  f2=0.1)
SIMG(sims=10000,R=100,S=300,p=0.5,f0=0.1,f1=0.1,  f2=0.1)
SIMG(sims=10000,R=300,S=100,p=0.5,f0=0.1,f1=0.1,  f2=0.1)

SIMG(sims=10000,R=200,S=200,p=0.5,f0=0.1,f1=0.187,  f2=0.187)
SIMG(sims=10000,R=200,S=200,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=10000,R=200,S=200,p=0.5,f0=0.1,f1=0.1,  f2=0.175)

SIMG(sims=10000,R=100,S=100,p=0.2,f0=0.16,f1=0.26,  f2=0.36)
SIMG(sims=10000,R=100,S=100,p=0.2,f0=0.2,f1=0.2,  f2=0.7)
SIMG(sims=10000,R=100,S=100,p=0.2,f0=0.2,f1=0.35,  f2=0.35)

SIMG(sims=10000,R=300,S=100,p=0.5,f0=0.1,f1=0.187,  f2=0.187)
SIMG(sims=10000,R=300,S=100,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=10000,R=300,S=100,p=0.5,f0=0.1,f1=0.1,  f2=0.175)
SIMG(sims=10000,R=100,S=300,p=0.5,f0=0.1,f1=0.187,  f2=0.187)
SIMG(sims=10000,R=100,S=300,p=0.5,f0=0.1,f1=0.15, f2=0.2)
SIMG(sims=10000,R=100,S=300,p=0.5,f0=0.1,f1=0.1,  f2=0.175)

