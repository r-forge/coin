\documentclass[a4paper]{report}

%\VignetteIndexEntry{Stratified K Sample Inference}
%\VignetteDepends{strKsam}
%\VignetteKeywords{conditional inference, conditional Monte Carlo}}
%\VignettePackage{strKsam}

%% packages
\usepackage{amsfonts,amstext,amsmath,amssymb,amsthm}

\usepackage[utf8]{inputenc}

\newif\ifshowcode
\showcodetrue

\usepackage{latexsym}
%\usepackage{html}

\usepackage{listings}

\usepackage{color}
\definecolor{linkcolor}{rgb}{0, 0, 0.7}



\usepackage[round]{natbib}


\usepackage[%
backref,%
pageanchor=false,%
raiselinks,%
pdfhighlight=/O,%
pagebackref,%
hyperfigures,%
breaklinks,%
colorlinks,%
pdfpagemode=UseNone,%
pdfstartview=FitBH,%
linkcolor={linkcolor},%
anchorcolor={linkcolor},%
citecolor={linkcolor},%
filecolor={linkcolor},%
menucolor={linkcolor},%
urlcolor={linkcolor}%
]{hyperref}

%%% ATTENTION: no bib keys with _ allowed!
\usepackage{underscore}

\usepackage[top=25mm,bottom=25mm,left=25mm,right=25mm]{geometry}

\usepackage{lmodern}

\newcommand{\pkg}[1]{\textbf{#1}}
\newcommand{\proglang}[1]{\textsf{#1}}
\newcommand{\code}[1]{\texttt{#1}}
\newcommand{\file}[1]{\texttt{#1}}

\newcommand{\R}{\mathbb{R} }
\newcommand{\samY}{\mathcal{Y} }
\newcommand{\Prob}{\mathbb{P} }
\newcommand{\N}{\mathbb{N} }
%\newcommand{\C}{\mathbb{C} }
\newcommand{\V}{\mathbb{V}} %% cal{\mbox{\textnormal{Var}}} }
\newcommand{\E}{\mathbb{E}} %%mathcal{\mbox{\textnormal{E}}} }
\newcommand{\Var}{\mathbb{V}} %%mathcal{\mbox{\textnormal{Var}}} }
\newcommand{\argmin}{\operatorname{argmin}\displaylimits}
\newcommand{\argmax}{\operatorname{argmax}\displaylimits}
\newcommand{\LS}{\mathcal{L}_n}
\newcommand{\TS}{\mathcal{T}_n}
\newcommand{\LSc}{\mathcal{L}_{\text{comb},n}}
\newcommand{\LSbc}{\mathcal{L}^*_{\text{comb},n}}
\newcommand{\F}{\mathcal{F}}
\newcommand{\A}{\mathcal{A}}
\newcommand{\yn}{y_{\text{new}}}
\newcommand{\z}{\mathbf{z}}
\newcommand{\X}{\mathbf{X}}
\newcommand{\Z}{\mathbf{Z}}
\newcommand{\Y}{\mathbf{Y}}
\newcommand{\mH}{\mathbf{H}}
\newcommand{\mA}{\mathbf{A}}
\newcommand{\mL}{\mathbf{L}}
\newcommand{\mU}{\mathbf{U}}
\newcommand{\sX}{\mathcal{X}}
\newcommand{\sY}{\mathcal{Y}}
\newcommand{\T}{\mathbf{T}}
\newcommand{\x}{\mathbf{x}}
\renewcommand{\a}{\mathbf{a}}
\newcommand{\xn}{\mathbf{x}_{\text{new}}}
\newcommand{\y}{\mathbf{y}}
\newcommand{\uvec}{\mathbf{u}}
\newcommand{\vvec}{\mathbf{v}}
\newcommand{\w}{\mathbf{w}}
\newcommand{\sbullet}{\mathbin{\vcenter{\hbox{\scalebox{0.5}{$\bullet$}}}}}
\newcommand{\wdot}{\mathbf{w}_{\sbullet}}
\renewcommand{\t}{\mathbf{t}}
\newcommand{\M}{\mathbf{M}}
\renewcommand{\vec}{\text{vec}}
\newcommand{\B}{\mathbf{B}}
\newcommand{\K}{\mathbf{K}}
\newcommand{\W}{\mathbf{W}}
\newcommand{\D}{\mathbf{D}}
\newcommand{\I}{\mathbf{I}}
\newcommand{\bS}{\mathbf{S}}
\newcommand{\cellx}{\pi_n[\x]}
\newcommand{\partn}{\pi_n(\mathcal{L}_n)}
\newcommand{\err}{\text{Err}}
\newcommand{\ea}{\widehat{\text{Err}}^{(a)}}
\newcommand{\ecv}{\widehat{\text{Err}}^{(cv1)}}
\newcommand{\ecvten}{\widehat{\text{Err}}^{(cv10)}}
\newcommand{\eone}{\widehat{\text{Err}}^{(1)}}
\newcommand{\eplus}{\widehat{\text{Err}}^{(.632+)}}
\newcommand{\eoob}{\widehat{\text{Err}}^{(oob)}}
\newcommand{\mub}{\boldsymbol{\mu}}
\newcommand{\Sigmab}{\boldsymbol{\Sigma}}
\def \thetavec        {\text{\boldmath$\theta$}}
\newcommand{\rT}{T}
\newcommand{\rt}{t}


\author{Torsten Hothorn \\ Universit\"at Z\"urich}

\title{Distribution-free Stratified $K$ Sample Inference}

\begin{document}

\pagenumbering{roman}
\maketitle
\tableofcontents

\chapter{Model and Parameterisation}
\label{ch:model}
\pagenumbering{arabic}

Treatment group $\rT \in \{1, \dots, K\}, K \ge 2$, outcome $Y \in \samY$ at least ordered,
stratum $S \in \{1, \dots, B\}$ in $B \ge 1$ blocks with conditional
cumulative distribution function (cdf)
$F_Y(y \mid S = b, \rT = k) = \Prob(Y \le y \mid S = b, \rT = k)$.

\paragraph{Model}

With $g: [0,1] \times \R \rightarrow [0,1]$, we model
$F(y \mid S = b, \rT = k) = g(F(y \mid S = b, \rT = 1), \beta_k)$ for all $b = 1,
\dots, B$, $k = 2, \dots, K$, and $y \in \samY$ based on parameters
$\beta_2, \dots, \beta_K \in \R$. For notational convenience: $\beta_1 := 0$. For
example, $g_\text{L}(p, \beta) = p^{\exp(-\beta)}$ (Lehmann alternative),
$g_\text{PH}(p, \beta) = 1 - (1 -
p)^{\exp(-\beta)}$ proportional hazards, or
$g_\text{PO}(p, \beta) = \text{expit}(\text{logit}(p) - \beta)$ proportional odds, or
$g_\text{Cd}(p, \beta) =
\Phi(\Phi^{-1}(p) - \beta)$ (generalised Cohen's $d$).

For some absolute continuous cdf $F$ with log-concave density $f = F^\prime$
and corresponding derivative $f^\prime$, we write
$$F_Y(y \mid S = b, \rT = k) = F\left(F^{-1}\left(F_Y(y \mid S = b, \rT =
1)\right) - \beta_k\right), \quad k = 2, \dots, K.$$
The negative shift term ensures that positive values of $\beta_k$ correspond
to the situation of outcomes being stochastically larger in group $k$ than
in group one.

Note that $F(z) = \exp(-\exp(-z))$ gives rise to $g_\text{L}$, 
$F(z) = 1 - \exp(-\exp(z))$ to $g_\text{PH}$, $F = \text{expit}$ to
$g_\text{PO}$, and $F = \Phi$ to $g_\text{Cd}$.

\paragraph{Hypthesis}

We are interested in inference for $\beta_2, \dots, \beta_K$, in terms of
confidence intervals and hypothesis tests of the form
\begin{eqnarray*}
& & H_0: \beta_k - \mu_k = 0, \text{``two.sided''}, \quad k = 2, \dots, K, \\
& & H_0: \beta_k - \mu_k \ge 0, \text{``less''}, \quad k = K = 2, \\
& & H_0: \beta_k - \mu_k \le 0, \text{``greater''}, \quad k = K = 2,
\end{eqnarray*}
with the latter two options only for the two-sample case.

\paragraph{Likelihood}

For an ordered categorical outcome $Y$ from sample space $\samY = \{y_1 < y_2 < \cdots <
y_C\}$, we parameterise the model in terms of intercept ($\vartheta_\cdot$) and
shift ($\beta_\cdot$) parameters

$$F_Y(y_c \mid S = b, \rT = k) = F(\vartheta_{c,b} - \beta_k), \quad c = 1, \dots, C.$$

The $C - 1$ intercept parameters are block-specific and monotone increasing
$\vartheta_{0,b} = -\infty < \vartheta_{1,b} < \cdots < \vartheta_{C,b} = \infty$
within each block $b = 1, \dots, B$.
We collect all parameters in a vector
\begin{eqnarray*}
\thetavec = (\theta_1 & := & \beta_2, \\
               & \dots & , \\
               \theta_{K - 1} & := & \beta_K, \\
               \theta_{K} & := & \vartheta_{1,1}, \\
               \theta_{K + 1} & := & \vartheta_{2,1} - \vartheta_{1,1} > 0, \\
               &  \dots, & \\
               \theta_{K + C - 2} & := & \vartheta_{C-1,1} - \vartheta_{C-2,1} > 0, \\
               \theta_{K + C - 1} & := & \vartheta_{1,2}, \\
               & \dots &, \\
               \theta_{B (C - 1) + K - 1} & := & \vartheta_{C-1,B} - \vartheta_{C-2,B} >
               0)
\end{eqnarray*}

For the $i$th observation $(y_i = y_c, s_i = b, \rt_i = k)$ from block $b$
under treatment $k$, the log-likelihood contribution is
$$\log(\Prob(y_{c - 1} < Y \le y_c \mid S = b, \rT = k)) = \log(F(\vartheta_{c,b} - \beta_k) - F(\vartheta_{c - 1,b} - \beta_k)).$$

For an absolutely continuous outcome $Y \in \R$, we define $y_c := y_{(c)}$,
the $c$th distinct ordered observation in the sample. The log-likelihood
above is then the empirical or nonparametric log-likelihood.

We represent the data in form of a $C \times K \times B$ contingency table,
whose element $(c, k, b)$ is the number of observations with configuration $(y = y_c, s = b,
\rt = k)$.
	
\chapter{Parameter Estimation}
\label{ch:est}

@o strKsam_src.R -cp
@{
@<cumsumrev@>
@<negative logLik@>
@<negative score@>
@<negative score residuals@>
@<Hessian@>
@<stratified negative logLik@>
@<stratified negative score@>
@<stratified Hessian@>
@<stratified negative score residual@>
@<R wcrossprod@>
@}


We start implementing the log-likelihood function for parameters \code{parm}
$= \thetavec$ (assuming only a single block) with data from a two-way $C
\times K$ contingency table \code{x}. 

From $\thetavec$, we first extract the shift parameters $\beta_\cdot$ and
then the intercept parameters $\vartheta_\cdot$, compute the differences
$\vartheta_{c,1} - \beta_k$ and evaluate the probabilities
\code{prb} $ = \Prob(y_{c - 1} < Y \le y_c \mid S = 1, \rT = k)$ for all
groups:

@d parm to prob
@{
bidx <- seq_len(ncol(x) - 1L)
beta <- c(0, mu + parm[bidx])
intercepts <- c(-Inf, cumsum(parm[- bidx]), Inf)
tmb <- intercepts - matrix(beta, nrow = length(intercepts),  
                                 ncol = ncol(x),
                                 byrow = TRUE)
Ftmb <- F(tmb)
prb <- Ftmb[- 1L, , drop = FALSE] - 
       Ftmb[- nrow(Ftmb), , drop = FALSE]
@}

With default null values $\mu_k = 0, k = 2, \dots, K$, we define the
negative log-likelihood function as the weighted (by number of observations) sum of
the log-probabilities

@d negative logLik
@{
.nll <- function(parm, x, mu = 0) {
    @<parm to prob@>
    - sum(x * log(prb))
}
@}

It is important to note that, with $F$ corresponding to distribution with
log-concave density $f$, the negative log-likelihood is a convex function of
the parameters $\thetavec$, and thus we can solve the corresponding
constrained minimisation problem quickly and reliably.

Next, we implement the gradient of the negative
log-likelihood, the negative score function for the parameters in
$\thetavec$. The score function for the empirical likelihood, evaluated at
parameters $\vartheta_\cdot$ and $\beta_\cdot$ is given in many places
\citep[for example in][Formula~(2)]{HothornMoestBuehlmann2017}. 
We begin computing the ratio of $f(\vartheta_{c,1} -
\beta_k)$ and the corresponding likelihood

@d d p ratio
@{
ftmb <- f(tmb)
zu <- x * ftmb[- 1, , drop = FALSE] / prb
zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
@}

and compute the negative score function

@d negative score
@{
.nsc <- function(parm, x, mu = 0) {
    @<parm to prob@>
    @<d p ratio@>

    ret <- numeric(length(parm))
    ret[bidx] <- colSums(zl)[-1L] -
                 colSums(zu[-nrow(zu),,drop = FALSE])[-1L]
    ret[-bidx] <- Reduce("+", 
                         lapply(1:ncol(x), 
                             function(j) {
                                 .rcr(zu[-nrow(zu),j]) - 
                                 .rcr(zl[-1,j])
                             })
                         )
    - ret
}
@}

Adjustment for the parameterisation in terms of differences between
intercepts needs this small helper function

@d cumsumrev
@{
.rcr <- function(z)
    Reduce('+', z, accumulate = TRUE, right = TRUE)
    ### rev(cumsum(rev(z)))
@}

In addition, we define negative score residuals, that is, the derivative of the
negative log-likelihood with respect to an intercept term constrained to
zero:

@d negative score residuals
@{
.nsr <- function(parm, x, mu = 0) {
    @<parm to prob@>
    @<d p ratio@>

    rowSums(zu - zl) / rowSums(x)
}
@}

We also need access to the observed Fisher information of the shift
parameters. We proceed by implementing the Hessian for the intercept
($\vartheta_\cdot$) and shift ($\beta_\cdot$) parameters, as given in Formula~(4) of
\cite{HothornMoestBuehlmann2017} first. This partitioned matrix
\begin{eqnarray*}
\mH(\vartheta_1, \dots, \vartheta_{C - 1}, \beta_2, \dots, \beta_K) = 
\left(\begin{array}{ll}
\mA & \X \\
\X^\top & \Z
\end{array} \right)
\end{eqnarray*}
consists of a symmetric tridiagonal $\mA \sim (C-1,C-1)$, a diagonal $\Z \sim (K - 1, K -
1)$, and a full $\X \sim (C - 1, K - 1)$ matrix. In a second step, we
compute the Fisher information matrix for the shift parameters only by means
of the Schur complement $\Z - \X^\top \mA^{-1} \X$.

In addition to probabilities \code{prb}, the Hessian necessitates the
computation of $f(\vartheta_{c,1} - \beta_k)$ and $f^\prime(\vartheta_{c,1} -
\beta_k)$. We start preparing these objects


@d Hessian prep
@{
ftmb <- f(tmb)
fptmb <- fp(tmb)

dl <- ftmb[- nrow(ftmb), , drop = FALSE]
du <- ftmb[- 1, , drop = FALSE]
dpl <- fptmb[- nrow(ftmb), , drop = FALSE]
dpu <- fptmb[- 1, , drop = FALSE]
dlm1 <- dl[,-1L, drop = FALSE]
dum1 <- du[,-1L, drop = FALSE]
dplm1 <- dpl[,-1L, drop = FALSE]
dpum1 <- dpu[,-1L, drop = FALSE]
prbm1 <- prb[,-1L, drop = FALSE]

i1 <- length(intercepts) - 1L
i2 <- 1L
@}

The off-diagonal elements of $\mA$ are now available as
@d Aoffdiag
@{
Aoffdiag <- -rowSums(x * du * dl / prb^2)[-i2]
Aoffdiag <- Aoffdiag[-length(Aoffdiag)]
@}

and the diagonal elements of $\mA$ as
@d Adiag
@{
Adiag <- -rowSums((x * dpu / prb)[-i1,,drop = FALSE] - 
                  (x * dpl / prb)[-i2,,drop = FALSE] - 
                  ((x * du^2 / prb^2)[-i1,,drop = FALSE] + 
                   (x * dl^2 / prb^2)[-i2,,drop = FALSE]
                  )
                 )
                  
@}

For the computation of $\X$ and $\Z$, the observations corresponding to the
control group ($k = 1$) are irrelevant, we remove these first

@d X and Z
@{
xm1 <- x[,-1L,drop = FALSE] 
X <- ((xm1 * dpum1 / prbm1)[-i1,,drop = FALSE] - 
      (xm1 * dplm1 / prbm1)[-i2,,drop = FALSE] - 
      ((xm1 * dum1^2 / prbm1^2)[-i1,,drop = FALSE] - 
       (xm1 * dum1 * dlm1 / prbm1^2)[-i2,,drop = FALSE] -
       (xm1 * dum1 * dlm1 / prbm1^2)[-i1,,drop = FALSE] +
       (xm1 * dlm1^2 / prbm1^2)[-i2,,drop = FALSE]
      )
     )

Z <- -colSums(xm1 * (dpum1 / prbm1 - 
                     dplm1 / prbm1 -
                     (dum1^2 / prbm1^2 - 
                      2 * dum1 * dlm1 / prbm1^2 +
                      dlm1^2 / prbm1^2
                     )
                    )
             )
if (length(Z) > 1L) Z <- diag(Z)
@}

We return the Fisher information for $\beta_2, \dots, \beta_K$ as the Schur
complement $\Z - \X^\top \mA^{-1} \X$ by means of a weighted crossproduct
implemented in Chapter~\ref{ch:schur}.

@d Hessian
@{
.hes <- function(parm, x, mu = 0) {
    @<parm to prob@>

    @<Hessian prep@>

    @<Aoffdiag@>
    @<Adiag@>
    @<X and Z@>

    return(Z - wcrossprod(x = X, list(a = Adiag, b = Aoffdiag)))
}
@}

We start with an example involving $K = 3$ groups for a binary outcome and
use a binary logistic regression model to estimate the two log-odds ratios
$\beta_2$ and $\beta_3$ along with their estimated covariance
<<glm>>=
source("strKsam_src.R")
source("ML.R")
source("linkfun.R")
dyn.load("Schur_src.so")
(x <- matrix(c(10, 5, 7, 11, 8, 9), nrow = 2))
d <- expand.grid(y = relevel(gl(2, 1), "2"), t = gl(3, 1))
d$x <- c(x)
m <- glm(y ~ t, data = d, weights = x, family = binomial())
logLik(m)
(cf <- coef(m))
vcov(m)[-1,-1]
@@

Replicating these results requires specification of the inverse link
function $F = \text{expit}$ and the density function $f$ of the standard
logistic.

<<glm-op>>=
F <- plogis
f <- dlogis
(op <- optim(par = c("mt2" = 0, "mt3" = 0, "(Intercept)" = 0), 
             fn = .nll, gr = .nsc, 
             x = x, method = "BFGS", hessian = TRUE))
c(cf[-1] * -1, cf[1]) - op$par
logLik(m) + op$value
.nsr(op$par, x)
obj <- .free1wayML(as.table(x), link = logit())
obj$coefficients
obj$value
@@

Parameter estimates and the in-sample log-likelihood are practically
identical. We now turn to the inverse Hessian of the shift terms, first
defining the derivative of the density of the standard logistic distribtion
<<glm-H>>=
fp <- function(x) {
    p <- plogis(x)
    p * (1 - p)^2 - p^2 * (1 - p)
}
solve(H <- .hes(op$par, x))
vcov(m)[-1,-1]
solve(op$hessian)[1:2,1:2]
obj$vcov
@@
Also here we see practically identical results.

In the next step, we extend our results to the stratified case. We iterate
over all blocks and evaluate the negative log-likelihood for the same values
of the shift parameters but block-specific values of the intercept
parameters.

@d stratum prep
@{
if (is.table(x)) {
    C <- dim(x)[1]
    K <- dim(x)[2]
    B <- dim(x)[3]
    sidx <- gl(B, C - 1)
    x <- lapply(seq_len(B), function(b) x[,,b,drop = TRUE])
} else {
    C <- sapply(x, nrow)
    K <- unique(sapply(x, ncol))
    stopifnot(length(K) == 1L)
    B <- length(x)
    sidx <- factor(rep(seq_len(B), times = C - 1L), levels = seq_len(B))
}
bidx <- seq_len(K - 1L)
beta <- parm[bidx]
intercepts <- split(parm[-bidx], sidx)
@}

@d stratified negative logLik
@{
.snll <- function(parm, x, mu = 0) {
    @<stratum prep@>
    ret <- 0
    for (b in seq_len(B))
        ret <- ret + .nll(c(beta, intercepts[[b]]), x[[b]], mu = mu)
    return(ret)
}
@}

In a similar way, we evaluate the gradients for each block and sum-up the
contributions by the shift parameters whereas the gradients for the
intercept parameters are only concatenated.

@d stratified negative score
@{
.snsc <- function(parm, x, mu = 0) {
    @<stratum prep@>
    ret <- numeric(length(bidx))
    for (b in seq_len(B)) {
        nsc <- .nsc(c(beta, intercepts[[b]]), x[[b]], mu = mu)
        ret[bidx] <- ret[bidx] + nsc[bidx]
        ret <- c(ret, nsc[-bidx])
    }
    return(ret)
}
@}

The score residuum is zero for an observation with weight zero, that is, a
row of zeros in the table.

@d stratified negative score residual
@{
.snsr <- function(parm, x, mu = 0) {
    @<stratum prep@>
    ret <- c()
    for (b in seq_len(B)) {
        idx <- attr(x[[b]], "idx")
        sr <- numeric(length(idx))
        sr[idx] <- .nsr(c(beta, intercepts[[b]]), x[[b]], mu = mu)
        ret <- c(ret, sr)
    }
    return(ret)
}
@}

<<glm-stratum>>=
(x <- as.table(array(c(10, 5, 7, 11, 8, 9,
                        9, 4, 8, 15, 5, 4), dim = c(2, 3, 2))))
d <- expand.grid(y = relevel(gl(2, 1), "2"), t = gl(3, 1), s = gl(2, 1))
d$x <- c(x)
m <- glm(y ~ s + t, data = d, weights = x, family = binomial())
logLik(m)
(cf <- coef(m))
vcov(m)[-(1:2),-(1:2)]
@@

<<glm-op-stratum>>=
(op <- optim(par = c("mt2" = 0, "mt3" = 0, "(Intercept 1)" = 0, "(Intercept 2)" = 0), 
             fn = .snll, gr = .snsc, 
             x = x, 
             method = "BFGS", 
             hessian = TRUE))
c(cf[-(1:2)] * -1, cf[1:2]) - op$par
logLik(m) + op$value
.snsr(op$par, x)
@@

@d stratified Hessian
@{
.shes <- function(parm, x, mu = 0) {
    @<stratum prep@>
    ret <- matrix(0, nrow = length(bidx), ncol = length(bidx))
    for (b in seq_len(B))
        ret <- ret + .hes(c(beta, intercepts[[b]]), x[[b]], mu = mu)
    ret
}
@}

<<glm-H-stratum>>=
H <- .shes(op$par, x)
solve(H)
vcov(m)[-(1:2),-(1:2)]
solve(op$hessian)[1:2,1:2]
@@

	
\chapter{Schur Complement}
\label{ch:schur}

@o Schur_src.c -cc
@{
#define STRICT_R_HEADERS
#include <R.h>
#include <Rinternals.h>
@<C_symtrisolve@>
@<NROW@>
@<NCOL@>
@<wcrossprod@>
@}

For a symmetric tridiagonal quadratic $N \times N$ matrix $\mA$ we compute $\X^\top \mA^{-1} \X$
utilising that the inverse $\mA^{-1}_{ij} = u_i v_j$ for $1 \le i \le j \le N$
can be characterised by two vectors $\uvec$ and $\vvec$, each of length $N$
\citep{Meurant1992}.

We begin with the diagonal $(a_1, \dots, a_N)^\top = \text{diag}(\mA)$ and the
negative lower- and upper off-diagonal $(-b_1, \dots, -b_{N - 1}) =
\text{diag}(\mA_{-N,-1})$. \cite{Meurant1992} starts with a decomposition
$\mA = \mU \D_U^{-1} \mU^\top$ with $(d_1, \dots, d_N)^\top =
\text{diag}(\D_U)$ where $\mU$ is upper triangular. The
decomposition also allows to compute the determinant of $\mA$ as $\prod_{i =
1}^N d_i$.

@d d vec
@{
d[n] = a[n];
det = d[n];
for (i = n - 1; i >= 0; i--) {
    d[i] = a[i] - pow(b[i], 2) / d[i + 1];
    /* DOI:10.1137/0613045 page 710: T = U D^-1 U^t with
       diag(U) = d (upper triangular),
       diag(D) = d (diagonal) => det(T) = prod(d)
    */
    det *= d[i];
}
@}

Following Proposition 1 in \cite{Meurant1992}, we compute $\uvec$

@d u vec
@{
u[0] = 1 / d[0];
prodb = 1.0;
for (i = 1; i <= n; i++) {
    prodb *= -b[i - 1] / d[i - 1];
    u[i] = prodb / d[i];
}
@}

Based on the next decomposition $\mA = \mL \D_L^{-1} \mL^\top$ with lower
triangular $\mL$ and $(\delta_1, \dots, \delta_N)^\top = \text{diag}(\D_L)$, we
compute

@d delta vec
@{
delta[0] = a[0];
for (i = 1; i <= n; i++)
    delta[i] = a[i] - pow(b[i - 1], 2) / delta[i - 1];
@}

and then, following Proposition 2, $\vvec$

@d v vec
@{
v[n] = 1 / (ans[n] * delta[n]);
v[0] = 1.0;
prodb = 1.0;
for (i = 1; i < n; i++) {
    prodb *= -b[n - i] / delta[n - i];
    v[n - i] = prodb * v[n];
}
@}

We wrap everything up in a function with arguments \code{a} and \code{b} of
length \code{n + 1} $=N$ and \code{n} $=N - 1$, respectively. The two
vectors $\uvec$ and $\vvec$ are stored in an $N \times 2$ real matrix \code{ans}.
We check if the determinant is larger than a small tolerance \code{tol}
before computing $\uvec$ and $\vvec$. The memory allocated for the $d$'s is
reused for computing $\delta$'s.

@d C_symtrisolve
@{
void C_symtrisolve (double *a, double *b, R_xlen_t n, double tol, double *ans)
{
    SEXP Rd;
    double *d, *delta, *u, *v, prodb, det;
    R_xlen_t i;

    /* output vectors */
    u = ans;
    v = ans + n + 1;

    /* n = N - 1 */
    PROTECT(Rd = allocVector(REALSXP, n + 1));
    d = REAL(Rd);

    @<d vec@>

    if (fabs(det) < tol) {
        error("Matrix not invertible");
    } else {
        @<u vec@>
        delta = d;
        @<delta vec@>
        @<v vec@>
    }
    UNPROTECT(1);
}
@}

In the next step, we compute the weighted crossproduct 
$\X^\top \mA^{-1} \X$$\X^\top \mA^{-1} \X$
without memory allocation for the full $N \times N$ matrix $\mA^{-1}$.
Because the resulting matrix is symmetric, we first compute the lower
triangular elements only.

@d lower wcrossprod
@{
for (p = 0; p < P; p++) {
    i = 0;
    dcs[i] = dx[p * N + i] * dvu[i];
    dcs[N + i] = dx[p * N + i] * dvu[N + i];
    for (i = 1; i < N; i++) {
        dcs[i] = dcs[i - 1] + dx[p * N + i] * dvu[i];
        dcs[N + i] = dcs[N + i - 1] + dx[p * N + i] * dvu[N + i];
    }
    for (j = 0; j < N; j++) {
        dxA = 0.0;
        dxA1 = dcs[N + j];
        dxA2 = dcs[N - 1] - dcs[j];
        dxA = dxA1 * dvu[j] + dxA2 * dvu[N + j];
        for (pp = p; pp < P; pp++)
            dans[p * P + pp] += dxA * dx[pp * N + j];
    }
}
@}

@d upper wcrossprod
@{
for (p = 0; p < P; p++) {
    for (pp = p + 1; pp < P; pp++)
        dans[pp * P + p] = dans[p * P + pp];
}
@}

The \proglang{R} interface requires access to the number of rows and columns
of matrices

@d NROW
@{
R_xlen_t NROW
(
    SEXP x
) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(XLENGTH(x));
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[0]);
    return((R_xlen_t) INTEGER(a)[0]);
}
@}

@d NCOL
@{
R_xlen_t NCOL
(
    SEXP x
) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[1]);
    return((R_xlen_t) INTEGER(a)[1]);
}
@}

@d wcrossprod
@{
SEXP R_wcrossprod (SEXP a, SEXP b, SEXP X, SEXP tol)
{

    SEXP ans, vu, cumsumvux;
    double *dans, *dx, dxA, dxA1, dxA2, *dvu, *dcs;
    R_xlen_t N, i, j;
    int p, pp, P;
    
    N = XLENGTH(a);

    if (NROW(X) != N)
        error("incorrect number of rows in X");
    if (!isReal(X))
        error("incorrect type of X");
    dx = REAL(X);
    P = (int) NCOL(X);

    PROTECT(ans = allocMatrix(REALSXP, P, P));
    dans = REAL(ans);

    if (XLENGTH(b) != N - 1)
        error("incorrect length of b");

    if (!isReal(a))
        error("incorrect type of a");

    if (!isReal(b))
        error("incorrect type of b");

    if (!isReal(tol))
        error("incorrect type of tol");

    PROTECT(vu = allocMatrix(REALSXP, N, 2));
    dvu = REAL(vu);
    C_symtrisolve(REAL(a), REAL(b), N - 1, REAL(tol)[0], dvu);
    PROTECT(cumsumvux = allocMatrix(REALSXP, N, 2));
    dcs = REAL(cumsumvux);

    for (p = 0; p < P * P; p++)
        dans[p] = 0.0;

    @<lower wcrossprod@>
    @<upper wcrossprod@>

    UNPROTECT(3);
    return(ans);
}
@}

@d R wcrossprod
@{
wcrossprod <- function(x, A, tol = .Machine$double.eps) {
    storage.mode(x) <- "double"
    .Call("R_wcrossprod", a = as.double(A$a), 
                          b = as.double(A$b),
                          X = x,
                          tol = tol)
}
@}

\chapter{Link Functions}
\label{ch:link}

@o linkfun.R -cp
@{
@<linkfun@>
@<logit@>
@<probit@>
@<cloglog@>
@<loglog@>
@}

@d linkfun
@{
.p <- function(link, q, ...)
    link$linkinv(q = q, ...)

.q <- function(link, p, ...)
    link$link(p = p, ...)

.d <- function(link, x, ...)
    link$dlinkinv(x = x, ...)

.dd <- function(link, x, ...)
    link$ddlinkinv(x = x, ...)

.ddd <- function(link, x, ...)
    link$dddlinkinv(x = x, ...)

.dd2d <- function(link, x, ...)
    link$dd2dlinkinv(x = x, ...)

linkfun <- function(alias, 
                    model, 
                    parm, 
                    link, 
                    linkinv,
                    dlinkinv, 
                    ddlinkinv,
                    ...) {

    ret <- list(alias = alias,
                model = model,
                parm = parm,
                link = link,
                linkinv = linkinv,
                dlinkinv = dlinkinv,
                ddlinkinv = ddlinkinv)
    if (is.null(ret$dd2d)) 
        ret$dd2d <- function(x) 
            ret$ddlinkinv(x) / ret$dlinkinv(x)
    ret <- c(ret, list(...))
    class(ret) <- "linkfun"
    ret
}
@}

@d logit
@{
logit <- function()
    linkfun(alias = "Wilcoxon",
            model = "proportional odds", 
            parm = "log-odds ratio",
            link = qlogis,
            linkinv = plogis,
            dlinkinv = dlogis,
            ddlinkinv = function(x) {
                p <- plogis(x)
                p * (1 - p)^2 - p^2 * (1 - p)
            },
            dddlinkinv = function(x) {
                ex <- exp(x)
                ifelse(is.finite(x), (ex - 4 * ex^2 + ex^3) / (1 + ex)^4, 0.0)
            },
            dd2d = function(x) {
                ex <- exp(x)
                (1 - ex) / (1 + ex)
            },
            parm2PI = function(x) {
               OR <- exp(x)
               ret <- OR * (OR - 1 - x)/(OR - 1)^2
               ret[abs(x) < .Machine$double.eps] <- 0.5
               return(ret)
            },
            PI2parm = function(p) {
               f <- function(x, PI)
                   x + (exp(-x) * (PI + exp(2 * x) * (PI - 1) + exp(x)* (1 - 2 * PI)))
               ret <- sapply(p, function(p) 
                   uniroot(f, PI = p, interval = 50 * c(-1, 1))$root)
               return(ret)
            },
            parm2OVL = function(x) 2 * plogis(-abs(x / 2))
    )
@}

@d loglog
@{
loglog <- function()
    linkfun(alias = "Lehmann", 
            model = "Lehmann", 
            parm = "log-reverse time hazard ratio",
            link = function(p, log.p = FALSE) {
                if (!log.p) p <- log(p)
                -log(-p)
            },
            linkinv = function(q, lower.tail = TRUE, log.p = FALSE) {
                ### p = exp(-exp(-q))
                if (log.p) {
                    if (lower.tail)
                        return(-exp(-q))
                    return(log1p(-exp(-exp(-q))))
                }
                if (lower.tail)
                    return(exp(-exp(-q)))
                -expm1(-exp(-q))
            },
            dlinkinv = function(x) 
                ifelse(is.finite(x), exp(- x - exp(-x)), 0.0),
            ddlinkinv = function(x) {
               ex <- exp(-x)
               ifelse(is.finite(x), exp(-ex - x) * (ex - 1.0), 0.0)
            },
            dddlinkinv = function(x) {
               ex <- exp(-x)
               ifelse(is.finite(x), exp(-x - ex) * (ex - 1)^2 - exp(-ex - 2 * x), 0.0)
            },
            dd2d = function(x) 
                expm1(-x),
            parm2PI = plogis,
            PI2parm = qlogis,
            parm2OVL = function(x) {
                x <- abs(x)
                rt <- exp(-x / (exp(x) - 1))
                ret <- rt^exp(x) + 1 - rt
                ret[abs(x) < .Machine$double.eps] <- 1
                x[] <- ret
                return(x)
            }
    )
@}

@d cloglog
@{
cloglog <- function()
    linkfun(alias = "Savage",
            model = "proportional hazards", 
            parm = "log-hazard ratio",
            link = function(p, log.p = FALSE) {
                if (log.p) p <- exp(p)
                log(-log1p(- p))
            },
            linkinv = function(q, lower.tail = TRUE, log.p = FALSE) {
                ### p = 1 - exp(-exp(q))
                ret <- exp(-exp(q))
                if (log.p) {
                    if (lower.tail)
                        return(log1p(-ret))
                    return(-exp(q))
                }
                if (lower.tail)
                    return(-expm1(-exp(q)))
                return(ret)
            },
            dlinkinv = function(x) 
                ifelse(is.finite(x), exp(x - exp(x)), 0.0),
            ddlinkinv = function(x) {
                ex <- exp(x)
                ifelse(is.finite(x), (ex - ex^2) / exp(ex), 0.0)
            },
            dddlinkinv = function(x) {
                ex <- exp(x)
                ifelse(is.finite(x), (ex - 3*ex^2 + ex^3) / exp(ex), 0.0)
            },
            dd2d = function(x)
               -expm1(x),
            parm2PI = plogis,
            PI2parm = qlogis,
            parm2OVL = function(x) {
                x <- abs(x)
                ret <- exp(x / (exp(-x) - 1)) - exp(-x / (exp(x) - 1)) + 1 
                ret[abs(x) < .Machine$double.eps] <- 1
                x[] <- ret
                return(x)
            }
    )
@}

@d probit
@{
probit <- function()
    linkfun(alias = "van der Waerden",
            model = "latent normal shift", 
            parm = "generalised Cohen's d",
            link = qnorm,
            linkinv = pnorm,
            dlinkinv = dnorm,
            ddlinkinv = function(x) 
                ifelse(is.finite(x), -dnorm(x = x) * x, 0.0), 
            dddlinkinv = function(x) 
                ifelse(is.finite(x), dnorm(x = x) * (x^2 - 1), 0.0),
            dd2d = function(x) -x,
            parm2PI = function(x) pnorm(x, sd = sqrt(2)),
            PI2parm = function(p) qnorm(p, sd = sqrt(2)),
            parm2OVL = function(x) 2 * pnorm(-abs(x / 2))
    )
@}



\chapter{ML Estimation}
\label{ch:ML}

@o ML.R -cp
@{
@<ML estimation@>
@<Strasser Weber@>
@<free1way@>
@<free1way methods@>
@<free1way summary@>
@<free1way confint@>
@<free1way interfaces@>
@}

@d setup
@{
K <- dim(x)[2L]
B <- dim(x)[3L]
xlist <- vector(mode = "list", length = B)
if (NS <- is.null(start))
    start <- rep.int(0, K - 1)
lwr <- rep(-Inf, times = K - 1)
for (b in seq_len(B)) {
    xb <- matrix(x[,,b, drop = TRUE], ncol = K)
    xw <- rowSums(abs(xb)) > tol
    xlist[[b]] <- xb[xw,,drop = FALSE]
    attr(xlist[[b]], "idx") <- xw
    lwr <- c(lwr, -Inf, rep.int(tol, times = sum(xw) - 2L))
    if (NS) {
        ecdf0 <- cumsum(rowSums(xlist[[b]]))
        ecdf0 <- ecdf0[-length(ecdf0)] / ecdf0[length(ecdf0)]
        Qecdf <- Q(ecdf0)
        start <- c(start, Qecdf[1], diff(Qecdf))
    }
}
upr <- rep(Inf, times = length(lwr))
@}

The profile negative log-likelihood can be evaluated for some of the
parameters in $\thetavec$ (denoted as \code{fix}), the remaining parameters
are updated. Note that \code{start} must contain the full parameter vector
$\thetavec$.

@d profile
@{
.profile <- function(start, fix = seq_len(K - 1)) {
    stopifnot(all(fix %in% seq_len(K - 1)))
    beta <- start[fix]
    ret <- optim(par = start[-fix], fn = function(par) {
                     p <- numeric(length(par) + length(fix))
                     p[fix] <- beta
                     p[-fix] <- par
                     .snll(p, x = xlist, mu = mu)
                 },
                 gr = function(par) {
                     p <- numeric(length(par) + length(fix))
                     p[fix] <- beta
                     p[-fix] <- par
                     .snsc(p, x = xlist, mu = mu)[-fix]
                 },
                 lower = lwr[-fix], upper = upr[-fix], method = "L-BFGS-B", 
                 hessian = FALSE, ...)
    p <- numeric(length(start))
    p[fix] <- beta
    p[-fix] <- ret$par
    ret$par <- p
    ret
}
@}

@d ML estimation
@{
.free1wayML <- function(x, link, mu = 0, start = NULL, fix = NULL, 
                   residuals = TRUE, score = TRUE, hessian = TRUE, 
                   tol = sqrt(.Machine$double.eps), ...) {

    stopifnot(is.table(x))
    dx <- dim(x)
    if (length(dx) == 2L)
        x <- as.table(array(c(x), dim = dx <- c(dx, 1L)))
    stopifnot(length(dx) == 3L)
    stopifnot(dim(x)[1L] > 1L)
    K <- dim(x)[2L]
    stopifnot(K > 1L)

    F <- function(q) .p(link, q = q)
    Q <- function(p) .q(link, p = p)
    f <- function(q) .d(link, x = q)
    fp <- function(q) .dd(link, x = q)

    @<setup@>

    @<cumsumrev@>
    @<negative logLik@>
    @<negative score@>
    @<negative score residuals@>
    @<Hessian@>
    @<stratified negative logLik@>
    @<stratified negative score@>
    @<stratified Hessian@>
    @<stratified negative score residual@>
    @<profile@>
 
    if (!length(fix)) {
        ret <- optim(par = start, 
                     fn = function(parm)
                         .snll(parm, x = xlist, mu = mu),
                     gr = function(parm)
                         .snsc(parm, x = xlist, mu = mu),
                     lower = lwr, upper = upr, method = "L-BFGS-B", 
                     hessian = FALSE, ...)
    } else if (length(fix) == length(start)) {
        ret <- list(par = start, 
                    value = .snll(start, x = xlist, mu = mu))
    } else {
        ret <- .profile(start, fix = fix)
    }

    if (is.null(fix) || (length(fix) == length(start)))
        parm <- seq_len(K - 1)
    else 
        parm <- fix
    if (any(parm >= K)) return(ret)

    ret$coefficients <- ret$par[parm]
    dn2 <- dimnames(x)[2L]
    names(ret$coefficients) <- cnames <- paste0(names(dn2), dn2[[1L]][1L + parm])

    if (score)
        ret$negscore <- .snsc(ret$par, x = xlist, mu = mu)[parm]
    if (hessian) {
        ret$hessian <- .shes(ret$par, x = xlist, mu = mu)
        if (length(parm) != nrow(ret$hessian))
            ret$hessian <- solve(ret$vcov <- solve(ret$hessian)[parm,parm])
        ret$vcov <- solve(ret$hessian)
        rownames(ret$vcov) <- colnames(ret$vcov) <- rownames(ret$hessian) <-
            colnames(ret$hessian) <-  cnames
    }
    if (residuals)
        ret$residuals <- .snsr(ret$par, x = xlist, mu = mu)

    ret$profile <- function(start, fix)
        .free1wayML(x, link = link, mu = mu, start = start, fix = fix, tol = tol, ...) 
    ret$table <- x

    class(ret) <- "free1wayML"
    ret
}
print.free1wayML <- function(x, ...)
    x$coefficients
@}

<<workhorse>>=

.free1wayML(x, logit())
op
N <- 10
a <- matrix(c(5, 6, 4,
                    3, 5, 7,
                    3, 4, 5,
                    3, 5, 6,
                    0, 0, 0,
                    4, 6, 5), ncol = 3, byrow = TRUE)
x <- as.table(array(c(a[1:3,], a[-(1:3),]), dim = c(3, 3, 2)))
x
(ret <- .free1wayML(x, logit()))
.free1wayML(x, logit(), start = ret$par, fix = 1:2)
.free1wayML(x, logit(), start = ret$par, fix = 2)
(ret <- .free1wayML(x, logit(), start = ret$par, fix = seq_along(ret$par)))
@@

\chapter{ML Inference}
\label{ch:MLinf}

\section{Wald}

@d statistics
@{
if (test == "Wald") {
    if (alternative == "two.sided") {
        STATISTIC <- c("Wald X-squared" = c(crossprod(cf, x$hessian %*% cf)))
        DF <- c("df" = length(parm))
        PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
    } else {
        STATISTIC <- c("Wald Z" = c(cf * sqrt(x$hessian)))
        PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
    }
} else if (test == "LRT") {
    par <- x$par
    par[parm] <- value
    unll <- x$value ### neg logLik
    rnll <- x$profile(par, parm)$value ### neg logLik
    STATISTIC <- c("logLR X-squared" = - 2 * (unll - rnll))
    DF <- c("df" = length(parm))
    PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
} else if (test == "Rao") {
    par <- x$par
    par[parm] <- value
    ret <- x$profile(par, parm)
    if (alternative == "two.sided") {
        STATISTIC <- c("Rao X-squared" = c(crossprod(ret$negscore, ret$vcov %*% ret$negscore)))
        DF <- c("df" = length(parm))
        PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
    } else {
        STATISTIC <- c("Rao Z" = -ret$negscore * sqrt(ret$vcov))
        PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
    }
} else if (test == "Permutation") {
    par <- x$par
    par[parm] <- value
    ret <- x$profile(par, parm)
    sc <- -ret$negscore
    if (length(cf) == 1L)
       sc <- sc / sqrt(ret$hessian)
    Esc <- sc - x$perm$Expectation
    if (alternative == "two.sided" && length(cf) > 1L) {
        STATISTIC <- c("Perm X-squared" = sum(Esc %*% solve(x$perm$Covariance) * Esc))
        ps <- x$perm$permStat
        if (!is.null(x$perm$permStat))
            PVAL <- mean(ps > STATISTIC + tol)
        else {
            DF <- c("df" = x$perm$DF)
            PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
        }
    } else {
        STATISTIC <- c("Perm Z" = Esc / sqrt(x$perm$Covariance))
        if (!is.null(x$perm$permStat)) {
            if (alternative == "two.sided")
                PVAL <- mean(abs(x$perm$permStat) < abs(STATISTIC) - tol)
            else if (alternative == "less")
                PVAL <- mean(x$perm$permStat < STATISTIC - tol)
            else
                PVAL <- mean(x$perm$permStat > STATISTIC + tol)
        } else {
            if (alternative == "two.sided")
                PVAL <- pchisq(STATISTIC^2, df = 1, lower.tail = FALSE)
            else
                PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
        }
    }
}
@}

\chapter{Permutation Inference}

@d Strasser Weber
@{
.SW <- function(res, xt) {

    if (length(dim(xt)) == 3L) {
        res <- matrix(res, nrow = dim(xt)[1L], ncol = dim(xt)[3])
        STAT <-  Exp <- Cov <- 0
        for (b in seq_len(dim(xt)[3L])) {
            sw <- .SW(res[,b, drop = TRUE], xt[,,b, drop = TRUE])
            STAT <- STAT + sw$Statistic
            Exp <- Exp + sw$Expectation
            Cov <- Cov + sw$Covariance
        }
        return(list(Statistic = STAT, Expectation = as.vector(Exp),
                    Covariance = Cov))
    }

    Y <- matrix(res, ncol = 1, nrow = length(xt))
    weights <- c(xt)
    x <- gl(ncol(xt), nrow(xt))
    X <- model.matrix(~ x, data = data.frame(x = x))[,-1L,drop = FALSE]

    w. <- sum(weights)
    wX <- weights * X
    wY <- weights * Y
    ExpX <- colSums(wX)
    ExpY <- colSums(wY) / w.
    CovX <- crossprod(X, wX)
    Yc <- t(t(Y) - ExpY)
    CovY <- crossprod(Yc, weights * Yc) / w.
    Exp <- kronecker(ExpY, ExpX)
    Cov <- w. / (w. - 1) * kronecker(CovY, CovX) -
           1 / (w. - 1) * kronecker(CovY, tcrossprod(ExpX))
    STAT <- crossprod(X, wY)
    list(Statistic = STAT, Expectation = as.vector(Exp),
         Covariance = Cov)
}

.perm <- function(res, xt, B = 10000) {

    if (length(dim(xt)) == 2L)
        xt <- as.table(array(xt, dim = c(dim(xt), 1)))

    res <- matrix(res, nrow = dim(xt)[1L], ncol = dim(xt)[3L])
    stat <- 0
    ret <- .SW(res, xt)
    if (dim(xt)[2L] == 2L) {
        ret$testStat <- c((ret$Statistic - ret$Expectation) / sqrt(ret$Covariance))
    } else {
        ES <- t(ret$Statistic - ret$Expectation)
        ret$testStat <- sum(ES %*% solve(ret$Covariance) * ES)
    }
    ret$DF <- dim(xt)[2L] - 1L

    if (B) {
        for (j in 1:dim(xt)[3L]) {
           rt <- r2dtable(B, r = rowSums(xt[,,j]), c = colSums(xt[,,j]))
           stat <- stat + sapply(rt, function(x) colSums(x[,-1L, drop = FALSE] * res[,j]))
        }
        if (dim(xt)[2L] == 2L) {
             ret$permStat <- (stat - ret$Expectation) / sqrt(ret$Covariance)
        } else {
            ES <- t(matrix(stat, ncol = B) - ret$Expectation)
            ret$permStat <- rowSums(ES %*% solve(ret$Covariance) * ES)
        }
    }
    ret
}
@}

<<SW>>=
w <- gl(2, 15)
(s <- .SW(r <- rank(u <- runif(length(w))), model.matrix(~ 0 + w)))
ps <- .perm(r, model.matrix(~ 0 + w), B = 100000)
ps$testStat
mean(abs(ps$permStat) > abs(ps$testStat) - .Machine$double.eps)
pchisq(ps$testStat^ifelse(ps$DF == 1, 2, 1), df = ps$DF, lower.tail = FALSE)
kruskal.test(u ~ w)
library("coin")
kruskal_test(u ~ w, distribution = approximate(100000))
@@


\chapter{Test for Equal Distributions in a Stratified One-way Layout}

Distribution-free

@d free1way
@{
free1way.test <- function(x, ...)
    UseMethod("free1way.test")
free1way.test.table <- function(object, link = c("logit", "probit", "cloglog", "loglog"), mu = 0, B = 0, ...)
{
    d <- dim(object)
    dn <- dimnames(object)
    DNAME <- NULL
    if (!is.null(dn)) {
        DNAME <- paste(names(dn)[1], "by", names(dn)[2])
        if (length(dn) == 3L)
            DNAME <- paste(DNAME, "with strata", names(dn)[3])
    }

    if (!inherits(link, "linkfun")) {
        link <- match.arg(link)
        link <- do.call(link, list())
    }

    ret <- .free1wayML(object, link = link, mu = mu, ...)
    ret$link <- link
    ret$DNAME <- DNAME
    ret$METHOD <- paste(link$alias, "test against", link$model, "alternatives")

    cf <- ret$par
    cf[idx <- seq_len(d[2L] - 1L)] <- 0
    pr <- ret$profile(cf, idx)
    if (d[2L] == 2L)
        res <- pr$residuals / sqrt(pr$hessian)
    else
        res <- pr$residuals
    ret$perm <- .perm(res, object, B = B)

    class(ret) <- "free1way.test"
    return(ret)
}
@}

@d free1way methods
@{
coef.free1way.test <- function(object, ...)
    object$coefficients
vcov.free1way.test <- function(object, ...)
    object$vcov
logLik.free1way.test <- function(object, ...)
    -object$value
@}

@d free1way summary
@{
print.free1way.test <- function(x, test = c("Permutation", "Wald", "LRT", "Rao"), 
                                alternative = c("two.sided", "less", "greater"), 
                                tol = .Machine$double.eps, ...)
{

    test <- match.arg(test)
    alternative <- match.arg(alternative)

    ### global
    cf <- coef(x)
    if ((length(cf) > 1L || test == "LRT") && alternative != "two.sided") 
        stop("Cannot compute one-sided p-values")

    DF <- NULL
    parm <- seq_along(cf)
    value <- 0

    @<statistics@>

    RVAL <- list(statistic = STATISTIC, parameter = DF, p.value = PVAL, 
        null.value = ret$mu, alternative = alternative, method = x$METHOD, 
        data.name = x$DNAME, estimate = cf)
    class(RVAL) <- "htest"
    RVAL

}
summary.free1way.test <- function(object, alternative = c("two.sided", "less", "greater"), ...)
{

    alternative <- match.arg(alternative)

    ESTIMATE <- coef(object)
    SE <- sqrt(diag(vcov(object)))
    STATISTIC <- unname(ESTIMATE / SE)
    if (alternative == "less") {
        PVAL <- pnorm(STATISTIC)
    } else if (alternative == "greater") {
        PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
    } else {
        PVAL <- 2 * pnorm(-abs(STATISTIC))
    }
    cfmat <- cbind(ESTIMATE, SE, STATISTIC, PVAL)
    colnames(cfmat) <- c(object$link$parm, "Std. Error", "z value",
                         switch(alternative, "two.sided" = "P(>|z|)",
                                             "less" = "P(<z)",
                                             "greater" = "P(>z)"))
    ret <- list(coefficients = cfmat)
    class(ret) <- "summary.free1way.test"
    return(ret)
}
print.summary.free1way.test <- function(x, ...) {
    cat("Coefficients:\n")
    printCoefmat(x$coefficients)
}
@}

@d free1way confint
@{
confint.free1way.test <- function(object, parm,
    level = .95, test = c("Permutation", "Wald", "LRT", "Rao"), ...)
{

    test <- match.arg(test)
    conf.level <- 1 - (1 - level) / 2

    cf <- coef(object)
    if (missing(parm)) 
        parm <- seq_along(cf)

    ESTIMATE <- cf[parm]
    qSE <- qnorm(conf.level) * sqrt(diag(vcov(object)))[parm]
    CINT <- cbind(ESTIMATE - qSE, ESTIMATE + qSE)
    colnames(CINT) <- paste(100 * c(1 - conf.level, conf.level), "%")
    if (test == "Wald")
        return(CINT)

    CINT[] <- cbind(ESTIMATE - qSE * 5, ESTIMATE + qSE * 5)

    sfun <- function(value, parm, quantile) {
        x <- object
        alternative <- "two.sided"
        tol <- .Machine$double.eps
        @<statistics@>
        return(STATISTIC - quantile)
    }

    if (test == "Permutation") {
        stopifnot(length(cf) == 1L)
        if (is.null(object$perm$permStat)) {
            qu <- qnorm(conf.level) * c(-1, 1)
        } else {
            qu <- quantile(object$perm$permStat, probs = c(1 - conf.level, conf.level))
            att.level <- mean(object$perm$permStat > qu[1] & object$perm$permStat < qu[2])
        }
    } else {
        qu <- rep.int(qchisq(conf.level, df = 1), 2)
    }

    for (p in parm) {
        CINT[p, 1] <- uniroot(sfun, interval = c(CINT[p,1], cf[p]), parm = p, quantile = qu[2])$root
        CINT[p, 2] <- uniroot(sfun, interval = c(cf[p], CINT[p, 2]), parm = p, quantile = qu[1])$root
    }
    if (test == "Permutation") attr(CINT, "Attained level") <- att.level
    return(CINT)
}
# power.free1way.test <- function()
@}

<<free>>=
ft <- free1way.test(x)
coef(ft)
vcov(ft)
summary(ft)
library("multcomp")
summary(glht(ft), test = Chisqtest())
summary(ft, test = "Rao")
summary(ft, test = "permRao")
summary(ft, test = "LRT")
confint(glht(ft), calpha = univariate_calpha())
confint(ft, test = "Rao")
confint(ft, test = "Wald")
confint(ft, test = "LRT")
@@

@d free1way interfaces
@{
free1way.test.formula <- function(formula, data, weights, subset, na.action = na.pass, ...)
{
    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    group <- 2
    if (length(attr(terms(formula[-2L]), "term.labels")) > 2L)
        stop("'formula' missing or incorrect")
    strata <- NA
    if (length(attr(terms(formula[-2L]), "term.labels")) == 2L)
       strata <- 3
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ") # works in all cases
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    w <- as.vector(model.weights(mf))
    y <- mf[[response]]
    g <- factor(mf[[group]])
    if(nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
    if (is.na(strata)) {
        ## Call the default method.
        RVAL <- free1way.test(y = y, x = g, weights = w, ...)
    } else {
        st <- factor(mf[[strata]])
        RVAL <- free1way.test(y = y, x = g, z = st, weights = w, ...)
    }
    RVAL$data.name <- DNAME
    RVAL
}
 
free1way.test.numeric <- function(y, x, z = NULL, weights = NULL, nbins = 0, ...) {

    DNAME <- paste(deparse1(substitute(x)), "and",
                   deparse1(substitute(y)))
    if (!is.null(z))
        DNAME <- paste(DNAME, "stratified by", deparse1(substitute(z)))

    uy <- unique(y)
    if (nbins && nbins < length(uy)) {
        nbins <- ceiling(nbins)
        breaks <- c(-Inf, quantile(y, prob = seq_len(nbins) / (nbins + 1L)), Inf)
    } else {
        breaks <- c(-Inf, sort(uy), Inf)
    }
    r <- cut(y, breaks = breaks)[, drop = TRUE]
    RVAL <- free1way.test(y = r, x = x, z = z, weights = weights, ...)
    RVAL$data.name <- DNAME
    RVAL$method <- gsub("Ordinal|Binary", "Semiparametric", RVAL$method)
    RVAL
}

free1way.test.factor <- function(y, x, z = NULL, weights = NULL, ...) {

    DNAME <- paste(deparse1(substitute(x)), "and",
                   deparse1(substitute(y)))
    if (!is.null(z))
        DNAME <- paste(DNAME, "stratified by", deparse1(substitute(z)))

    stopifnot(is.factor(x))
    d <- data.frame(y = y, x = x, w = 1)
    if (!is.null(weights)) d$w <- weights
    if (!is.null(z)) {
        d$z <- z
        tab <- xtabs(w ~ y + x + z, data = d)
    } else {
        tab <- xtabs(w ~ y + x, data = d)
    }
    RVAL <- free1way.test(tab, ...)
    RVAL$data.name <- DNAME
    RVAL
}
@}

<<formula>>=
set.seed(29)
N <- 25
w <- gl(2, N)
y <- rlogis(length(w), location = c(0, 1)[w])
ft <- free1way.test(y ~ w, B = 10000)
summary(ft)
print(ft, test = "Permutation", alternative = "less")
print(ft, test = "Permutation", alternative = "greater")
print(ft, test = "Permutation")
print(ft, test = "Wald", alternative = "less")
print(ft, test = "Wald", alternative = "greater")
print(ft, test = "Wald")
print(ft, test = "LRT")
print(ft, test = "Rao", alternative = "less")
print(ft, test = "Rao", alternative = "greater")
print(ft, test = "Rao")
library("lehmann")
trafo.test(y ~ w)
summary(ft, test = "LRT")
summary(ft, test = "permRao")
summary(ft, test = "Rao")
confint(ft)
confint(ft, test = "LRT")
confint(ft, test = "Wald")
confint(ft, test = "Rao")
wilcox.test(y ~ w)
library("rms")
rev(coef(or <- orm(y ~ w)))[1]
logLik(or)
logLik(ft)
ci <- confint(or)
ci[nrow(ci),]
vcov(or)
vcov(ft)
@@

\chapter*{Index}

\section*{Files}

@f

\section*{Fragments}

@m

\section*{Identifiers}

@u

\bibliographystyle{plainnat}
\bibliography{refs}

\end{document}
