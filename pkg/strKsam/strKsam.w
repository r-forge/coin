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


%\usepackage{underscore}

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
\newcommand{\sX}{\mathcal{X}}
\newcommand{\sY}{\mathcal{Y}}
\newcommand{\T}{\mathbf{T}}
\newcommand{\x}{\mathbf{x}}
\renewcommand{\a}{\mathbf{a}}
\newcommand{\xn}{\mathbf{x}_{\text{new}}}
\newcommand{\y}{\mathbf{y}}
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

\title{Stratified $K$ Sample Inference}

\begin{document}

\pagenumbering{roman}
\maketitle
\tableofcontents

\chapter{Model and Parameterisation}
\pagenumbering{arabic}

Treatment group $\rT \in \{1, \dots, K\}, K \ge 2$, outcome $Y \in \samY$ at least ordered,
stratum $S \in \{1, \dots, B\}$ in $B \ge 1$ blocks with conditional
cumulative distribution function (cdf)
$F_Y(y \mid S = b, \rT = k) = \Prob(Y \le y \mid S = b, \rT = k)$.

\paragraph{Model}

With $g: [0,1] \times \R \rightarrow \R$, we model
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
$$F_Y(y \mid S = b, \rT = k) = F(F^{-1}(F_Y(y \mid S = b, \rT = 1)) - \beta_k).$$
The negative shift term ensures that positive values of $\beta_k$ correspond
to the situation of outcomes being stochastically larger in group $k$ than
in group one.

Note that $F(z) = \exp(-\exp(-z))$ gives rise to $g_\text{L}$, 
$F(z) = 1 - \exp(-\exp(z))$ to $g_\text{PH}$, $F = \text{expit}$ to
$g_\text{PO}$, and $F = \Phi$ to $g_\text{Cd}$.

\paragraph{Hypthesis}

We are interested in inference for $\beta_2, \dots, \beta_K$, in terms of
confidence intervals and hypothesis tests of the form

$$H_0: \beta_k = \mu_k, \quad k = 2, \dots, K.$$

\paragraph{Likelihood}

For an ordered categorical outcome $Y$ from sample space $\samY = \{y_1 < y_2 < \cdots <
y_C\}$, we parameterise the model in terms of intercept ($\vartheta_\cdot$) and
shift ($\beta_\cdot$) parameters

$$F_Y(y_c \mid S = b, \rT = k) = F(\vartheta_{c,b} - \beta_k), \quad c = 1, \dots, C.$$

The $C - 1$ intercept parameters are block-specific and monotone increasing
$\vartheta_{0,b} = -\infty < \vartheta_{1,b} < \cdots < \vartheta_{C,b} = \infty$.
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

For the $i$th observation $(y_i = y_c, s_i = b, \rt_i = k)$ from block $b$, the log-likelihood contribution is
$$\log(\Prob(y_{c - 1} < Y \le y_c \mid S = b, \rT = k)) = \log(F(\vartheta_{c,b} - \beta_k) - F(\vartheta_{c - 1,b} - \beta_k)).$$

For an absolutely continuous outcome $Y \in \R$, we define $y_c := y_{(c)}$,
the $c$th distinct ordered observation in the sample. The log-likelihood
above is then the empirical or nonparametric log-likelihood.

We represent the data in form of a $C \times K \times B$ contingency table,
whose element $(c, k, b)$ is the number of observations $(y = y_c, s = b,
\rt = k)$.
	
\chapter{Parameter Estimation}

@o strKsam_src.R -cp
@{
@<cumsumrev@>
@<negative logLik@>
@<negative score@>
@<Hessian@>
@<stratified negative logLik@>
@<stratified negative score@>
@<stratified Hessian@>
@}


We start implementing the log-likelihood function for parameters \code{parm}
$= \thetavec$ (assuming only a single block) with data from a two-way $C
\times K$ contingency table \code{x}. 

From $\thetavec$, we first extract the shift parameters $\beta_\cdot$ and
then the intercept parameters $\vartheta_\cdot$, compute the differences
$\vartheta_{c,1} - \beta_k$ and evaluate the probability
\code{prb} $ = \Prob(y_{c - 1} < Y \le y_c \mid S = 1, \rT = k)$:

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
.nll <- function(parm, x, mu = numeric(ncol(x) - 1L)) {
    @<parm to prob@>
    - sum(x * log(prb))
}
@}

It is important to note that, with $F$ corresponding to distribution with
log-concave density $f$, the negative log-likelihood is a convex function of
the parameters $\thetavec$, and thus we can solve the corresponding
constrained minimisation problem quickly and reliably.

To speed things up, we implement the gradient of the negative
log-likelihood, the negative score function for the parameters in
$\thetavec$. The score function for the empirical likelihood, evaluated at
parameters $\vartheta_\cdot$ and $\beta_\cdot$ is given in many places
\citep[for example in][Formula~(2)]{Hothorn_Moest_Buehlmann_2017}. The score
involves $f = F^\prime$:

@d negative score
@{
.nsc <- function(parm, x, mu = numeric(ncol(x) - 1L)) {
    @<parm to prob@>
    ftmb <- f(tmb)

    ret <- numeric(length(parm))
    zu <- x * ftmb[- 1, , drop = FALSE] / prb
    zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
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
    rev(cumsum(rev(z)))
@}

We also need access to the observed Fisher information of the shift
parameters. We proceed by implementing the Hessian for the intercept
($\vartheta_\cdot$) and shift ($\beta_\cdot$) parameters, as given in Formula~(4) of
\cite{Hothorn_Moest_Buehlmann_2017} first. This partitioned matrix
\begin{eqnarray*}
\mH(\vartheta_1, \dots, \vartheta_{C - 1}, \beta_2, \dots, \beta_K) = 
\left(\begin{array}{ll}
\mA & \X \\
\X^\top & \Z
\end{array} \right)
\end{eqnarray*}
consists of a tridiagonal $\mA \sim (C-1,C-1)$, a diagonal $\Z \sim (K - 1, K -
1)$, and a full $\X \sim (C - 1, K - 1)$ matrix. In a second step, we
compute the Fisher information matrix for the shift parameters only by means
of the Schur complement $\Z - \X^\top \mA^{-1} \X$.

In addition to probabilities \code{prb}, the Hessian necessitates the
computation of $f()$ and $f^\prime()$ for linear functions of the
parameters. We start preparing these objects


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

i1 <- length(intercepts) - 1
i2 <- 1
@}

The off-diagonal elements of $\mA$ are now available as
@d Aoffdiag
@{
Aoffdiag <- -rowSums(x * dpu * dpl / prb^2)[-i2]
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

Instead of computing the inverse of $\mH$ directly, we return a list
containing the components of this block matrix for later processing
elsewhere.

@d Hessian
@{
.hes <- function(parm, x, mu = numeric(ncol(x) - 1L)) {
    @<parm to prob@>

    @<Hessian prep@>

    @<Aoffdiag@>
    @<Adiag@>
    @<X and Z@>

    return(list(a = Adiag, b = Aoffdiag, X = X, Z = Z))
}
@}

We start with an example involving $K = 3$ groups for a binary outcome and
use a binary logistic regression model to estimate the two log-odds ratios
$\beta_2$ and $\beta_3$ along with their estimated covariance
<<glm>>=
source("strKsam_src.R")
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
@@

Parameter estimates and the in-sample log-likelihood are practically
identical. We now turn to the inverse Hessian of the shift terms, first
defining the derivative of the density of the standard logistic distribtion
<<glm-H>>=
fp <- function(x) {
    p <- plogis(x)
    p * (1 - p)^2 - p^2 * (1 - p)
}
H <- .hes(op$par, x)
solve(H$Z - crossprod(H$X,  H$X / H$a))
vcov(m)[-1,-1]
solve(op$hessian)[1:2,1:2]
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
    sidx <- factor(rep(seq_len(B), times = C), levels = seq_len(B))
}
bidx <- seq_len(K - 1L)
beta <- parm[bidx]
intercepts <- split(parm[-bidx], sidx)
if (is.null(mu)) mu <- numeric(K - 1L)
@}

@d stratified negative logLik
@{
.snll <- function(parm, x, mu = NULL) {
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
.snsc <- function(parm, x, mu = NULL) {
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
@@

@d stratified Hessian
@{
.shes <- function(parm, x, mu = NULL) {
    @<stratum prep@>
    ret <- matrix(0, nrow = length(bidx), ncol = length(bidx))
    for (b in seq_len(B)) {
        H <- .hes(c(beta, intercepts[[b]]), x[[b]], mu = mu)
        ret <- ret + do.call(lehmann:::Schur_symtri, H)
    }
    ret
}
@}

<<glm-H-stratum>>=
H <- .shes(op$par, x)
solve(H)
vcov(m)[-(1:2),-(1:2)]
solve(op$hessian)[1:2,1:2]
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
