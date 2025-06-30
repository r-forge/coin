\documentclass[a4paper]{report}

%\VignetteIndexEntry{Stratified K-sample Inference}
%\VignetteDepends{free1way,rms,coin,multcomp,survival}
%\VignetteKeywords{conditional inference, conditional Monte Carlo}}
%\VignettePackage{free1way}

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

\title{Semiparametrically Efficient Population and Permutation Inference in 
       Distribution-free Stratified $K$-sample Oneway Layouts}

\begin{document}

\pagenumbering{roman}
\maketitle
\tableofcontents

\chapter{Model and Parameterisation}
\label{ch:model}
\pagenumbering{arabic}

We consider $K$ treatment groups $\rT \in \{1, \dots, K\}, K \ge 2$ for an
at least ordered outcome $Y \in \samY$ observed in
stratum $S \in \{1, \dots, B\}$ out of $B \ge 1$ blocks with conditional
cumulative distribution function (cdf)
$F_Y(y \mid S = b, \rT = k) = \Prob(Y \le y \mid S = b, \rT = k)$. Detecting
and describing differential distributions arising from different treatments
is our main objetive. We refer to the first treatment $\rT = 1$ as
``control''.

\paragraph{Model}

With model function $g: [0,1] \times \R \rightarrow [0,1]$, we describe
the conditional distribution under treatment $k$ as a function of the 
conditional distribution under control and a scalar parameter
$\beta_k$:
\begin{eqnarray*}
F(y \mid S = b, \rT = k) = g(F(y \mid S = b, \rT = 1), \beta_k).
\end{eqnarray*}
The model is assumed to hold for all blocks $b = 1,
\dots, B$, treatments $k = 2, \dots, K$, and outcome values $y \in \samY$ based on parameters
$\beta_2, \dots, \beta_K \in \R$. For notational convenience, we define $\beta_1 := 0$. 

This model formulation gives rise to several specific models, for
example, $g_\text{L}(p, \beta) = p^{\exp(-\beta)}$ (Lehmann alternatives),
$g_\text{PH}(p, \beta) = 1 - (1 -
p)^{\exp(-\beta)}$ (proportional hazards),
$g_\text{PO}(p, \beta) = \text{expit}(\text{logit}(p) - \beta)$ (proportional
odds), or $g_\text{Cd}(p, \beta) =
\Phi(\Phi^{-1}(p) - \beta)$ (generalised Cohen's $d$).

Instead of directly working with $g$, we parameterise the model in terms of
some absolute continuous cdf $F$ with log-concave density $f = F^\prime$
and corresponding derivative $f^\prime$. The location model 
\begin{eqnarray} \label{model}
F_Y(y \mid S = b, \rT = k) = F\left(F^{-1}\left(F_Y(y \mid S = b, \rT = 1)\right) - \beta_k\right), \quad k = 2, \dots, K
\end{eqnarray}
describes different distributions by means of shift parameter on a latent
scale defined by $F$. The negative shift term ensures that positive values of $\beta_k$ correspond
to the situation of outcomes being stochastically larger in group $k$
compared to control.

The choise $F(z) = \exp(-\exp(-z))$ gives rise to $g_\text{L}$, 
$F(z) = 1 - \exp(-\exp(z))$ corresponds to $g_\text{PH}$, $F = \text{expit}$
leads to  $g_\text{PO}$, and $F = \Phi$ results in $g_\text{Cd}$. The choice
of $F$ is made a priori and determines the interpretation of $\beta_k$. 

This document describes the implementation of estimators of these shift parameters,
as well as of confidence intervals and formal hypothesis tests for contrasts thereof under
the permutation and population model.

\paragraph{Hypthesis}

We are interested in inference for $\beta_2, \dots, \beta_K$, in terms of
confidence intervals and hypothesis tests of the form
\begin{eqnarray*}
& & H_0: \beta_k - \mu_k = 0, \text{``two.sided''}, \quad k = 2, \dots, K, \\
& & H_0: \beta_k - \mu_k \ge 0, \text{``less''}, \quad k = K = 2, \\
& & H_0: \beta_k - \mu_k \le 0, \text{``greater''}, \quad k = K = 2,
\end{eqnarray*}
with the latter two options only for the two-sample case ($K = 2$).

\paragraph{Likelihood}

For an ordered categorical outcome $Y$ from sample space $\samY = \{y_1 < y_2 < \cdots <
y_C\}$, we parameterise the model in terms of intercept ($\vartheta_\cdot$) and
shift ($\beta_\cdot$) parameters
\begin{eqnarray*}
F_Y(y_c \mid S = b, \rT = k) = F(\vartheta_{c,b} - \beta_k), \quad c = 1, \dots,
C,
\end{eqnarray*}
that is we replace the transformed control outcome $F^{-1}\left(F_Y(y_c \mid S = b, \rT = 1)\right) =
\vartheta_{c,b}$ with a corresponding intercept parameter.
These $C - 1$ intercept parameters are block-specific and monotone increasing
$\vartheta_{0,b} = -\infty < \vartheta_{1,b} < \cdots < \vartheta_{C,b} = \infty$
within each block $b = 1, \dots, B$.

We collect all model parameters in a vector
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
featuring contrasts of the intercept parameters $\vartheta_{\cdot}$ such that monotonicity of
the intercept parameters can be ensured by box constraints for $\thetavec$.

For the $i$th observation $(y_i = y_c, s_i = b, \rt_i = k)$ from block $b$
under treatment $k$, the log-likelihood contribution is
\begin{eqnarray*}
\log(\Prob(y_{c - 1} < Y \le y_c \mid S = b, \rT = k)) = \log(F(\vartheta_{c,b} - \beta_k) - F(\vartheta_{c - 1,b} - \beta_k)).
\end{eqnarray*}

For an absolutely continuous outcome $Y \in \R$, we define $y_c := y_{(c)}$,
the $c$th distinct ordered observation in the sample. The log-likelihood
above is then the empirical or nonparametric log-likelihood.

If observations were independently right-censored, the contribution of the
event $Y > \tilde{y}$ to the log-likelihood is
\begin{eqnarray*}
\log(\Prob(Y > \tilde{y} \mid S = b, \rT = k)) = \log(1 - F(\vartheta_{c - 1,b} - \beta_k))
\end{eqnarray*}
where $y_{c - 1} = \max \{y \in \samY \mid y \le \tilde{y}\}$, that is,
observations right-censored between $y_{c - 1}$ and $y_c$ correspond to the
parameter $\vartheta_{c - 1,b}$.

Maximising this form of the log-likelihood leads to semiparametrically efficient
estimators \citep[Chapter 15.5][]{vdVaart1998}. In this framework, tests
against deviations from the hyptheses $H_0$ above are locally most
powerful rank tests, for example against proportional odds ($F =
\text{expit})$ or proportional hazards alternatives 
\citep[$F(z) = 1 - \exp(-\exp(z))$,][Example 15.16]{vdVaart1998}.

We represent the data in form of a $C \times K \times B$ contingency table,
whose element $(c, k, b)$ is the number of observations with configuration $(y = y_c, s = b,
\rt = k)$. In the presence of right-censoring, a fourth dimension is added 
($C \times K \times B \times 2)$ whose first $C \times K \times B$ table presents
right-censoring and the second table contains numbers of events.


	
\chapter{Parameter Estimation}
\label{ch:est}

<<localfun, echo = FALSE>>=
@<cumsumrev@>
@<table2list@>
@<negative logLik@>
@<negative score@>
@<negative score residuals@>
@<Hessian@>
@<stratified negative logLik@>
@<stratified negative score@>
@<stratified Hessian@>
@<stratified negative score residual@>
@<ML estimation@>
@<Strasser Weber@>
@<resampling@>
@<linkfun@>
@<logit@>
@<probit@>
@<cloglog@>
@<loglog@>
@@

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
if (rightcensored) {
    prb <- pmax(1 - Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
} else {
    prb <- pmax(Ftmb[- 1L, , drop = FALSE] - 
                Ftmb[- nrow(Ftmb), , drop = FALSE], sqrt(.Machine$double.eps))
} 
@}

If the table \code{x} represents right-censored observations, we compute
\code{prb} $ = 1 - \Prob(Y \le y_c \mid S = 1, \rT = k)$.

With default null values $\mu_k = 0, k = 2, \dots, K$, we define the
negative log-likelihood function as the weighted (by number of observations) sum of
the log-probabilities

@d negative logLik
@{
.nll <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<parm to prob@>
    return(- sum(x * log(prb)))
}
@}

The code assumes that all elements of the margins of the table \code{x} are larger than
zero; otherwise, the corresponding parameter is not identified. We will
handle such situation at a higher level later on.

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

@d density prob ratio
@{
ftmb <- f(tmb)
zu <- x * ftmb[- 1, , drop = FALSE] / prb
if (rightcensored) zu[] <- 0 ### derivative of a constant
zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
@}

and then compute the negative score function:

@d negative score
@{
.nsc <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<parm to prob@>
    @<density prob ratio@>

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
intercepts needs this small helper function:

@d cumsumrev
@{
.rcr <- function(z)
    # Reduce('+', z, accumulate = TRUE, right = TRUE)
    rev.default(cumsum(rev.default(z)))
@}

(<TH>maybe add \code{rev = TRUE} to \code{cumsum}?</TH>).

In addition, we define negative score residuals, that is, the derivative of the
negative log-likelihood with respect to an intercept term constrained to
zero:

@d negative score residuals
@{
.nsr <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<parm to prob@>
    @<density prob ratio@>

    ret <- rowSums(zu - zl) / rowSums(x)
    ret[!is.finite(ret)] <- 0
    ret
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
\beta_k)$. We start preparing these objects, keeping in mind to remove terms
not being present under right-censoring:

@d Hessian prep
@{
ftmb <- f(tmb)
fptmb <- fp(tmb)

dl <- ftmb[- nrow(ftmb), , drop = FALSE]
du <- ftmb[- 1, , drop = FALSE]
if (rightcensored) du[] <- 0
dpl <- fptmb[- nrow(ftmb), , drop = FALSE]
dpu <- fptmb[- 1, , drop = FALSE]
if (rightcensored) dpu[] <- 0
dlm1 <- dl[,-1L, drop = FALSE]
dum1 <- du[,-1L, drop = FALSE]
dplm1 <- dpl[,-1L, drop = FALSE]
dpum1 <- dpu[,-1L, drop = FALSE]
prbm1 <- prb[,-1L, drop = FALSE]

i1 <- length(intercepts) - 1L
i2 <- 1L
@}

The off-diagonal elements of $\mA$ are now available as
@d off-diagonal elements for Hessian of intercepts
@{
Aoffdiag <- -rowSums(x * du * dl / prb^2)[-i2]
Aoffdiag <- Aoffdiag[-length(Aoffdiag)]
@}

and the diagonal elements of $\mA$ as
@d diagonal elements for Hessian of intercepts
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

@d intercept / shift contributions to Hessian
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

We return the three matrices $\mA$, $\X$, and $\Z$ necessary for the
computation of the Fisher information for $\beta_2, \dots, \beta_K$ as the Schur
complement $\Z - \X^\top \mA^{-1} \X$. Because the matrix $\mA$ is symmetric
tridiagonal, we use infrastructure from the \pkg{Matrix} package to
represent this matrix:

@d Hessian
@{
.hes <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<parm to prob@>

    @<Hessian prep@>

    @<off-diagonal elements for Hessian of intercepts@>
    @<diagonal elements for Hessian of intercepts@>
    @<intercept / shift contributions to Hessian@>

    if (length(Adiag) > 1L) {
        if(is.null(tryCatch(loadNamespace("Matrix"), error = function(e)NULL)))
                stop(gettextf("%s needs package 'Matrix' correctly installed",
                              "free1way.test"),
                     domain = NA)
        A <- Matrix::bandSparse(length(Adiag), k = 0:1, diagonals = list(Adiag, Aoffdiag), 
                                symmetric = TRUE)
    } else {
        A <- matrix(Adiag)
    }
    return(list(A = A, X = X, Z = Z))
}
@}

We start with an example involving $K = 3$ groups for a binary outcome and
use a binary logistic regression model to estimate the two log-odds ratios
$\beta_2$ and $\beta_3$ along with their estimated covariance
<<glm>>=
library("free1way")
(x <- matrix(c(10, 5, 7, 11, 8, 9), nrow = 2))
d <- expand.grid(y = relevel(gl(2, 1), "2"), t = gl(3, 1))
d$x <- c(x)
m <- glm(y ~ t, data = d, weights = x, family = binomial())
(cf <- coef(m))
@@

Replicating these results requires specification of the inverse link
function $F = \text{expit}$ and the density function $f$ of the standard
logistic. Note that \code{glm} operates with a positive linear predictor, so
we need to change the sign of the log-odds ratios:

<<glm-op>>=
F <- plogis
f <- dlogis
op <- optim(par = c("mt2" = 0, "mt3" = 0, "(Intercept)" = 0), 
            fn = .nll, gr = .nsc, 
            x = x, method = "BFGS", hessian = TRUE)
cbind(c(cf[-1] * -1, cf[1]), op$par)
logLik(m)
-op$value
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
### analytical covariance of parameters
solve(H$Z - crossprod(H$X, solve(H$A, H$X)))
### numerical covariance
solve(op$hessian)[1:2,1:2]
### from glm
vcov(m)[-1,-1]
@@
Also here we see practically identical results. We will later implement a
low-level function \code{.free1way} taking a table and an object describing the inverse link
$F$ as arguments; these results are also in line with \code{glm}:
<<glm-free1way>>=
obj <- .free1wayML(as.table(x), link = logit())
obj$coefficients
-obj$value
### analytical covariance
obj$vcov
@@

In the next step, we extend our results to the stratified case. We iterate
over all blocks and evaluate the negative log-likelihood for the same values
of the shift parameters but block-specific values of the intercept
parameters. Before we begin, we convert the table $C \times K \times B
(\times 2)$ table \code{x} into a list of non-empty $C^\prime \times K$
tables with non-zero row sums:

@d table2list body
@{
dx <- dim(x)
if (length(dx) == 1L)
    stop("")
if (length(dx) == 2L)
    x <- as.table(array(x, dim = c(dx, 1)))
ms <- c(list(x), lapply(seq_along(dx), function(j) marginSums(x, j) > 0))
ms$drop <- FALSE
x <- do.call("[", ms)
dx <- dim(x)
stopifnot(length(dx) >= 3L)
K <- dim(x)[2L]
B <- dim(x)[3L]
stopifnot(dx[1L] > 1L)
stopifnot(K > 1L)
xrc <- NULL
if (length(dx) == 4L) {
    if (dx[4] == 2L) {
        xrc <- array(x[,,,"FALSE", drop = TRUE], dim = dx[1:3])
        x <- array(x[,,,"TRUE", drop = TRUE], dim = dx[1:3])
    } else {
        stop("")
    }
}

xlist <- xrclist <- vector(mode = "list", length = B)

lwr <- rep(-Inf, times = K - 1)
for (b in seq_len(B)) {
    xb <- matrix(x[,,b, drop = TRUE], ncol = K)
    xw <- rowSums(abs(xb)) > 0
    ### do not remove last parameter if there are corresponding
    ### right-censored observations
    if (!is.null(xrc) && any(xrc[dx[1],,b,drop = TRUE] > 0))
        xw[length(xw)] <- TRUE
    if (sum(xw) > 1L) {
        xlist[[b]] <- xb[xw,,drop = FALSE]
        attr(xlist[[b]], "idx") <- xw
        if (!is.null(xrc)) {
            xrclist[[b]] <- matrix(xrc[xw,,b,drop = TRUE], ncol = K)
            attr(xrclist[[b]], "idx") <- xw
        }
    }
}
strata <- !sapply(xlist, is.null)
@}

@d table2list
@{
.table2list <- function(x) {

    @<table2list body@>

    ret <- list(xlist = xlist[strata])
    if (!is.null(xrc))
        ret$xrclist <- xrclist[strata]
    ret$strata <- strata
    ret
}
@}

We first extract the shift parameters $\beta_{\cdot}$ and then, separately
for each stratum, the corresponding contrasts of the intercept parameters:

@d stratum prep
@{
C <- sapply(x, NROW) ### might differ by stratum
K <- unique(do.call("c", lapply(x, ncol))) ### the same
B <- length(x)
sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), levels = seq_len(B))
bidx <- seq_len(K - 1L)
beta <- parm[bidx]
intercepts <- split(parm[-bidx], sidx)
@}

before we loop over the non-empty strata and return the sum of the
corresponding log-likelihoods:

@d stratified negative logLik
@{
.snll <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<stratum prep@>
    ret <- 0
    for (b in seq_len(B))
        ret <- ret + .nll(c(beta, intercepts[[b]]), x[[b]], mu = mu,
                          rightcensored = rightcensored)
    return(ret)
}
@}

In a similar way, we evaluate the gradients for each block and sum-up the
contributions by the shift parameters whereas the gradients for the
intercept parameters are only concatenated:

@d stratified negative score
@{
.snsc <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<stratum prep@>
    ret <- numeric(length(bidx))
    for (b in seq_len(B)) {
        nsc <- .nsc(c(beta, intercepts[[b]]), x[[b]], mu = mu,
                    rightcensored = rightcensored)
        ret[bidx] <- ret[bidx] + nsc[bidx]
        ret <- c(ret, nsc[-bidx])
    }
    return(ret)
}
@}

The score residuum is zero for an observation with weight zero, that is, a
row of zeros in the table:

@d stratified negative score residual
@{
.snsr <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<stratum prep@>
    ret <- c()
    for (b in seq_len(B)) {
        idx <- attr(x[[b]], "idx")
        sr <- numeric(length(idx))
        sr[idx] <- .nsr(c(beta, intercepts[[b]]), x[[b]], mu = mu,
                        rightcensored = rightcensored)
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
m <- glm(y ~ 0 + s + t, data = d, weights = x, family = binomial())
logLik(m)
(cf <- coef(m))
@@

<<glm-op-stratum>>=
xl <- .table2list(x)$xlist
op <- optim(par = c("mt2" = 0, "mt3" = 0, "(Intercept 1)" = 0, "(Intercept 2)" = 0), 
            fn = .snll, gr = .snsc, 
            x = xl, 
            method = "BFGS", 
            hessian = TRUE)
cbind(c(cf[-(1:2)] * -1, cf[1:2]), op$par)
logLik(m)
-op$value
@@

For the analytical Hessian, we sum-up over the stratum-specific
Hessians of the shift parameters. For right-censored observations, we need
to compute the contributions by the events and obtain the joint Hessian for
shift- and intercept parameters first:

@d stratified Hessian
@{
.shes <- function(parm, x, mu = 0, xrc = NULL) {
    @<stratum prep@>
    ret <- matrix(0, nrow = length(bidx), ncol = length(bidx))
    for (b in seq_len(B)) {
        H <- .hes(c(beta, intercepts[[b]]), x[[b]], mu = mu)
        if (!is.null(xrc)) {
            Hrc <- .hes(c(beta, intercepts[[b]]), xrc[[b]], mu = mu, 
                        rightcensored = TRUE)
            H$X <- H$X + Hrc$X
            H$A <- H$A + Hrc$A
            H$Z <- H$Z + Hrc$Z
        }
        ret <- ret + (H$Z - crossprod(H$X, solve(H$A, H$X)))
    }
    ret
}
@}

<<glm-H-stratum>>=
### analytical covariance of parameters
solve(.shes(op$par, xl))
### numerical covariance
solve(op$hessian)[1:2,1:2]
### from glm
vcov(m)[-(1:2),-(1:2)]
@@

<<glm-free1way-strata>>=
obj <- .free1wayML(as.table(x), link = logit())
obj$coefficients
-obj$value
### analytical covariance
obj$vcov
@@
	

\chapter{Link Functions}
\label{ch:link}

Similar to \code{family} objects, we provide some infrastructure for
\code{link} functions $F^{-1}$ and derived quantities (\code{linkinv} $F$,
\code{dlinkinv} $f$, and \code{ddlinkinv} $f^\prime$). If not provided, we also
set-up the ratio $f^\prime / f$ in the constructor.

Although there is some overlap with \code{family} objects for binomial
outcomes, it doesn't seem beneficial to extend this richer class.

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

We start with the logit link, that is $F(z) = (1 + \exp(-z))^{-1}$, giving rise
to Wilcoxon or Kruskal-Wallis type score residuals:

@d logit
@{
logit <- function()
    linkfun(alias = c("Wilcoxon", "Kruskal-Wallis"),
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

The \code{parm2PI} function converts log-odds ratios to probabilistic
indices (or AUCs) and the inverse operation is implemented by
\code{PI2parm}. The overlap coefficient can be obtained from a log-odds
ratio via \code{parm2OVL}.

The log-log link, with $F(z) = \exp(-\exp(-z))$, is used to construct tests
against Lehmann alternatives:

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

The complementary log-log link, with $F(z) = 1 - \exp(-\exp(z))$, provides
log-rank or Savage score residuals against proportional hazards
alternatives:

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

The probit link, with $F(z) = \Phi$, leads to normal scores tests, where the
shift effect can be interpreted as a generalised version of Cohen's $d$,
that is, differences on a latent normal scale with variance one:

@d probit
@{
probit <- function()
    linkfun(alias = "van der Waerden normal scores",
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

@o free1way.R -cp
@{
@<ML estimation@>
@<free1way@>
@<free1way methods@>
@<free1way print@>
@<free1way summary@>
@<free1way confint@>
@<free1way formula@>
@<free1way numeric@>
@<free1way factor@>
@<ppplot@>
@<r2dsim@>
@<power@>
@}

We now put together a low-level function for parameter estimation and
evaluation of scores, Hessians, and residuals. We also set-up a profile
likelihood function for later re-use. 

Assuming all shift effects been zero, we compute starting values for the
intercept parameters from the empirical cumulative distribution function
after merging all treatment groups:

@d setup and starting values
@{
@<table2list body@>
if (NS <- is.null(start))
    start <- rep.int(0, K - 1)
lwr <- rep(-Inf, times = K - 1)
for (b in seq_len(length(xlist))) {
    lwr <- c(lwr, -Inf, rep.int(tol, times = nrow(xlist[[b]]) - 2L))
    if (NS) {
        ecdf0 <- cumsum(rowSums(xlist[[b]]))
        ecdf0 <- ecdf0[-length(ecdf0)] / ecdf0[length(ecdf0)]
        Qecdf <- Q(ecdf0)
        start <- c(start, Qecdf[1], diff(Qecdf))
        start[!is.finite(start)] <- 0
    }
}
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
                     ret <- .snll(p, x = xlist, mu = mu)
                     if (!is.null(xrc))
                         ret <- ret + .snll(p, x = xrclist, mu = mu, 
                                            rightcensored = TRUE)
                     ret
                 },
                 gr = function(par) {
                     p <- numeric(length(par) + length(fix))
                     p[fix] <- beta
                     p[-fix] <- par
                     ret <- .snsc(p, x = xlist, mu = mu)[-fix]
                     if (!is.null(xrc))
                         ret <- ret + .snsc(p, x = xrclist, mu = mu, 
                                            rightcensored = TRUE)[-fix]
                     ret
                 },
                 lower = lwr[-fix], method = "L-BFGS-B", 
                 hessian = FALSE, ...)
    p <- numeric(length(start))
    p[fix] <- beta
    p[-fix] <- ret$par
    ret$par <- p
    ret
}
@}

The heart of the function is a call to \code{optim}, trying to obtain
parameter estimates of $\thetavec$ by minimising the negative
log-likelihood. We allow some (or all) parameters to be fixed at some
constants, and provide a profile version of the likelihood:

@d optim
@{
fn <- function(par) {
    ret <- .snll(par, x = xlist, mu = mu)
    if (!is.null(xrc))
        ret <- ret + .snll(par, x = xrclist, mu = mu, 
                           rightcensored = TRUE)
    return(ret)
}
gr <- function(par) {
    ret <- .snsc(par, x = xlist, mu = mu)
    if (!is.null(xrc))
        ret <- ret + .snsc(par, x = xrclist, mu = mu, 
                           rightcensored = TRUE)
    return(ret)
}
if (!length(fix)) {
    ret <- optim(par = start, fn = fn, gr = gr,
                 lower = lwr, method = "L-BFGS-B", 
                 hessian = FALSE, ...)
} else if (length(fix) == length(start)) {
    ret <- list(par = start, 
                value = fn(start))
} else {
    ret <- .profile(start, fix = fix)
}
@}

After parameter estimation, we evaluate scores, the Hessian, and residuals
as requested:

@d post processing
@{
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
    if (!is.null(xrc))
        ret$negscore <- ret$negscore + .snsc(ret$par, x = xrclist, mu = mu, 
                                             rightcensored = TRUE)[parm]
if (hessian) {
    if (!is.null(xrc)) {
        ret$hessian <- .shes(ret$par, x = xlist, mu = mu, xrc = xrclist)
    } else {
        ret$hessian <- .shes(ret$par, x = xlist, mu = mu)
    }
    if (length(parm) != nrow(ret$hessian))
       ret$hessian <- solve(ret$vcov <- solve(ret$hessian)[parm,parm])
    ret$vcov <- solve(ret$hessian)
    rownames(ret$vcov) <- colnames(ret$vcov) <- rownames(ret$hessian) <-
        colnames(ret$hessian) <-  cnames
}
if (residuals) {
    ret$residuals <- .snsr(ret$par, x = xlist, mu = mu)
    if (!is.null(xrc)) {
        rcr <- .snsr(ret$par, x = xrclist, mu = mu, rightcensored = TRUE)
        C <- sapply(xlist, NROW) 
        ret$residuals <- c(rbind(matrix(ret$residuals, nrow = C),
                                 matrix(rcr, nrow = C)))
     }
}
ret$profile <- function(start, fix)
    .free1wayML(x, link = link, mu = mu, start = start, fix = fix, tol = tol, 
               ...) 
ret$table <- x
ret$mu <- mu
ret$strata <- strata
names(ret$mu) <- link$parm
@}

Finally, we put everything into one function which returns an object of
class \code{free1wayML} for later use:

@d ML estimation
@{
.free1wayML <- function(x, link, mu = 0, start = NULL, fix = NULL, 
                        residuals = TRUE, score = TRUE, hessian = TRUE, 
                        tol = sqrt(.Machine$double.eps), ...) {

    ### convert to three-way table
    stopifnot(is.table(x))
    dx <- dim(x)
    dn <- dimnames(x)
    if (length(dx) == 2L) {
        x <- as.table(array(c(x), dim = dx <- c(dx, 1L)))
        dimnames(x) <- dn <- c(dn, list(A = "A"))
    }

    ### short-cuts for link functions
    F <- function(q) .p(link, q = q)
    Q <- function(p) .q(link, p = p)
    f <- function(q) .d(link, x = q)
    fp <- function(q) .dd(link, x = q)

    @<setup and starting values@>
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
    @<optim@> 
    @<post processing@>

    class(ret) <- "free1wayML"
    ret
}
@}

As an example, consider a stratified (two stata) $3 \times 3$ problem where
outcome category B is missing from the second stratum:

<<workhorse>>=
N <- 10
a <- matrix(c(5, 6, 4,
                    3, 5, 7,
                    3, 4, 5,
                    3, 5, 6,
                    0, 0, 0,
                    4, 6, 5), ncol = 3, byrow = TRUE)
x <- as.table(array(c(a[1:3,], a[-(1:3),]), dim = c(3, 3, 2)))
x
ret <- .free1wayML(x, logit())
ret[c("value", "par")]
cf <- ret$par
cf[1:2] <- cf[1:2] + .5
cf
### profile for cf[1:2]
.free1wayML(x, logit(), start = cf, fix = 1:2)[c("value", "par")]
### profile for cf[2]
.free1wayML(x, logit(), start = cf, fix = 2)[c("value", "par")]
### evaluate log-likelihood at cf
.free1wayML(x, logit(), start = cf, 
            fix = seq_along(ret$par))[c("value", "par")]
@@

\chapter{ML Inference}
\label{ch:MLinf}

Based on an object of class \code{free1wayML}, we can setup different test
statistics and obtain the limiting null distribution based on classical ML
theory under the population model:

@d statistics
@{
if (test == "Wald") {
    @<Wald statistic@>
} else if (test == "LRT") {
    @<LRT@>
} else if (test == "Rao") {
    @<Rao@>
} else if (test == "Permutation") {
    @<Permutation@>
}
@}


\section{Wald}

We only need access to the parameter estimates $\hat{\beta}_2, \dots,
\hat{\beta}_K$ and the corresponding Hessian:

@d Wald statistic
@{
if (alternative == "two.sided") {
    STATISTIC <- c("Wald chi-squared" = c(crossprod(cf, x$hessian %*% cf)))
    DF <- c("df" = length(parm))
    PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
} else {
    STATISTIC <- c("Wald Z" = c(cf * sqrt(c(x$hessian))))
    PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
}
@}

\section{Likelihood-ratio}

In addition to the log-likelihood evaluated at the ML estimates, we need to
evaluate the profile log-likelihood at some value corresponding the null
hypothesis to be tested:

@d LRT
@{
par <- x$par
par[parm] <- value
unll <- x$value ### neg logLik
rnll <- x$profile(par, parm)$value ### neg logLik
STATISTIC <- c("logLR chi-squared" = - 2 * (unll - rnll))
DF <- c("df" = length(parm))
PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
@}

\section{Rao Score}

For the Rao score test, the inverse of the Hessian as well as the score
function of the shift parameters evaluated for some null values need to be
computed:

@d Rao
@{
par <- x$par
par[parm] <- value
ret <- x$profile(par, parm)
if (alternative == "two.sided") {
    STATISTIC <- c("Rao chi-squared" = c(crossprod(ret$negscore, ret$vcov %*% ret$negscore)))
    DF <- c("df" = length(parm))
    PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
} else {
    STATISTIC <- c("Rao Z" = -ret$negscore * sqrt(c(ret$vcov)))
    PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
}
@}

\chapter{Permutation Inference}
\label{ch:Perminf}

Under the permutation model, that is, in randomised experiments where the
random treatment allocation is the only relevant source of randomness, we
compute a permuation variant of the Rao score test, based on the conditional
asymptotic distribution or based on a Monte-Carlo estimate of the reference
distribution:

@d Permutation
@{
par <- x$par
par[parm] <- value
ret <- x$profile(par, parm)
sc <- -ret$negscore
if (length(cf) == 1L)
   sc <- sc / sqrt(c(ret$hessian))
Esc <- sc - x$perm$Expectation
if (alternative == "two.sided" && length(cf) > 1L) {
    STATISTIC <- c("Perm chi-squared" = sum(Esc %*% solve(x$perm$Covariance) * Esc))
    ps <- x$perm$permStat
    if (!is.null(x$perm$permStat))
        PVAL <- mean(ps > STATISTIC + tol)
    else {
        DF <- c("df" = x$perm$DF)
        PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
    }
} else {
    STATISTIC <- c("Perm Z" = Esc / sqrt(c(x$perm$Covariance)))
    if (!is.null(x$perm$permStat)) {
        if (alternative == "two.sided")
            PVAL <- mean(abs(x$perm$permStat) > abs(STATISTIC) + tol)
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
@}

The mean and variance of the linear permutation statistic under the null was
given by \cite{strasserweber1999}:

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
@}

For small samples, we used the \code{r2dtable} function to sample from
tables with fixed marginal distributions:

@d resampling
@{
.resample <- function(res, xt, B = 10000) {

    if (length(dim(xt)) == 2L)
        xt <- as.table(array(xt, dim = c(dim(xt), 1)))

    res <- matrix(res, nrow = dim(xt)[1L], ncol = dim(xt)[3L])
    stat <- 0
    ret <- .SW(res, xt)
    if (dim(xt)[2L] == 2L) {
        ret$testStat <- c((ret$Statistic - ret$Expectation) / sqrt(c(ret$Covariance)))
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
             ret$permStat <- (stat - ret$Expectation) / sqrt(c(ret$Covariance))
        } else {
            ES <- t(matrix(stat, ncol = B) - ret$Expectation)
            ret$permStat <- rowSums(ES %*% solve(ret$Covariance) * ES)
        }
    }
    ret
}
@}

As an example, consider the Wilcoxon rank sum test, where the scores under
the null are a linear function of the ranks of the data. We compute the
asymptotic and approximated reference distribution and corresponding
p-values for a test statistics in quadratic form:

<<SW>>=
w <- gl(2, 15)
(s <- .SW(r <- rank(u <- runif(length(w))), model.matrix(~ 0 + w)))
ps <- .resample(r, model.matrix(~ 0 + w), B = 100000)
ps$testStat^2
mean(abs(ps$permStat) > abs(ps$testStat) - .Machine$double.eps)
pchisq(ps$testStat^ifelse(ps$DF == 1, 2, 1), df = ps$DF, lower.tail = FALSE)
### exactly the same
kruskal.test(u ~ w)
library("coin")
### almost the same
kruskal_test(u ~ w, distribution = approximate(100000))
@@


\chapter{Distribution-free Tests in Stratified $K$-sample Oneway Layouts}

\section{\code{free1way.test}}

We provide a new test procedure in a generic \code{free1way.test}, featuring
a method for tables (the main workhorse) and additional user interfaces. 

@d link2fun
@{
if (!inherits(link, "linkfun")) {
    link <- match.arg(link)
    link <- do.call(link, list())
}
@}

@d free1way
@{
free1way.test <- function(y, ...)
    UseMethod("free1way.test")

free1way.test.table <- function(y, link = c("logit", "probit", "cloglog", "loglog"), 
                                mu = 0, B = 0, ...)
{

    cl <- match.call()

    d <- dim(y)
    dn <- dimnames(y)
    DNAME <- NULL
    if (!is.null(dn)) {
        DNAME <- paste(names(dn)[1], "by", names(dn)[2], 
                       paste0("(", paste0(dn[2], collapse = ", "), ")"))
        if (length(dn) == 3L)
            DNAME <- paste(DNAME, "\n\t stratified by", names(dn)[3])
    }

    @<link2fun@>

    ret <- .free1wayML(y, link = link, mu = mu, ...)
    ret$link <- link
    ret$data.name <- DNAME
    ret$call <- cl

    alias <- link$alias
    if (length(link$alias) == 2L) alias <- alias[1L + (d[2] > 2L)]
    ret$method <- paste(ifelse(d[3L] > 1L, "Stratified", ""), 
                        paste0(d[2L], "-sample"), alias, 
                        "test against", link$model, "alternatives")

    cf <- ret$par
    cf[idx <- seq_len(d[2L] - 1L)] <- 0
    pr <- ret$profile(cf, idx)
    if (d[2L] == 2L)
        res <- pr$residuals / sqrt(c(pr$hessian))
    else
        res <- pr$residuals

    @<Strasser Weber@>
    @<resampling@>

    if (length(dim(y)) == 3L) y <- y[,,ret$strata, drop = FALSE]
    if (length(dim(y)) == 4L) {
        y <- y[,,ret$strata,, drop = FALSE]
        dy <- dim(y)
        dy[1] <- dy[1] * 2
        y <- apply(y, 3, function(x) rbind(x[,,2], x[,,1]))
        y <- array(y, dim = dy[1:3])
    }
    ret$perm <- .resample(res, y, B = B)

    if (!is.null(names(dn))) {
        fm <- as.formula(paste(names(dn)[1:2], collapse = "~"))
        ret$terms <- terms(fm, data = as.data.frame(y))
    }

    class(ret) <- "free1way"
    return(ret)
}
@}

The \code{formula} method allows formulae \code{outcome ~ treatment +
strata(s)} for model specification

<TH>strata is only defined in \pkg{survival}, import? </TH>

@d free1way formula
@{
free1way.test.formula <- function(formula, data, weights, subset, na.action = na.pass, ...)
{

    cl <- match.call()

    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")

    strata <- function(object) object
    formula <- terms(formula, specials = "strata")

    stratum <- attr(formula, "specials")$strata
    if (is.null(stratum)) stratum <- 0L
    
    if (length(attr(formula, "term.labels")) > 1L + stratum)
        stop("'formula' missing or incorrect")
    group <- attr(formula, "term.labels") 
    if (stratum) group <- group[-(stratum - 1L)]

    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    DNAME <- paste(vn <- c(names(mf)[response], group), collapse = " by ") # works in all cases
    w <- as.vector(model.weights(mf))
    y <- mf[[response]]
    g <- mf[[group]]
    stopifnot(is.factor(g))
    lev <- levels(g)
    DNAME <- paste(DNAME, paste0("(", paste0(lev, collapse = ", "), ")"))
    if (nlevels(g) < 2L)
        stop("grouping factor must have at least 2 levels")
    if (stratum) {
        st <- factor(mf[[stratum]], levels = )
        if (nlevels(st) < 2L)
            stop("at least two strata must be present")
        vn <- c(vn, names(mf)[stratum])
        RVAL <- free1way.test(y = y, x = g, z = st, weights = w, 
                              varnames = vn, ...)
        DNAME <- paste(DNAME, paste("\n\t stratified by", names(mf)[stratum]))
    } else {
        ## Call the corresponding method
        RVAL <- free1way.test(y = y, x = g, weights = w, varnames = vn, ...)
    }
    RVAL$data.name <- DNAME
    RVAL$call <- cl
    RVAL
}
@}

The method for numeric outcomes provides a discretisation at the unique
observed outcome values, or (for very large sample sizes), for binned
outcomes. The \code{event} argument is a logical where \code{TRUE} is
interpreted as an event and \code{FALSE} as right-censored observation

<TH>add event to formula interface</TH>

@d free1way numeric
@{
free1way.test.numeric <- function(y, x, z = NULL, event = NULL, weights = NULL, nbins = 0, 
    varnames = c(deparse1(substitute(y)), 
                 deparse1(substitute(x)), 
                 deparse1(substitute(z))), ...) {

    cl <- match.call()
    DNAME <- paste(varnames[1], "by", varnames[2])
    DNAME <- paste(DNAME, paste0("(", paste0(levels(x), collapse = ", "), ")"))

    if (!is.null(z))
        DNAME <- paste(DNAME, "\n\t stratified by", varnames[3])
    varnames <- varnames[varnames != "NULL"]

    if (!is.null(event)) {
        stopifnot(is.logical(event))
        uy <- sort(unique(y[event]))
        if (all(y[!event] < uy[length(uy)]))
            uy <- uy[-length(uy)]
    } else {
        uy <- sort(unique(y))
    }
    if (nbins && nbins < length(uy)) {
        nbins <- ceiling(nbins)
        breaks <- c(-Inf, quantile(y, probs = seq_len(nbins) / (nbins + 1L)), Inf)
    } else {
        breaks <- c(-Inf, uy, Inf)
    }
    r <- cut(y, breaks = breaks, ordered_result = TRUE)[, drop = TRUE]
    RVAL <- free1way.test(y = r, x = x, z = z, event = event, weights = weights, 
                          varnames = varnames, ...)
    RVAL$data.name <- DNAME
    RVAL$call <- cl
    RVAL
}
@}

The \code{factor} method also allows right-censoring but otherwise is just a
call to \code{xtabs}:

@d free1way factor
@{
free1way.test.factor <- function(y, x, z = NULL, event = NULL, weights = NULL, 
    varnames = c(deparse1(substitute(y)), 
                 deparse1(substitute(x)), 
                 deparse1(substitute(z))), ...) {

    cl <- match.call()
    DNAME <- paste(varnames[1], "by", varnames[2])
    DNAME <- paste(DNAME, paste0("(", paste0(levels(x), collapse = ", "), ")"))

    if (!is.null(z))
        DNAME <- paste(DNAME, "\n\t stratified by", varnames[3])
    varnames <- varnames[varnames != "NULL"]

    stopifnot(is.factor(x))
    if (nlevels(y) > 2L)
        stopifnot(is.ordered(y))
    d <- data.frame(w = 1, y = y, x = x)
    if (!is.null(weights)) d$w <- weights
    if (is.null(z)) z <- gl(1, nrow(d))
    d$z <- z 
    if (!is.null(event)) {
        stopifnot(is.logical(event))
        d$event <- event
    }
    tab <- xtabs(w ~ ., data = d)
    dn <- dimnames(tab)
    names(dn)[seq_along(varnames)] <- varnames
    dimnames(tab) <- dn
    RVAL <- free1way.test(tab, ...)
    RVAL$data.name <- DNAME
    RVAL$call <- cl
    RVAL
}
@}

\section{\code{free1way} Methods}

We start with \code{coef}, \code{vcov}, and
\code{model.frame}/\code{model.matrix} methods such that multiple comparison
procedures from \pkg{multcomp} will work out of the box. The \code{coef}
method allows to obtain effects at alternative scales: probabilistic indices
(\code{AUC} = \code{PI}) or the overlap coefficient:

@d free1way methods
@{
coef.free1way <- function(object, what = c("shift", "PI", "AUC", "OVL"), ...)
{
    what <- match.arg(what)
    cf <- object$coefficients
    return(switch(what, "shift" = cf,
                        "PI" = object$link$parm2PI(cf),
                        "AUC" = object$link$parm2PI(cf),	### same as PI
                        "OVL" = object$link$parm2OVL(cf)))
}
vcov.free1way <- function(object, ...)
    object$vcov
logLik.free1way <- function(object, ...)
    -object$value
### the next two could go into multcomp
model.frame.free1way <- function(formula, ...)
    as.data.frame(formula$table)
model.matrix.free1way <- function (object, ...) 
{
    mm <- model.matrix(delete.response(terms(object)), data = model.frame(object))
    at <- attributes(mm)
    mm <- mm[, -1]
    at$dim[2] <- at$dim[2] - 1
    at$dimnames[[2]] <- at$dimnames[[2]][-1]
    at$assign <- at$assign[-1]
    attributes(mm) <- at
    mm
}
@}

We use the \code{print} method to report different test statistics and
corresponding $p$-values via the \code{test} and \code{alternative}
arguments. The reason for doing so is that the parameter estimation only
needs to be performed once in cases users are interested in different
tests or (see below) confidence intervals. By default, an asymptotic  
permutation test is performed, mainly because the $p$-values coincide with
some special cases (\code{wilcox,kruskal,friedman.test}):

@d free1way print
@{
.print.free1way <- function(x, test = c("Permutation", "Wald", "LRT", "Rao"), 
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
        null.value = x$mu, alternative = alternative, method = x$method, 
        data.name = x$data.name)
    class(RVAL) <- "htest"
    return(RVAL)
}

print.free1way <- function(x, ...) {
    print(ret <- .print.free1way(x, ...))
    return(invisible(x))
}
@}

The \code{summary} method performs population Wald inference unless the
\code{test} argument is specified:

@d free1way summary
@{
summary.free1way <- function(object, test, alternative = c("two.sided", "less", "greater"), 
                             tol = .Machine$double.eps, ...)
{

    if (!missing(test))
        return(.print.free1way(object, test = test, alternative = alternative, tol = tol))
   
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
    ret <- list(call = object$call, coefficients = cfmat)
    class(ret) <- "summary.free1way"
    return(ret)
}
print.summary.free1way <- function(x, ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Coefficients:\n")
    printCoefmat(x$coefficients)
}
@}

Confidence intervals are computed by inversion of the corresponding test
statistics. Because LRT and Rao confidence intervals are invariant wrt to
transformations, proper LRT or Rao confidence intervals for probabilistic
indices or overlap coefficients can also be computed. The computation always
starts with Wald intervals, which are either returned or used as starting
values for the inversion:

@d free1way confint
@{
confint.free1way <- function(object, parm,
    level = .95, test = c("Permutation", "Wald", "LRT", "Rao"), 
    what = c("shift", "PI", "AUC", "OVL"), ...)
{

    test <- match.arg(test)
    conf.level <- 1 - (1 - level) / 2

    cf <- coef(object)
    if (missing(parm)) 
        parm <- seq_along(cf)

    CINT <- confint.default(object, level = level)
    if (test == "Wald")
        return(CINT)
    wlevel <- level
    wlevel <- 1 - (1 - level) / 10
    CINT[] <- confint.default(object, level = wlevel)

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
            qu <- quantile(object$perm$permStat, 
                           probs = c(1 - conf.level, conf.level))
            att.level <- mean(object$perm$permStat > qu[1] & 
                              object$perm$permStat < qu[2])
            attr(CINT, "Attained level") <- att.level
        }
    } else {
        qu <- rep.int(qchisq(level, df = 1), 2) ### always two.sided
    }

    for (p in parm) {
        CINT[p, 1] <- uniroot(sfun, interval = c(CINT[p,1], cf[p]), 
                              parm = p, quantile = qu[2])$root
        CINT[p, 2] <- uniroot(sfun, interval = c(cf[p], CINT[p, 2]), 
                              parm = p, quantile = qu[1])$root
    }

    what <- match.arg(what)
    CINT <- switch(what, "shift" = CINT,
                         "PI" = object$link$parm2PI(CINT),
                         "AUC" = object$link$parm2PI(CINT), ### same as PI 
                         "OVL" = object$link$parm2OVL(CINT))
    return(CINT)
}
@}

As an example, we compute log-odds ratios for the table introduced above and
report some tests and confidence intervals:

<<free>>=
x
ft <- free1way.test(x)
coef(ft)
vcov(ft)
### Wald per parameter
summary(ft)
library("multcomp")
summary(glht(ft), test = univariate())

### global Wald
summary(ft, test = "Wald")
summary(glht(ft), test = Chisqtest())

### Rao score, Permutation score, LRT
summary(ft, test = "Rao")
summary(ft, test = "Permutation")
summary(ft, test = "LRT")

### Wald confidence intervals, unadjusted
confint(glht(ft), calpha = univariate_calpha())
confint(ft, test = "Wald")

### Rao and LRT intervals
confint(ft, test = "Rao")
confint(ft, test = "LRT")
@@

\section{Wilcoxon Test}

The second example is a Wilcoxon test for a single log-odds ratio
comparing to treatment groups:

<<formula>>=
set.seed(29)
N <- 25
w <- gl(2, N)
y <- rlogis(length(w), location = c(0, 1)[w])

#### link = logit is default
ft <- free1way.test(y ~ w)

### Wald 
summary(ft)

### Permutation test
wilcox.test(y ~ w, alternative = "greater", correct = FALSE)$p.value
pvalue(wilcox_test(y ~ w, alternative = "greater"))
summary(ft, test = "Permutation", alternative = "less")$p.value
wilcox.test(y ~ w, alternative = "less", correct = FALSE)$p.value
pvalue(wilcox_test(y ~ w, alternative = "less"))
summary(ft, test = "Permutation", alternative = "greater")$p.value
wilcox.test(y ~ w, correct = FALSE)$p.value
kruskal.test(y ~ w)$p.value
pvalue(wilcox_test(y ~ w))
summary(ft, test = "Permutation")$p.value

### Wald tests
summary(ft, test = "Wald", alternative = "less")
summary(ft, test = "Wald", alternative = "greater")
summary(ft, test = "Wald")

### Rao score tests
summary(ft, test = "Rao", alternative = "less")
summary(ft, test = "Rao", alternative = "greater")
summary(ft, test = "Rao")

### LRT (only two-sided)
summary(ft, test = "LRT")

### confidence intervals for log-odds ratios
confint(ft, test = "Permutation")
confint(ft, test = "LRT")
confint(ft, test = "Wald")
confint(ft, test = "Rao")

### confidence interval for "Wilcoxon Parameter" = PI = AUC
confint(ft, test = "Rao", what = "AUC")

### comparison with rms::orm
library("rms")
rev(coef(or <- orm(y ~ w)))[1]
coef(ft)
logLik(or)
logLik(ft)
vcov(or)[2,2]
vcov(ft)
ci <- confint(or)
ci[nrow(ci),]
confint(ft, test = "Wald")
@@


\section{Mantel-Haenszel test}

<<mh>>=
example(mantelhaen.test, echo = FALSE)
mantelhaen.test(UCBAdmissions, correct = FALSE)
ft <- free1way.test(UCBAdmissions)
summary(ft, test = "Wald")
exp(coef(ft))
exp(confint(ft, test = "Wald"))
exp(sapply(dimnames(UCBAdmissions)[[3L]], function(dept)
       confint(free1way.test(UCBAdmissions[,,dept]), test = "Permutation")))
sapply(dimnames(UCBAdmissions)[[3L]], function(dept)
       fisher.test(UCBAdmissions[,,dept], conf.int = TRUE)$conf.int)
@@

\section{Kruskal-Wallis}

<<kw>>=
example(kruskal.test, echo = FALSE)
kruskal.test(x ~ g)
free1way.test(x ~ g)
@@

\section{Savage / Log-rank}

We start without censoring (Savage test) and add strata

<<sw>>=
library("survival")
N <- 10
nd <- expand.grid(g = gl(3, N), s = gl(3, N))
nd$tm <- rexp(nrow(nd))
nd$ev <- TRUE
survdiff(Surv(tm, ev) ~ g + strata(s), data = nd, rho = 0)$chisq
cm <- coxph(Surv(tm, ev) ~ g + strata(s), data = nd)


(ft <- free1way.test(tm ~ g + strata(s), data = nd, link = "cloglog"))
coef(cm)
coef(ft)
vcov(cm)
vcov(ft)
summary(ft)
summary(cm)$sctest
summary(ft, test = "Rao")
summary(cm)$logtest
summary(ft, test = "LRT")
summary(cm)$waldtest
summary(ft, test = "Wald")
summary(ft, test = "Permutation")

library("coin")
independence_test(Surv(tm, ev) ~ g | s, data = nd, ytrafo = function(...)
                  trafo(..., numeric_trafo = logrank_trafo, block = nd$s), teststat = "quad")
@@

Wilcoxon against proportional odds

<<Peto>>=
survdiff(Surv(tm, ev) ~ g + strata(s), data = nd, rho = 1)$chisq
(ft <- free1way.test(tm ~ g + strata(s), data = nd, link = "logit"))
summary(ft)
summary(ft, test = "Rao")
summary(ft, test = "LRT")
summary(ft, test = "Wald")
summary(ft, test = "Permutation")
@@

\section{van der Waerden}

Normal scores test against a generalised Cohen's $d$:

<<normal>>=
nd$y <- rnorm(nrow(nd))
free1way.test(y ~ g + strata(s), data = nd, link = "probit")
independence_test(y ~ g | s, data = nd, ytrafo = function(...)
                  trafo(..., numeric_trafo = normal_trafo, block = nd$s), teststat = "quad")
@@

\section{Friedman}

Each observation is a block

<<friedman>>=
example(friedman.test, echo = FALSE)
rt <- expand.grid(str = gl(22, 1),
                  trt = gl(3, 1, labels = c("Round Out", "Narrow Angle", "Wide Angle")))
rt$tm <- c(RoundingTimes)
friedman.test(RoundingTimes)
(ft <- free1way.test(tm ~ trt + strata(str), data = rt))
summary(ft)
@@

\chapter{Model Diagnostics}

The classical shift model $F_Y(y \mid T = 2) = F_Y(y - \mu \mid T = 1)$
can be critisised using confidence bands for QQ-plots in \code{qqplot},
because the parameter $\mu$ shows up as a vertical shift of the diagonal
if the model is appropriate.

Likewise, model~(\ref{model}) can be graphically assessed using the PP-plot.
We concentrate on the two-sample case. The shift parameter $\beta_2$ gives
rise to the model-based PP graph $(p, F(F^{-1}(p) - \beta_2))$ and a
confidence \emph{band} can be obtained from a confidence \emph{interval} for
$\beta_2$. The PP-plot is, up to rescalings, identical to the ROC curve.

@d ROC bands
@{
 if (!is.null(conf.level)) {
    prb <- seq_len(1000) / 1001
    res <- c(x, y)
    grp <- gl(2, 1, labels = c(xlab, ylab))
    grp <- grp[rep(1:2, c(length(x), length(y)))]
    args <- conf.args
    args$y <- res
    args$x <- grp
    args$border <- args$col <- args$type <- NULL
    f1w <- do.call("free1way.test", args)

    ci <- confint(f1w, level = conf.level, type = args$type)
    lwr <- .p(f1w$link, .q(f1w$link, prb) - ci[1,1])
    upr <- .p(f1w$link, .q(f1w$link, prb) - ci[1,2])
    x <- c(prb, rev(prb))
    y <- c(lwr, rev(upr))
    xn <- c(x[1L], rep(x[-1L], each = 2))
    yn <- c(rep(y[-length(y)], each = 2), y[length(y)])
    polygon(x = xn, y = yn, col = conf.args$col, border = conf.args$border)
    lines(prb, .p(f1w$link, .q(f1w$link, prb) - coef(f1w)))
}
@}

We introduce a new function \code{ppplot}, closely following the
implementation of \code{qqplot}, allowing to plot the empirical and
corresponding model-based PP-plot, the latter for a certain choice of link
function:

@d ppplot
@{
ppplot <- function(x, y, plot.it = TRUE,
                   xlab = deparse1(substitute(x)),
                   ylab = deparse1(substitute(y)), 
                   ..., conf.level = NULL, 
                   conf.args = list(link = "logit", type = "Wald", 
                                    col = NA, border = NULL)) {

    force(xlab)
    force(ylab)
    if (xlab == ylab) {
        xlab <- paste0("x = ", xlab)
        ylab <- paste0("y = ", ylab)
    }

    ex <- ecdf(x)
    if (interpolate) {
        vals <- sort(unique(x))
        ex <- splinefun(vals, ex(vals), method = "hyman")
    }
    sy <- sort(unique(y))
    py <- ecdf(y)(sy)
    px <- ex(sy)
    ret <- list(x = px, y = py)
    if (!plot.it)
        return(ret)

    plot(px, py, xlim = c(0, 1), ylim = c(0, 1), 
         xlab = xlab, ylab = ylab, type = "n", ...)

    @<ROC bands@>

    points(px, py, ...)
    return(invisible(ret)) 
}
@}

Correct logistic model with log-odds ratio three:

\begin{figure}
<<ppplot, fig = TRUE>>=
y <- rlogis(50)
x <- rlogis(50, location = 3)
ppplot(y, x, conf.level = .95)
@@
\end{figure}

Incorrect proportional hazards alternative:

\begin{figure}
<<ppplot-savage, fig = TRUE>>=
ppplot(y, x, conf.args = list(link = "cloglog", type = "Wald", 
                              col = NA, border = NULL),
       conf.level = .95)
@@
\end{figure}


\chapter{Power and Sample Size}

The term ``distribution-free'' refers to the invariance of the reference
distribution with respect to the distribution of an absolutely continuous
outcome under control. Unfortunately, this is no longer true for
non-continuous outcomes (due to ties) and under the alternative. That means
that sample size assessments always take place under certain assumptions
regarding the outcome distribution.

We start implementing a function for simulating $C \times K$ tables. We need
to specify the number of observations in each treatment group (\code{c}),
the discrete distribution of the control (\code{r}), a model (\code{link}), and a
treatment effect (\code{delta}, in line with \code{power.XYZ.test}). In
essence, we draw samples from the multinomial distribution after computing
the relevant discrete density.

@d r2dsim
@{
r2dsim <- function(n, r, c, delta = 0,
                   link = c("logit", "probit", "cloglog", "loglog")) 
{

    if (length(n <- as.integer(n)) == 0L || (n < 0) || is.na(n)) 
        stop("invalid argument 'n'")
    colsums <- c
    if (length(colsums[] <- as.integer(c)) <= 1L || 
        any(colsums < 0) || anyNA(colsums)) 
        stop("invalid argument 'c'")

    prob <- r
    if (length(prob[] <- as.double(r / sum(r))) <= 1L || 
        any(prob < 0) || anyNA(prob)) 
        stop("invalid argument 'r'")

    if (is.null(names(prob))) 
        names(prob) <- paste0("i", seq_along(prob))
    
    K <- length(colsums)
    if (is.null(names(colsums))) 
        names(colsums) <- LETTERS[seq_len(K)]
    delta <- rep_len(delta, K - 1L)

    @<link2fun@>

    p0 <- cumsum(prob)
    h0 <- .q(link, p0)

    h1 <- h0 - matrix(delta, nrow = length(prob), ncol = K - 1, byrow = TRUE)
    p1 <- .p(link, h1)
    p <- cbind(p0, p1)
    ret <- vector(mode = "list", length = n)

    for (i in seq_len(n)) {
        tab <- sapply(seq_len(K), function(k) 
            rmultinom(1L, size = colsums[k], 
                      prob = c(p[1,k], diff(p[,k]))))
        ret[[i]] <- as.table(array(tab, dim = c(length(prob), K), 
                          dimnames = list(names(prob), 
                                          names(colsums))))
    }
    return(ret)
}
@}

We are now ready to put together a function for power evaluation and sample
size assessment. The core idea is to draw samples from the relevant data
(under a specific model in the alternative) and to estimate the Fisher
information of the treatment effect parameters for this configuration. The
power of the global Wald test can than be approximated by a non-central
$\chi^2$ distribution. This is much faster than approximating the power
directly. Nevertheless, this is a random experiment, so we first make
computations reproducible:

@d random seed
@{
if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
if (is.null(seed)) 
    seed <- RNGstate <- get(".Random.seed", envir = .GlobalEnv)
else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
}
@}

@d power setup
@{

@<link2fun@>

### matrix means control distributions in different strata
if (!is.matrix(prob))
    prob <- matrix(prob, nrow = NROW(prob))
prob <- prop.table(prob, margin = 2L)
C <- nrow(prob)
K <- length(delta) + 1L
B <- ncol(prob)
if (is.null(colnames(prob))) 
    colnames(prob) <- paste0("stratum", seq_len(B))
if (is.null(names(delta))) 
    names(delta) <- LETTERS[seq_len(K)[-1]]
p0 <- apply(prob, 2, cumsum)
h0 <- .q(link, p0)
if (length(alloc_ratio) == 1L) 
    alloc_ratio <- rep_len(alloc_ratio, K - 1)
stopifnot(length(alloc_ratio) == K - 1)
if (length(strata_ratio) == 1L) 
    strata_ratio <- rep_len(strata_ratio, B - 1)
stopifnot(length(strata_ratio) == B - 1)
### sample size per group (columns) and stratum (rows)
N <- n * matrix(c(1, alloc_ratio), nrow = B, ncol = K, byrow = TRUE) * 
         matrix(c(1, strata_ratio), nrow = B, ncol = K)
rownames(N) <- colnames(prob)
ctrl <- "Control"
dn <- dimnames(prob)
if (!is.null(names(dn)[1L]))
    ctrl <- names(dn)[1L]
colnames(N) <- c(ctrl, names(delta))
@}

@d estimate Fisher information
@{
he <- 0
deltamu <- delta - mu
for (i in seq_len(nsim)) {
    parm <- deltamu
    x <- as.table(array(0, dim = c(C, K, B)))
    for (b in seq_len(B)) {
        x[,,b] <- r2dsim(1L, r = prob[, b], c = N[b,], delta = delta, link = link)[[1L]]
        rs <- rowSums(x[,,b]) > 0
        h <- h0[rs, b]
        theta <- c(h[1], diff(h[-length(h)]))
        parm <- c(parm, theta)
    }
    ### evaluate observed hessian for true parameters parm and x data
    he <- he + .free1wayML(x, link = link, mu = mu, start = parm, 
                           fix = seq_along(parm))$hessian
}
### estimate expected Fisher information
he <- he / nsim
@}

The power function now depends on sample size (\code{n}; the number of
control observations in the first stratum), a discrete control distribution
(\code{prob}, this can be a $C \times B$ matrix for stratum-specific control
distributions), a vector of allocation ratios (\code{alloc_ratio = 2} means
control:treatment = 1:2) and the sample size ratios between strata.

The treatment effects are contained in $K - 1$ vector \code{delta}:

@d power call
@{
power.free1way.test(n = n, prob = prob, alloc_ratio = alloc_ratio, 
                    strata_ratio = strata_ratio, delta = delta, mu = mu,
                    sig.level = sig.level, link = link, 
                    alternative = alternative, 
                    nsim = nsim, seed = seed, tol = tol)$power - power
@}

@d power htest output
@{
ss <- paste(colSums(N), paste0("(", colnames(N), ")"), collapse = " + ")
ret <- list(n = n, 
            "Total sample size" = paste(ss, "=", sum(N)),
            power = power, 
            sig.level = sig.level)
if (mu != 0) ret$mu <- mu
ret[[link$parm]] <- delta
ret$note <- "'n' is sample size in control group"
if (B > 1) ret$note <- paste(ret$note, "of first stratum")
alias <- link$alias
if (length(link$alias) == 2L) alias <- alias[1L + (K > 2L)]
ret$method <- paste(ifelse(B > 1L, "Stratified", ""), 
                    paste0(K, "-sample"), alias, 
                    "test against", link$model, "alternatives")
class(ret) <- "power.htest"
@}

@d power
@{
power.free1way.test <- function(n = NULL, prob = rep.int(1 / n, n), 
                                alloc_ratio = 1, strata_ratio = 1, 
                                delta = NULL, mu = 0, sig.level = .05, power = NULL,
                                link = c("logit", "probit", "cloglog", "loglog"),
                                alternative = c("two.sided", "less", "greater"), 
                                nsim = 100, seed = NULL, tol = .Machine$double.eps^0.25) 
{

    if (sum(vapply(list(n, delta, power, sig.level), is.null, 
        NA)) != 1) 
        stop("exactly one of 'n', 'delta', 'power', and 'sig.level' must be NULL")
    stats:::assert_NULL_or_prob(sig.level)
    stats:::assert_NULL_or_prob(power)

    @<random seed@>

    if (is.null(n)) 
        n <- ceiling(uniroot(function(n) {
                 @<power call@>
             }, interval = c(5, 1e+03), tol = tol, extendInt = "upX")$root)
    else if (is.null(delta)) {
        ### 2-sample only
        stopifnot(K == 2L)
        delta <- uniroot(function(delta) {
                 @<power call@>
    ### <TH> interval depending on alternative, symmetry? </TH>
            }, interval = c(0, 10), tol = tol, extendInt = "upX")$root
        }
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) {
                @<power call@>
            }, interval = c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
    
    @<power setup@>
    @<estimate Fisher information@>

    alternative <- match.arg(alternative)
    if (K == 2L) {
        se <- 1 / sqrt(c(he))
        power  <- switch(alternative, 
            "two.sided" = pnorm(qnorm(sig.level / 2) + deltamu / se) + 
                          pnorm(qnorm(sig.level / 2) - deltamu / se),
            "less" = pnorm(qnorm(sig.level) - deltamu / se),
            "greater" = pnorm(qnorm(sig.level) + deltamu / se)
        )
    } else {
        stopifnot(alternative == "two.sided")
        ncp <- sum((chol(he) %*% deltamu)^2)
        qsig <- qchisq(sig.level, df = K - 1L, lower.tail = FALSE)
        power <- pchisq(qsig, df = K - 1L, ncp = ncp, lower.tail = FALSE)
    }

    @<power htest output@>

    ret
}
@}

We start with the power of a binomial experiment with $N = 2 \times 25$
observations. In the control group, the odds of winning is 1. Under
treatment, we increase this odds by $50\%$. We compare the results with
\code{power.prop.test}:

<<power.prop.test>>=
delta <- log(1.5)
power.prop.test(n = 25, p1 = .5, p2 = plogis(qlogis(.5) - delta))
power.free1way.test(n = 25, prob = c(.5, .5), delta = delta)
@@

Under stratification (twice as many observations in the second stratum) 
and with an ordered outcome at four levels, we might want to compare four
groups, with $25\%$, $50\%$, and $75\%$ increase compared to the odds of the
control:
<<power.odds.test>>=
prb <- matrix(c(.25, .25, .25, .25,
                .10, .20, .30, .40), ncol = 2)
colnames(prb) <- c("s1", "s2")
power.free1way.test(n = 20, prob = prb, 
                    strata_ratio = 2,
                    alloc_ratio = c(1.5, 2, 2), 
                    delta = log(c("low" = 1.25, "med" = 1.5, "high" = 1.75)))
@@

We now estimate the power of a Wilcoxon test with, first by simulation from
a logistic distribution, and then by our power function:

<<wilcox>>=
Nsim <- 100
delta <- log(3)
N <- 15
w <- gl(2, N)
pw <- numeric(Nsim)
for (i in seq_along(pw)) {
    y <- rlogis(length(w), location = c(0, delta)[w])
    pw[i] <- wilcox.test(y ~ w)$p.value
}
mean(pw < .05)

power.free1way.test(n = N, delta = delta)
@@

The power of the Kruskal-Wallis test only needs one additional treatment
effect

<<kruskal>>=
delta <- c("B" = log(2), "C" = log(3))
N <- 15
w <- gl(3, N)
pw <- numeric(Nsim)
for (i in seq_along(pw)) {
    y <- rlogis(length(w), location = c(0, delta)[w])
    pw[i] <- kruskal.test(y ~ w)$p.value
}
mean(pw < .05)

power.free1way.test(n = N, delta = delta)
@@

We next use the \code{r2dsim} function to sample from $4 \times 3$ tables with odds ratios $2$ and $3$
and compare the resulting power with result obtained from the approximated
Fisher information. The plot shows the distribution of the parameter
estimates and the corresponding population values as red dots.	

<<table, fig = TRUE>>=
prb <- rep(1, 4)
x <- r2dsim(Nsim, r = prb, c = table(w), delta = delta)
pw <- numeric(length(x))
cf <- matrix(0, nrow = length(x), ncol = length(delta))
colnames(cf) <- names(delta)
for (i in seq_along(x)) {
    ft <- free1way.test(x[[i]])
    cf[i,] <- coef(ft)
    pw[i] <- summary(ft, test = "Permutation")$p.value
}
mean(pw < .05)
boxplot(cf)
points(c(1:2), delta, pch = 19, col = "red")
power.free1way.test(n = N, prob = rep(1, 4), delta = delta)
@@

In the last example, we sample from $4 \times 3$ tables with odds ratios $2$ and $3$ for three
strata with different control distributions, and again compare the
simulation results to the power function:

<<stable, fig = TRUE>>=
prb <- cbind(S1 = rep(1, 4), 
             S2 = c(1, 2, 1, 2), 
             S3 = 1:4)
dimnames(prb) <- list(Ctrl = paste0("i", seq_len(nrow(prb))),
                      Strata = colnames(prb))

x1 <- r2dsim(Nsim, r = prb[, "S1"], c = table(w), delta = delta)
x2 <- r2dsim(Nsim, r = prb[, "S2"], c = table(w), delta = delta)
x3 <- r2dsim(Nsim, r = prb[, "S3"], c = table(w), delta = delta)
stab <- function(...) {
    args <- list(...)
    as.table(array(unlist(args), dim = c(dim(args[[1]]), length(args))))
}
pw <- numeric(length(x1))
cf <- matrix(0, nrow = length(x1), ncol = length(delta))
colnames(cf) <- names(delta)
for (i in seq_along(x)) {
    ft <- free1way.test(stab(x1[[i]], x2[[i]], x3[[i]]))
    cf[i,] <- coef(ft)
    pw[i] <- summary(ft, test = "Permutation")$p.value
}
mean(pw < .05)
boxplot(cf)
points(c(1:2), delta, pch = 19, col = "red")

power.free1way.test(n = N, prob = prb, delta = delta, seed = 3)
power.free1way.test(power = .8, prob = prb, delta = delta, seed = 3)
power.free1way.test(n = 19, prob = prb, delta = delta, seed = 3)
@@

\chapter*{Index}

\section*{Files}

@f

\section*{Fragments}

@m

\section*{Identifiers}

@u

\bibliographystyle{plainnat}
\bibliography{\Sexpr{gsub("\\.bib", "", system.file("refs.bib", package = "free1way"))}}

\end{document}

\chapter{Schur Complement}
\label{ch:schur}

@o Schur.c -cc
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
    .Call(R_wcrossprod, a = as.double(A$a), 
                        b = as.double(A$b),
                        X = x,
                        tol = tol)
}
@}
