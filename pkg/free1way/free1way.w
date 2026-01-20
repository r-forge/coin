\documentclass[a4paper]{report}

%\VignetteIndexEntry{Stratified K-sample Inference}
%\VignetteDepends{free1way,multcomp,survival,Hmisc,coin,rms}
%\VignetteKeywords{semiparametric model,conditional inference}}
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
pageanchor=true,%
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
\newcommand{\rT}{G}
\newcommand{\rS}{S}
\newcommand{\rt}{t}


\author{Torsten Hothorn \\ Universit\"at Z\"urich \and
        Kurt Hornik \\ WU Wirtschaftsuniversit\"at Wien}

\title{Semiparametrically Efficient Population and Permutation Inference in 
       Distribution-free Stratified $K$-sample Oneway Layouts}

\begin{document}

\pagenumbering{roman}
\maketitle
\tableofcontents

\chapter*{Introduction}

Comparing two or more independent samples with respect to some outcome
measure is a common task. Many procedures are available in 
\pkg{stats} and other add-on packages, most of these implementations
making rather strict assumptions regarding the outcome distribution, the
number of samples, the presence of blocks or strata and typically offer
either conditional or unconditional (exact or asymptotic) inference.

This document presents a unified, dense, and yet holistic implementation
covering many classical procedures as special cases. Leveraging 
transformation models, likelihood-based parameter estimation as well as
permutation- and likelihood-based inference are formulated and implemented. One
can understand this contribution as a unification of many of
\code{stats::*.test} procedures, the models available in 
\code{MASS::polr}, \code{rms::orm}, \code{rms::lrm}, \code{survival::coxph}, 
or the \pkg{tram} add-on package (among many others), 
and permutation-based inference in \pkg{coin}.

This implementation is, however, free of any strong dependencies and only
uses functionality available in \proglang{R} itself and the \pkg{stats},
\pkg{graphics}, and \pkg{Matrix} recommended packages.

\chapter{Model and Parameterisation}
\label{ch:model}
\pagenumbering{arabic}

We consider $K$ treatment groups $\rT \in \{1, \dots, K\}, K \ge 2$ for an
at least ordered outcome $Y \in \samY$ observed in
stratum $\rS \in \{1, \dots, B\}$ out of $B \ge 1$ blocks with conditional
cumulative distribution function (cdf)
$F_Y(y \mid \rS = b, \rT = k) = \Prob(Y \le y \mid \rS = b, \rT = k)$. Detecting
and describing differential distributions arising from different treatments
is our main objective. We refer to the first treatment $\rT = 1$ as
``control''.

\paragraph{Model}

With model function $g: [0,1] \times \R \rightarrow [0,1]$, we describe
the conditional distribution under treatment $k$ as a function of the 
conditional distribution under control and a scalar parameter
$\delta_k$:
\begin{eqnarray*}
F(y \mid \rS = b, \rT = k) = g(F(y \mid \rS = b, \rT = 1), \delta_k).
\end{eqnarray*}
The model is assumed to hold for all blocks $b = 1,
\dots, B$, treatments $k = 2, \dots, K$, and outcome values $y \in \samY$ based on parameters
$\delta_2, \dots, \delta_K \in \R$. For notational convenience, we define $\delta_1 := 0$. 

This model formulation gives rise to several specific models, for
example, $g_\text{L}(p, \delta) = p^{\exp(-\delta)}$ (Lehmann alternatives),
$g_\text{PH}(p, \delta) = 1 - (1 -
p)^{\exp(-\delta)}$ (proportional hazards),
$g_\text{PO}(p, \delta) = \text{expit}(\text{logit}(p) - \delta)$ (proportional
odds), or $g_\text{Cd}(p, \delta) =
\Phi(\Phi^{-1}(p) - \delta)$ (generalised Cohen's $d$).

Instead of directly working with $g$, we parameterise the model in terms of
some absolute continuous cdf $F$ with log-concave density $f = F^\prime$
and corresponding derivative $f^\prime$. The location model 
\begin{eqnarray} \label{model}
F_Y(y \mid \rS = b, \rT = k) = F\left(F^{-1}\left(F_Y(y \mid \rS = b, \rT = 1)\right) - \delta_k\right), \quad k = 2, \dots, K
\end{eqnarray}
describes different distributions by means of shift parameter on a latent
scale defined by $F$. The negative shift term ensures that positive values of $\delta_k$ correspond
to the situation of outcomes being stochastically larger in group $k$
compared to control.
The shift parameters are invariant with respect to monotone transformations
of the response values, that is, transforming the observations of all
treatment groups by the same function does not affect the values of
$\delta_k$.

The choice $F(z) = \exp(-\exp(-z))$ gives rise to $g_\text{L}$, 
$F(z) = 1 - \exp(-\exp(z))$ corresponds to $g_\text{PH}$, $F = \text{expit}$
leads to  $g_\text{PO}$, and $F = \Phi$ results in $g_\text{Cd}$. The choice
of $F$ is made a priori and determines the interpretation of $\delta_k$. 

This document describes the implementation of estimators of these shift parameters,
as well as of confidence intervals and formal hypothesis tests for contrasts thereof under
the permutation and population model. Proportional odds models ($g_\text{PO}$)
are explained in-depths by \cite{Harrell2015RMS}, although the models are
presented in terms of survivor, not distribution, functions. 

\paragraph{Hypothesis}

We are interested in inference for $\delta_2, \dots, \delta_K$, in terms of
confidence intervals and hypothesis tests of the form
\begin{eqnarray*}
& & H_0: \delta_k - \mu_k = 0, \text{``two.sided''}, \quad k = 2, \dots, K, \\
& & H_0: \delta_k - \mu_k \ge 0, \text{``less''}, \quad k = K = 2, \\
& & H_0: \delta_k - \mu_k \le 0, \text{``greater''}, \quad k = K = 2,
\end{eqnarray*}
with the latter two options only for the two-sample case ($K = 2$).

\paragraph{Likelihood}

For an ordered categorical outcome $Y$ from sample space $\samY = \{y_1 < y_2 < \cdots <
y_C\}$, we parameterise the model in terms of intercept ($\vartheta_\cdot$) and
shift ($\delta_\cdot$) parameters
\begin{eqnarray*}
F_Y(y_c \mid \rS = b, \rT = k) = F(\vartheta_{c,b} - \delta_k), \quad c = 1, \dots,
C,
\end{eqnarray*}
that is we replace the transformed control outcome $F^{-1}\left(F_Y(y_c \mid
\rS = b, \rT = 1)\right) =
\vartheta_{c,b}$ with a corresponding intercept parameter.
These $C - 1$ intercept parameters are block-specific and monotone increasing
$\vartheta_{0,b} = -\infty < \vartheta_{1,b} < \cdots < \vartheta_{C,b} = \infty$
within each block $b = 1, \dots, B$.

We collect all model parameters in a vector
\begin{eqnarray*}
\thetavec = (\theta_1 & := & \delta_2, \\
               & \dots & , \\
               \theta_{K - 1} & := & \delta_K, \\
               \theta_{K} & := & \vartheta_{1,1}, \\
               \theta_{K + 1} & := & \vartheta_{2,1} > \vartheta_{1,1}, \\
               &  \dots, & \\
               \theta_{K + C - 2} & := & \vartheta_{C-1,1} > \vartheta_{C-2,1}, \\
               \theta_{K + C - 1} & := & \vartheta_{1,2}, \\
               & \dots &, \\
               \theta_{B (C - 1) + K - 1} & := & \vartheta_{C-1,B} >
               \vartheta_{C-2,B}).
\end{eqnarray*}
If there is no observation for level $c$ in block $b$, the corresponding
parameter is not identified and removed from $\thetavec$.
The parameter space is defined by all parameter vectors $\thetavec$
satisfying the monotonicity of the intercept parameters. Violations always
lead to the log-likelihood function being undefined and thus taking the
value $-\infty$. \cite{Harrell2024} evaluates unconstrained optimisation in
this context and recommends Newton-based algorithms leveraging the
analytically available Hessian (see below).

For the $i$th observation $(y_i = y_c, s_i = b, \rt_i = k)$ from block $b$
under treatment $k$, the log-likelihood contribution is
\begin{eqnarray*}
\log(\Prob(y_{c - 1} < Y \le y_c \mid \rS = b, \rT = k)) = \log(F(\vartheta_{c,b} - \delta_k) - F(\vartheta_{c - 1,b} - \delta_k)).
\end{eqnarray*}

For an absolutely continuous outcome $Y \in \R$, we define $y_c := y_{(c)}$,
the $c$th distinct ordered observation in the sample. The log-likelihood
above is then the empirical or nonparametric log-likelihood.

If observations were independently right-censored, the contribution of the
event $Y > \tilde{y}$ to the log-likelihood is
\begin{eqnarray*}
\log(\Prob(Y > \tilde{y} \mid \rS = b, \rT = k)) = \log(1 - F(\vartheta_{c - 1,b} - \delta_k))
\end{eqnarray*}
where $y_{c - 1} = \max \{y \in \samY \mid y \le \tilde{y}\}$, that is,
observations right-censored between $y_{c - 1}$ and $y_c$ correspond to the
parameter $\vartheta_{c - 1,b}$.

Maximising this form of the log-likelihood leads to semiparametrically efficient
estimators \citep[Chapter 15.5][]{vdVaart1998}. In this framework, tests
against deviations from the hypothesis $H_0$ above are locally most
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
Nsim <- 100
options(digits = 5)
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
@<NewtonRaphson@>
@@

We start implementing the log-likelihood function for parameters \code{parm}
$= \thetavec$ (assuming only a single block) with data from a two-way $C
\times K$ contingency table \code{x}. 

From $\thetavec$, we first extract the shift parameters $\delta_\cdot$ and
then the intercept parameters $\vartheta_\cdot$ and evaluate the probabilities
\code{prb} $ = \Prob(y_{c - 1} < Y \le y_c \mid \rS = 1, \rT = k)$ for all
groups:

@d parm to prob
@{
bidx <- seq_len(ncol(x) - 1L)
delta <- c(0, mu + parm[bidx])
intercepts <- c(-Inf, parm[- bidx], Inf)
tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                  ncol = ncol(x),
                                  byrow = TRUE)
Ftmb <- F(tmb)
if (rightcensored) {
    prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
} else {
    prb <- Ftmb[- 1L, , drop = FALSE] - 
           Ftmb[- nrow(Ftmb), , drop = FALSE]
} 
@}

If the table \code{x} represents right-censored observations, we compute
\code{prb} $ = 1 - \Prob(Y \le y_c \mid \rS = 1, \rT = k)$.

With default null values $\mu_k = 0, k = 2, \dots, K$, we define the
negative log-likelihood function as the weighted (by number of observations) sum of
the log-probabilities

@d negative logLik
@{
.nll <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<parm to prob@>
    if (any(prb < .Machine$double.eps^10)) 
        return(Inf)
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
parameters $\vartheta_\cdot$ and $\delta_\cdot$ is given in many places
\citep[for example in][Formula~(2)]{HothornMoestBuehlmann2017}. 
We begin computing the ratio of $f(\vartheta_{c,1} -
\delta_k)$ and the corresponding likelihood

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
    ret[bidx] <- .colSums(zl, m = nrow(zl), n = ncol(zl))[-1L] -
                 .colSums(zu[-nrow(zu),,drop = FALSE], 
                          m = nrow(zu) - 1L, n = ncol(zu))[-1L]
    ret[- bidx] <- .rowSums(zu[-nrow(zu),,drop = FALSE] - 
                            zl[-1,,drop = FALSE], 
                            m = nrow(zu) - 1L, n = ncol(zu))
    return(- ret)
}
@}

In addition, we define negative score residuals, that is, the derivative of the
negative log-likelihood with respect to an intercept term constrained to
zero:

@d negative score residuals
@{
.nsr <- function(parm, x, mu = 0, rightcensored = FALSE) {
    @<parm to prob@>
    @<density prob ratio@>

    ret <- .rowSums(zl - zu, m = nrow(zl), n = ncol(zl)) / 
           .rowSums(x, m = nrow(x), n = ncol(x))
    ret[!is.finite(ret)] <- 0
    return(- ret)
}
@}

We also need access to the observed Fisher information of the shift
parameters. We proceed by implementing the Hessian for the intercept
($\vartheta_\cdot$) and shift ($\delta_\cdot$) parameters, as given in Formula~(4) of
\cite{HothornMoestBuehlmann2017} first. This partitioned matrix
\begin{eqnarray*}
\mH(\vartheta_1, \dots, \vartheta_{C - 1}, \delta_2, \dots, \delta_K) = 
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
computation of $f(\vartheta_{c,1} - \delta_k)$ and $f^\prime(\vartheta_{c,1} -
\delta_k)$. We start preparing these objects, keeping in mind to remove terms
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
Aoffdiag <- - .rowSums(x * du * dl / prb^2, m = nrow(x), n = ncol(x))[-i2]
Aoffdiag <- Aoffdiag[-length(Aoffdiag)]
@}

and the diagonal elements of $\mA$ as
@d diagonal elements for Hessian of intercepts
@{
Adiag <- - .rowSums((x * dpu / prb)[-i1,,drop = FALSE] - 
                    (x * dpl / prb)[-i2,,drop = FALSE] - 
                    ((x * du^2 / prb^2)[-i1,,drop = FALSE] + 
                     (x * dl^2 / prb^2)[-i2,,drop = FALSE] ), 
                    m = nrow(x) - length(i1), n = ncol(x)
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

Z <- - .colSums(xm1 * (dpum1 / prbm1 - 
                       dplm1 / prbm1 -
                       (dum1^2 / prbm1^2 - 
                        2 * dum1 * dlm1 / prbm1^2 +
                        dlm1^2 / prbm1^2
                       )
                      ),
                m = nrow(xm1), n = ncol(xm1)
                )
if (length(Z) > 1L) Z <- diag(Z)
@}

We return the three matrices $\mA$, $\X$, and $\Z$ necessary for two
different purposes: We need the \code{full} Hessian for all parameters $\thetavec$ as a
dense \code{matrix} such that \code{nlminb} can compute updates from this
object. In addition, the
computation of the Fisher information for $\delta_2, \dots, \delta_K$ as the Schur
complement $\Z - \X^\top \mA^{-1} \X$. Because the matrix $\mA$ is symmetric
tridiagonal, we use infrastructure from the \pkg{Matrix} package to
represent this matrix:

@d Hessian
@{
.hes <- function(parm, x, mu = 0, rightcensored = FALSE, full = FALSE) {
    @<parm to prob@>

    @<Hessian prep@>

    @<off-diagonal elements for Hessian of intercepts@>
    @<diagonal elements for Hessian of intercepts@>
    @<intercept / shift contributions to Hessian@>

    if (length(Adiag) > 1L) {
        if (!isFALSE(full)) {
            A <- list(Adiag = Adiag, Aoffdiag = Aoffdiag)
        } else {
            A <- Matrix::bandSparse(length(Adiag), 
                k = 0:1, diagonals = list(Adiag, Aoffdiag), 
                symmetric = TRUE)
        }
    } else {
        if (!isFALSE(full)) {
            A <- list(Adiag = Adiag, Aoffdiag = NULL)
        } else {
            A <- matrix(Adiag)
        }
    }
    return(list(A = A, X = X, Z = Z))
}
@}

We start with an example involving $K = 3$ groups for a binary outcome and
use a binary logistic regression model to estimate the two log-odds ratios
$\delta_2$ and $\delta_3$ along with their estimated covariance
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
logistic. We use \code{optim} with numerically approximated Hessian to be
able to check the correctness of the analytical Hessian. 
Note that \code{glm} operates with a positive linear predictor, so
we need to change the sign of the log-odds ratios:

<<glm-op>>=
F <- plogis
f <- dlogis
op <- optim(par = c("mt2" = 0, "mt3" = 0, "(Intercept)" = 0), 
            fn = .nll, gr = .nsc, 
            x = x, method = "BFGS", hessian = TRUE)
cbind(glm = c(cf[-1] * -1, cf[1]), free1way = op$par)
logLik(m)
-op$value
@@

Parameter estimates and the in-sample log-likelihood are practically
identical. We now turn to the inverse Hessian of the shift terms, first
defining the derivative of the density of the standard logistic distribution
<<glm-H>>=
fp <- function(x) {
    p <- plogis(x)
    p * (1 - p)^2 - p^2 * (1 - p)
}
H <- .hes(op$par, x)
### analytical covariance of parameters
solve(H$Z - crossprod(H$X, Matrix::solve(H$A, H$X)))
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
tables (yet still allowing zero row sums):

@d table2list body
@{
dx <- dim(x)
if (length(dx) == 1L)
    stop("")
if (length(dx) == 2L)
    x <- as.table(array(x, dim = c(dx, 1)))
dx <- dim(x)
if (length(dx) < 3L)
    stop("Incorrect dimensions")
C <- dim(x)[1L]
K <- dim(x)[2L]
B <- dim(x)[3L]
if (C < 2L)
    stop("At least two response categories required")
if (K < 2L)
    stop("At least two groups required")
xrc <- NULL
if (length(dx) == 4L) {
    if (dx[4] == 2L) {
        xrc <- array(x[,,,"FALSE", drop = TRUE], dim = dx[1:3])
        x <- array(x[,,,"TRUE", drop = TRUE], dim = dx[1:3])
    } else {
        stop(gettextf("%s currently only allows independent right-censoring",
                              "free1way"),
                     domain = NA)
    }
}

xlist <- xrclist <- vector(mode = "list", length = B)

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
strata <- !vapply(xlist, is.null, NA)
xlist <- xlist[strata]
xrclist <- xrclist[strata]
@}

@d table2list
@{
.table2list <- function(x) {

    @<table2list body@>

    ret <- list(xlist = xlist)
    if (!is.null(xrc))
        ret$xrclist <- xrclist
    ret$strata <- strata
    ret
}
@}

We first extract the shift parameters $\delta_{\cdot}$ and then, separately
for each stratum, the corresponding contrasts of the intercept parameters:

@d stratum prep
@{
C <- vapply(x, NROW, 0L) ### might differ by stratum
K <- unique(do.call("c", lapply(x, ncol))) ### the same
B <- length(x)
sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), levels = seq_len(B))
bidx <- seq_len(K - 1L)
delta <- parm[bidx]
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
        ret <- ret + .nll(c(delta, intercepts[[b]]), x[[b]], mu = mu,
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
        nsc <- .nsc(c(delta, intercepts[[b]]), x[[b]], mu = mu,
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
        sr[idx] <- .nsr(c(delta, intercepts[[b]]), x[[b]], mu = mu,
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
cbind(glm = c(cf[-(1:2)] * -1, cf[1:2]), free1way = op$par)
logLik(m)
-op$value
@@

For the analytical Hessian, we sum-up over the stratum-specific
Hessians of the shift parameters. For right-censored observations, we need
to compute the contributions by the events and obtain the joint Hessian for
shift- and intercept parameters first. We differentiate between computing
the null Hessian for $\thetavec$ as a dense \code{matrix}:

@d full Hessian
@{
for (b in seq_len(B)) {
    H <- .hes(c(delta, intercepts[[b]]), x[[b]], mu = mu, full = full)
    if (!is.null(xrc)) {
        Hrc <- .hes(c(delta, intercepts[[b]]), xrc[[b]], mu = mu, 
                    rightcensored = TRUE, full = full)
        H$X <- H$X + Hrc$X
        H$A$Adiag <- H$A$Adiag + Hrc$A$Adiag
        H$A$Aoffdiag <- H$A$Aoffdiag + Hrc$A$Aoffdiag
        H$Z <- H$Z + Hrc$Z
    }
    if (b == 1L) {
        Adiag <- H$A$Adiag
        Aoffdiag <- H$A$Aoffdiag
        X <- H$X
        Z <- H$Z
    } else {
        Adiag <- c(Adiag, H$A$Adiag) ### Matrix::bdiag(A, H$A)
        Aoffdiag <- c(Aoffdiag, 0, H$A$Aoffdiag) ### Matrix::bdiag(A, H$A)
        X <- rbind(X, H$X)
        Z <- Z + H$Z
    }
}

if (length(Adiag) > 1L) {
    A <- Matrix::bandSparse(length(Adiag),
                            k = 0:1, diagonals = list(Adiag, Aoffdiag),
                            symmetric = TRUE)
} else {
    A <- matrix(Adiag)
}

ret <- cbind(Z, t(X))
ret <- rbind(ret, cbind(X, A))
if (retMatrix) return(ret)
return(as.matrix(ret))
@}

and the computation of the Hessian for the shift parameters using
\code{Matrix} technology:

@d stratified Hessian
@{
.shes <- function(parm, x, mu = 0, xrc = NULL, full = FALSE, retMatrix = FALSE) {
    @<stratum prep@>
    if (!isFALSE(ret <- full)) {
        @<full Hessian@>
    }
    ret <- matrix(0, nrow = length(bidx), ncol = length(bidx))
    for (b in seq_len(B)) {
        H <- .hes(c(delta, intercepts[[b]]), x[[b]], mu = mu)
        if (!is.null(xrc)) {
            Hrc <- .hes(c(delta, intercepts[[b]]), xrc[[b]], mu = mu, 
                        rightcensored = TRUE)
            H$X <- H$X + Hrc$X
            H$A <- H$A + Hrc$A
            H$Z <- H$Z + Hrc$Z
        }
        sAH <- tryCatch(Matrix::solve(H$A, H$X), error = function(e) NULL)
        if (is.null(sAH))
            stop(gettextf("Error computing the Hessian in %s",
                          "free1way"),
                     domain = NA)
        ret <- ret + (H$Z - crossprod(H$X, sAH))
    }
    as.matrix(ret)
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
               ret <- vapply(p, function(p) 
                   uniroot(f, PI = p, interval = 50 * c(-1, 1))$root, 0)
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

\chapter{Optimisation}

\cite{Harrell2024} reports on experiments with a number of optimisers for the specific
optimisation problem arising here and recommends a Newton-Raphson algorithm
leveraging the sparse matrix structure of the observed Fisher information
matrix. The following code was adopted from his \pkg{rms} package, function
\code{rms:::lrm.fit}.

@d Newton update
@{
gradthe <- gradient(theta)     # Compute the gradient vector
hessthe <- hessian(theta)      # Compute the Hessian matrix

delta <- Matrix::solve(hessthe, gradthe, tol = control$tolsolve)

if (control$trace)
    cat(iter, ': ', theta, "\n", sep = "")

step_size <- 1L                # Initialize step size for step-halving
@}

@d Newton step halving
@{
new_theta <- theta - step_size * delta # Update parameter vector
objnew_the <- objective(new_theta)

if (control$trace)
    cat("Old, new, old - new objective:", 
        objthe, objnew_the, objthe - objnew_the, "\n")

# Objective function failed to be reduced or is infinite
if (!is.finite(objnew_the) || (objnew_the > objthe + 1e-6)) {
    step_size <- step_size / 2         # Reduce the step size

    if (control$trace) 
        cat("Step size reduced to", step_size, "\n")

    if (step_size <= control$minstepsize) {
        msg <- paste("Step size ", step_size, " has reduced below minstepsize")
        return(list(par = theta, objective = objthe, convergence = 1, message = msg)) 
    }
} else {
    theta  <- new_theta                   # Accept the new parameter vector
    oldobj <- objthe
    objthe <- objnew_the
    break
}
@}

@d Newton convergence
@{
# Convergence check - must meet 3 criteria
if ((objthe <= oldobj + 1e-6 && (oldobj - objthe < control$objtol)) &&
    (max(abs(gradthe)) < control$gradtol) &&
    (max(abs(delta)) < control$paramtol))

    return(list(par            = theta,
                objective      = objthe,
                convergence    = 0,
                message        = "Normal convergence"))
@}

@d NewtonRaphson
@{
### adopted from rms:::lrm.fit
.NewtonRaphson <- function(start, objective, gradient, hessian, 
                           control = list(iter.max = 150, trace = trace, 
                                          objtol = 5e-4, gradtol = 1e-5, 
                                          paramtol = 1e-5, minstepsize = 1e-2, 
                                          tolsolve = .Machine$double.eps),
                           trace = FALSE
                           )
{

    theta  <- start # Initialize the parameter vector
    oldobj <- Inf
    objthe <- objective(theta)
    if (!is.finite(objthe)) {
        msg <- "Infeasible starting values"
        return(list(par = theta, objective = objthe, convergence = 1, message = msg)) 
    }

    ### Note: This is done in the call to .free1wayML
    ### gradtol <- gradtol * n / 1000.

    for (iter in seq_len(control$iter.max)) {

        @<Newton update@>

        # Step-halving loop
        while (TRUE) {
            @<Newton step halving@>
        }

        @<Newton convergence@> 
    }

    msg <- paste("Reached", control$iter.max, "iterations without convergence")
    return(list(par = theta, objective = objthe, convergence = 1, message = msg)) 
}
@}

<<Newton, echo = FALSE>>=
@<Newton@>
@@

We can now test the optimiser on a least-squares problem

<<Newton-test>>=
N <- 10000
P <- 30
X <- matrix(rnorm(N * P), ncol = P)
y <- X %*% runif(P) + rnorm(nrow(X))

f <- function(par) sum((y - X %*% par)^2)
g <- function(par) colSums(- 2 * c(y - X %*% par) * X)
h <- function(par) 2 * crossprod(X)

start <- runif(P)

cf <- .NewtonRaphson(start = start, objective = f, gradient = g, hessian = h)

cf2 <- coef(m <- lm(y ~ 0 + X))
all.equal(sum((y - fitted(m))^2), cf$objective)
all.equal(unname(cf$par), unname(cf2))
@@

\chapter{ML Estimation}
\label{ch:ML}

@o free1way.R -cp
@{
@<NewtonRaphson@>
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
@<rfree1way@>
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
    bC <- nrow(xlist[[b]]) - 1L
    lwr <- c(lwr, -Inf, rep.int(0, times = bC - 1L))
    if (NS) {
        ecdf0 <- cumsum(rowSums(xlist[[b]]))
        ### ensure that 0 < ecdf0 < 1 such that quantiles exist
        ecdf0 <- pmax(1, ecdf0[-length(ecdf0)]) / (max(ecdf0) + 1)
        Qecdf <- Q(ecdf0)
        start <- c(start, Qecdf)
    }
}
@}

The profile negative log-likelihood can be evaluated for some of the
parameters in $\thetavec$ (denoted as \code{fix}), the remaining parameters
are updated. Note that \code{start} must contain the full and feasible
(meeting monotonicity constraints for the intercept parameters) parameter vector
$\thetavec$.

We call \code{nlminb} and will increase the \code{eval.max} and
\code{iter.max} control parameters if
we encounter optimisation issues and restart at the current solution:

@d do optim
@{
maxit <- control[[1L]]$iter.max
while(maxit < 10001) {
   ret <- do.call(names(control)[[1L]], opargs)
   maxit <- 5 * maxit
   if (ret$convergence > 0) {
       opargs$control$eval.max <- maxit
       opargs$control$iter.max <- maxit
       opargs$start <- ret$par
   } else {
       break()
   }
}

if (isTRUE(correctFirth)) {
    @<Firth correction@>
} else {
    if (ret$convergence > 0) {
        if (is.na(correctFirth)) { ### only after failure
            warning(gettextf(paste("Firth correction was applied in %s because initial optimisation failed with:", 
                             ret$message),
                            "free1way"),
                             domain = NA)
            correctFirth <- TRUE
            @<Firth correction@>
        }
   }
}
if (ret$convergence > 0)
    warning(gettextf(paste("Unsuccessful optimisation in %s:", ret$message),
                           "free1way"),
                           domain = NA)

ret$correctFirth <- correctFirth
ret$value <- ret$objective
ret$objective <- NULL
@}

We first set-up the target function (the negative log-likelihood, also
dealing with right-censoring) and the corresponding gradient. We then add
the profile negative log-likelihood, which in turn calls the two functions
defined first

@d profile
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

### allocate memory for hessian
Hess <- Matrix::Matrix(0, nrow = length(start), ncol = length(start))

he <- function(par) {
    if (!is.null(xrc)) {
        ret <- .shes(par, x = xlist, mu = mu, xrc = xrclist, full = Hess, 
                     retMatrix = names(control)[1L] == ".NewtonRaphson")
    } else {
        ret <- .shes(par, x = xlist, mu = mu, full = Hess, 
                     retMatrix = names(control)[1L] == ".NewtonRaphson")
    }
    return(ret)
}
.profile <- function(start, fix = seq_len(K - 1)) {
    if (!all(fix %in% seq_len(K - 1)))
        stop(gettextf("Incorrect argument 'fix' in %s",
                      "free1way"),
                     domain = NA)
    delta <- start[fix]
    opargs <- list(start = start[-fix], 
                     objective = function(par) {
                         p <- numeric(length(par) + length(fix))
                         p[fix] <- delta
                         p[-fix] <- par
                         fn(p)
                     },
                     gradient = function(par) {
                         p <- numeric(length(par) + length(fix))
                         p[fix] <- delta
                         p[-fix] <- par
                         gr(p)[-fix]
                     },
                     hessian = function(par) {
                         p <- numeric(length(par) + length(fix))
                         p[fix] <- delta
                         p[-fix] <- par
                         he(p)[-fix, -fix, drop = FALSE]
                     })
    opargs$control <- control[[1L]]
    correctFirth <- FALSE ### turn off Firth correction in .profile
    @<do optim@>
    p <- numeric(length(start))
    p[fix] <- delta
    p[-fix] <- ret$par
    ret$par <- p
    ret
}
@}

Chapter~\ref{ch:Firth} introduces a bias correction \citep{Firth1993},
essentially by adding a penalty term to the log-likelihood. The
\code{correctFirth} argument triggers this bias correction to by applied
when \code{TRUE}, not to be applied when \code{FALSE}, and to applied in
case the unpenalised ML estimation resulted in a convergance issue
(\code{NA}). This part is still experimental and needs more testing. It also
seems unclear if and how the Fisher information needs additional correction
and it is certainly unclear if one can proceed with permutation testing
after correcting the bias.

@d Firth correction
@{
.Firth_ll <- function(cf, start) {
    fix <- seq_along(cf)
    start[fix] <- cf
    ### compute profile likelihood w/o warnings
    ret <- suppressWarnings(.profile(start, fix = fix))
    Hfull <- he(ret$par)
    Hfix <- as.matrix(solve(solve(Hfull)[fix, fix]))
    ret$value - .5 * determinant(Hfix, logarithm = TRUE)$modulus
}
if (K == 2) {
    MLcf <- ret$par[seq_len(K - 1)]
    Fret <- optim(MLcf, fn = .Firth_ll, start = ret$par,
                  method = "Brent", lower = MLcf - 5, upper = MLcf + 5)
} else {
    ### Nelder-Mead
    Fret <- optim(ret$par[seq_len(K - 1)], fn = .Firth_ll, start = ret$par)
}
if (Fret$convergence == 0) {
    start <- ret$par
    start[seq_len(K - 1)] <- Fret$par
    ret <- .profile(start, fix = seq_len(K - 1))
    ret$objective <- ret$value
}
@}

The heart of the function is a call to \code{nlminb}, trying to obtain
parameter estimates of $\thetavec$ by minimising the negative
log-likelihood. We allow some (or all) parameters to be fixed at some
constants, and provide a profile version of the likelihood:

@d optim
@{
if (!length(fix)) {
    opargs <- list(start = start, 
                   objective = fn, 
                   gradient = gr,
                   hessian = he)
    opargs$control <- control[[1L]]
    @<do optim@>
} else if (length(fix) == length(start)) {
    ret <- list(par = start, 
                value = fn(start))
} else {
    ret <- .profile(start, fix = fix)
}
@}

After parameter estimation, we evaluate negative scores, the Hessian, and
negative residuals as requested:

@d post processing
@{
if (is.null(fix) || (length(fix) == length(start)))
    parm <- seq_len(K - 1)
else 
    parm <- fix
if (any(parm >= K)) return(ret)

ret$coefficients <- ret$par[parm]
dn2 <- dimnames(xt)[2L]
names(ret$coefficients) <- cnames <- paste0(names(dn2), dn2[[1L]][1L + parm])

par <- ret$par

if (score) {
    ret$negscore <- .snsc(par, x = xlist, mu = mu)[parm]
    if (!is.null(xrc))
        ret$negscore <- ret$negscore + .snsc(par, x = xrclist, mu = mu, 
                                             rightcensored = TRUE)[parm]
}
if (hessian) {
    if (!is.null(xrc)) {
        ret$hessian <- .shes(par, x = xlist, mu = mu, xrc = xrclist)
    } else {
        ret$hessian <- .shes(par, x = xlist, mu = mu)
    }
    ret$vcov <- solve(ret$hessian)
    if (length(parm) != nrow(ret$hessian))
       ret$hessian <- solve(ret$vcov <- ret$vcov[parm, parm, drop = FALSE])
    rownames(ret$vcov) <- colnames(ret$vcov) <- rownames(ret$hessian) <-
        colnames(ret$hessian) <-  cnames
}
if (residuals) {
    ret$negresiduals <- .snsr(par, x = xlist, mu = mu)
    if (!is.null(xrc)) {
        rcr <- .snsr(par, x = xrclist, mu = mu, rightcensored = TRUE)
        ret$negresiduals <- c(rbind(matrix(ret$negresiduals, nrow = C),
                                    matrix(rcr, nrow = C)))
     }
}
ret$profile <- function(start, fix)
    .free1wayML(xt, link = link, mu = mu, start = start, fix = fix, tol = tol, 
               ...) 
ret$table <- xt

ret$strata <- strata
ret$mu <- mu
if (length(ret$mu) == 1) {
    names(ret$mu) <- link$parm
} else {
    names(ret$mu) <- c(paste(link$parm, cnames[1L], sep = ":"), cnames[-1L])
}
@}

Finally, we put everything into one function which returns an object of
class \code{free1wayML} for later use. The control parameters for
\code{.NewtonRaphson} and \code{stats::nlminb} are the ones suggested by
\cite{Harrell2024}. By default, the internal Newton-Raphson implementation
is used, we can switch to \code{stats::nlminb} by specifying \code{dooptim =
"nlminb"}. The latter option cannot handle Fisher information matrices in
form of a \code{Matrix} object and thus computing the updates takes more
time whenever a larger number of intercept parameters in present in the
problem.

@d ML estimation
@{
.free1wayML <- function(x, link, mu = 0, start = NULL, fix = NULL, 
                        residuals = TRUE, score = TRUE, hessian = TRUE, 
                        correctFirth = FALSE,
                        ### use nlminb for small sample sizes
                        dooptim = c(".NewtonRaphson", "nlminb")[1 + (sum(x) < 20)],                         
                        control = list(
                            "nlminb" = list(trace = trace, iter.max = 200,
                                            eval.max = 200, rel.tol = 1e-10,
                                            abs.tol = 1e-20, xf.tol = 1e-16),
                            ".NewtonRaphson" = list(iter.max = 200, trace = trace, 
                                             objtol = 5e-4, 
                                             gradtol = 1e-5 * sum(x) / 1000, 
                                             paramtol = 1e-5, minstepsize = 1e-2, 
                                             tolsolve = .Machine$double.eps)
                        )[dooptim],
                        trace = FALSE, 
                        tol = sqrt(.Machine$double.eps), ...) {

    ### convert to three-way table
    xt <- x
    if (!is.table(x))
      stop(gettextf("Incorrect argument 'x' in %s",
                    "free1way"),
                     domain = NA)
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
### new2old parameterisation
c(cf[1:2], cf[3], log(cf[4] - cf[3]), cf[5])
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
    @<Permutation statistics@>
}
@}

\section{Wald Statistics}

We only need access to the parameter estimates $\hat{\delta}_2, \dots,
\hat{\delta}_K$ and the corresponding Hessian:

@d Wald statistic
@{
if (alternative == "two.sided") {
    STATISTIC <- c("Wald chi-squared" = c(crossprod(cf, x$hessian %*% cf)))
    DF <- c("df" = length(parm))
    PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
} else {
    STATISTIC <- c("Wald Z" = unname(c(cf * sqrt(c(x$hessian)))))
    PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
}
@}

\section{Likelihood-ratio Statistics}

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

\section{Rao Score Statistics}

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
    STATISTIC <- c("Rao Z" = unname(- ret$negscore * sqrt(c(ret$vcov))))
    PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
}
@}

\chapter{Permutation Inference}
\label{ch:Perminf}

Under the permutation model, that is, in randomised experiments where the
random treatment allocation is the only relevant source of randomness, we
compute a permutation variant of the Rao score test, based on the conditional
asymptotic distribution or based on a Monte-Carlo estimate of the reference
distribution:

@d Permutation statistics
@{
par <- x$par
par[parm] <- value
ret <- x$profile(par, parm)
sc <- - ret$negscore
if (length(cf) == 1L)
    sc <- sc / sqrt(c(ret$hessian))
if (!is.null(x$exact)) {
    STATISTIC = c("W" = sc)
} else {
    Esc <- sc - x$perm$Expectation

    if (alternative == "two.sided" && length(cf) > 1L) {
        STATISTIC <- c("Perm chi-squared" = sum(Esc * solve(x$perm$Covariance, Esc)))
    } else {
        STATISTIC <- c("Perm Z" = Esc / sqrt(c(x$perm$Covariance)))
    }
}
@}

In addition, we compute permutation p-values

@d Permutation p-values
@{
if (!is.null(x$exact)) {
    PVAL <- switch(alternative,
                   "two.sided" = 2 * min(c(x$exact$ple(sc), x$exact$pgr(sc))),
                   "less" = x$exact$ple(sc),
                   "greater" = x$exact$pgr(sc))
} else {
    .pm <- function(x) sum(x) / length(x) 
    ps <- x$perm$permStat

    .GE <- function(x, y)
        (y - x) <= sqrt(.Machine$double.eps)

    .LE <- function(x, y)
        (x - y) <= sqrt(.Machine$double.eps)

    if (alternative == "two.sided" && length(cf) > 1L) {
        if (!is.null(ps)) {
            PVAL <- .pm(.GE(ps, STATISTIC))
        } else {
            DF <- c("df" = x$perm$DF)
            PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
        }
    } else {
        if (!is.null(ps)) {
            PVALle <- .pm(.LE(ps, STATISTIC))
            PVALge <- .pm(.GE(ps, STATISTIC))
            if (alternative == "two.sided")
                PVAL <- 2 * min(c(PVALle, PVALge))
            else if (alternative == "less")
                PVAL <- PVALle
            else
                PVAL <- PVALge
        } else {
            if (alternative == "two.sided")
                PVAL <- pchisq(STATISTIC^2, df = 1, lower.tail = FALSE)
            else
                PVAL <- pnorm(STATISTIC, lower.tail = alternative == "less")
        }
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
        ES <- ret$Statistic - ret$Expectation
        ret$testStat <- sum(ES * solve(ret$Covariance, ES))
    }
    ret$DF <- dim(xt)[2L] - 1L

    if (B) {
        for (j in 1:dim(xt)[3L]) {
           rt <- r2dtable(B, r = rowSums(xt[,,j]), c = colSums(xt[,,j]))
           stat <- stat + vapply(rt, function(x) .colSums(x[,-1L, drop = FALSE] * res[,j], 
                                                          m = nrow(x), n = ncol(x) - 1L), 0)
        }
        if (dim(xt)[2L] == 2L) {
             ret$permStat <- (stat - ret$Expectation) / sqrt(c(ret$Covariance))
        } else {
            ES <- matrix(stat, ncol = B) - ret$Expectation
            ret$permStat <- rowSums(crossprod(ES, solve(ret$Covariance, ES)))
        }
    }
    ret
}
@}

For the special case of the unstratified Wilcoxon two-sample test, we can
also provide exact p-values computed via the Streitberg-R\"ohmel shift
algorithm, mainly because the scores can be mapped to integers:

@d exact proportional odds
@{
.exact <- function(z, grp, w = rep.int(1, length(z))) {

    z <- rep(z, times = w)
    grp <- rep(grp, times = w)
    x <- rank(z)
    f <- 2 - all(x == floor(x))
    x <- as.integer(x * f)
    x <- x - min(x) + 1L
    sx <- sort(x)

    m <- as.integer(sum(grp > 0))
    stopifnot(m > 1)
    stopifnot(m < length(x))

    d <- .Call(stats:::C_dpermdist2, sx, m)
    s <- seq.int(from = 1L, to = sum(rev(sx)[seq_len(m)]), by = 1L)

    STATISTIC <- sum(x[grp > 0])
    F <- cumsum(d)
    S <- rev(cumsum(rev(d)))
    cf <- lm.fit(x = cbind(1, x), y = as.double(z))$coefficients

    z2x <- function(z) round((z - m * cf[1]) / cf[2])

    c(ple = function(z) sum(d[s <= z2x(z)]),    ### s and STATISTIC are integers
      pgr = function(z) sum(d[s >= z2x(z)]), 
      qle = function(q) c(m, max(s[F < q + 1e-08])) %*% cf,
      qgr = function(q) c(m, min(s[S < q + 1e-08])) %*% cf)
}
@}


As an example, consider the Wilcoxon rank sum test, where the scores under
the null are a linear function of the ranks of the data. We compute the
asymptotic and approximated reference distribution and corresponding
p-values for a test statistics in quadratic form:

<<SW>>=
set.seed(29)
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
and the exact versions are
<<Wexact>>=
wilcox_test(u ~ w, distribution = "exact")
free1way(u ~ w, exact = TRUE)
@@

<<Wexact-le>>=
wilcox_test(u ~ w, distribution = "exact", alternative = "less")
print(free1way(u ~ w, exact = TRUE), alternative = "greater")
@@

<<Wexact-gr>>=
wilcox_test(u ~ w, distribution = "exact", alternative = "greater")
print(free1way(u ~ w, exact = TRUE), alternative = "less")
@@



<TH>Ordered alternatives: Use contrast based tests in multcomp</TH>

\chapter{Distribution-free Tests in Stratified $K$-sample Oneway Layouts}

\section{\code{free1way}}

We provide a new test procedure in a generic \code{free1way}, featuring
a method for tables (the main workhorse) and additional user interfaces. 

@d link2fun
@{
if (!inherits(link, "linkfun")) {
    link <- match.arg(link)
    link <- do.call(link, list())
}
@}

<TH>\code{B = 0} comes from \code{chisq.test} and means the default
asymptotic permutation distribution. \code{B = 1000} means 1000 random
permutations. Can we use \code{B = Inf} for the exact distribution once
available?</TH>

We use the positive residuals for defining a permutation test with treatment
effect coding using the first group as control, that is, the test statistic
is defined through the sum of the positive residuals in all but the control
group. Unfortunately, most \code{stats::*.test} procedures use the second
group as control, so factors need to be releveled to obtain identical
results (this is relevant for the one-sided case).

@d free1way
@{
free1way <- function(y, ...)
    UseMethod("free1way")

free1way.table <- function(y, link = c("logit", "probit", "cloglog", "loglog"), 
                           mu = 0, B = 0, exact = FALSE, ...)
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

    if (!(length(mu) == 1L || length(mu) == d[2L] - 1L)) {
        warning("Incompatible length of argument 'mu'")
        mu <- rep(mu, length.out = d[2L] - 1L)
    }

    ret <- .free1wayML(y, link = link, mu = mu, ...)
    ret$link <- link
    ret$data.name <- DNAME
    ret$call <- cl

    @<free1way permutation tests@>

    if (ret$correctFirth) 
        ret$method <- paste(ret$method, 
            "with Firth bias correction", sep = ", ")

    class(ret) <- "free1way"
    return(ret)
}
@}

where preparations for permutations tests are performed before returning the
object

@d free1way permutation tests
@{
alias <- link$alias
if (length(link$alias) == 2L) alias <- alias[1L + (d[2] > 2L)]
stratified <- FALSE
if (length(d) == 3L) stratified <- d[3L] > 1
ret$method <- paste(ifelse(stratified, "Stratified", ""), 
                    paste0(d[2L], "-sample"), alias, 
                    "test against", link$model, "alternatives")

cf <- ret$par
### compute the permutation distribution always
### for H0: delta = 0, not delta = mu
### otherwise, permutation confidence intervals
### are not aligned with permutation p-values
cf[idx <- seq_len(d[2L] - 1L)] <- -mu
pr <- ret$profile(cf, idx)
res <- - pr$negresiduals
if (d[2L] == 2L)
    res <- res / sqrt(c(pr$hessian))

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

### exact two-sample Wilcoxon w/o stratification
if (exact) {
    if (!stratified && link$model == "proportional odds") {
        @<exact proportional odds@>
        ret$exact <- .exact(c(res, res), grp = unclass(gl(2, d[1L])) - 1L,
                            w = c(y[,,1L]))
        B <- 0
    } else {
        warning("Cannot compute exact distribution")
    }
} 
ret$perm <- .resample(res, y, B = B)

if (!is.null(names(dn))) {
    fm <- as.formula(paste(names(dn)[1:2], collapse = "~"))
    ret$terms <- terms(fm, data = as.data.frame(y))
}
@}

The \code{formula} method allows formulae \code{outcome ~ treatment |
stratum)} for model specification

@d free1way formula
@{
free1way.formula <- function(formula, data, weights, subset, na.action = na.pass, ...)
{
    cl <- match.call()

    if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")

    if (stratum <- (length(formula[[3L]]) > 1)) {
      if ((length(formula[[3L]]) != 3L) || 
          (formula[[3L]][[1L]] != as.name("|")) || 
          (length(formula[[3L]][[2L]]) !=  1L) || 
          (length(formula[[3L]][[3L]]) != 1L)) 
        stop("incorrect specification for 'formula'")
      formula[[3L]][[1L]] <- as.name("+")
    }

    formula <- terms(formula)
    if (length(attr(formula, "term.labels")) > 1L + stratum)
        stop("'formula' missing or incorrect")
    group <- attr(formula, "term.labels")[1L]

    m <- match.call(expand.dots = FALSE)
    m$formula <- formula
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
    event <- NULL
    if (inherits(y, "Surv")) {
        if (attr(y, "type") != "right")
            stop(gettextf("%s currently only allows independent right-censoring",
                          "free1way"),
                          domain = NA)
        event <- (y[,2] > 0)
        y <- y[,1]
    }
    g <- factor(mf[[group]])
    lev <- levels(g)
    DNAME <- paste(DNAME, paste0("(", paste0(lev, collapse = ", "), ")"))
    if (nlevels(g) < 2L)
        stop(gettextf("Incorrect argument 'groups' in %s, at least two groups needed",
                      "free1way"),
                      domain = NA)
    if (stratum) {
        st <- factor(mf[[3L]])
        ### nlevels(st) == 1L is explicitly allowed
        vn <- c(vn, names(mf)[3L])
        RVAL <- free1way(y = y, groups = g, blocks = st, event = event, weights = w,
                         varnames = vn, ...)
        DNAME <- paste(DNAME, paste("\n\t stratified by", names(mf)[3L]))
    } else {
        ## Call the corresponding method
        RVAL <- free1way(y = y, groups = g, event = event, weights = w, varnames = vn, ...)
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

@d variable names and checks
@{
cl <- match.call()
DNAME <- paste(varnames[1], "by", varnames[2])
groups <- factor(groups)
if (nlevels(groups) < 2L)
    stop(gettextf("Incorrect argument 'groups' in %s, at least two groups needed",
                  "free1way"),
                  domain = NA)
DNAME <- paste(DNAME, paste0("(", paste0(levels(groups), collapse = ", "), ")"))

if (!is.null(blocks)) {
    if (length(unique(blocks)) < 2L) {
        blocks <- NULL
    } else {
        blocks <- factor(blocks)
        DNAME <- paste(DNAME, "\n\t stratified by", varnames[3])
    }
}
varnames <- varnames[varnames != "NULL"]
@}

Note that the return value of \code{unique} might differ between platforms.
Because users can decide about the unique values in the vector \code{y} (by using
\code{round} or \code{trunc}, for example), before calling this function, we
refrain from handling this issue internally. However, we offer an
\code{nbins} argument for binning response observations at sample quantiles 
in the absence of right-censoring.

@d free1way numeric
@{
free1way.numeric <- function(y, groups, blocks = NULL, event = NULL, weights = NULL, nbins = 0, 
    varnames = c(deparse1(substitute(y)), 
                 deparse1(substitute(groups)), 
                 deparse1(substitute(blocks))), ...) {

    @<variable names and checks@>

    if (!is.null(event)) {
        if (!is.logical(event))
            stop(gettextf("%s currently only allows independent right-censoring",
                          "free1way"),
                          domain = NA)
        uy <- sort(unique(y[event]))
        if (all(y[!event] < uy[length(uy)]))
            uy <- uy[-length(uy)]
    } else {
        uy <- sort(unique(y))
    }
    if (nbins && nbins < length(uy) && is.null(event)) {
        nbins <- ceiling(nbins)
        breaks <- c(-Inf, quantile(y, probs = seq_len(nbins) / (nbins + 1L)), Inf)
    } else {
        breaks <- c(-Inf, uy, Inf)
    }
    r <- ordered(cut(y, breaks = breaks, ordered_result = TRUE, 
                     labels = FALSE)) ### avoids costly formatC call
    RVAL <- free1way(y = r, groups = groups, blocks = blocks, 
                     event = event, weights = weights, 
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
free1way.factor <- function(y, groups, blocks = NULL, event = NULL, weights = NULL, 
    varnames = c(deparse1(substitute(y)), 
                 deparse1(substitute(groups)), 
                 deparse1(substitute(blocks))), ...) {

    @<variable names and checks@>

    if (nlevels(y) > 2L && !is.ordered(y))
        stop(gettextf("%s is not defined for unordered responses",
                              "free1way"),
                     domain = NA)
    d <- data.frame(w = 1, y = y, groups = groups)
    if (!is.null(weights)) d$w <- weights
    if (is.null(blocks)) blocks <- gl(1, nrow(d))
    d$blocks <- blocks 
    if (!is.null(event)) {
       if (!is.logical(event))
            stop(gettextf("%s currently only allows independent right-censoring",
                          "free1way"),
                          domain = NA)
        d$event <- event
    }
    tab <- xtabs(w ~ ., data = d)
    dn <- dimnames(tab)
    names(dn)[seq_along(varnames)] <- varnames
    dimnames(tab) <- dn
    RVAL <- free1way(tab, ...)
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
                            tol = sqrt(.Machine$double.eps), 
                            mu = 0, ### allow permutation testing non-null hypotheses
                                    ### in alignment with confint(free1way(, B > 0))
                            ...)
{

    test <- match.arg(test)
    alternative <- match.arg(alternative)

    ### global
    cf <- coef(x)
    if ((length(cf) > 1L || test == "LRT") && alternative != "two.sided") 
        stop("Cannot compute one-sided p-values")

    DF <- NULL
    parm <- seq_along(cf)
    value <- mu

    @<statistics@>

    if (test == "Permutation") {
        @<Permutation p-values@>
    }

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
        return(.print.free1way(object, test = test, alternative = alternative, tol = tol, ...))
   
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
indices or overlap coefficients can also be computed. 

We begin computing the critical values for permutation tests, making sure
the confidence intervals will be in line with one- and two-sided p-values:
@d permutation confint
@{
if (length(cf) > 1L)
    stop("Permutation confidence intervals only available for two sample comparisons")
if (!is.null(object$exact)) {
    qu <- c(object$exact$qle(1 - conf.level),
            object$exact$qgr(1 - conf.level))
} else {
    if (is.null(object$perm$permStat)) {
        qu <- qnorm(conf.level) * c(-1, 1)
    } else {
        .pq <- function(s, alpha) {
            su <- sort(unique(s)) 
            ### F = P(T <= t), S = P(T >= t)
            Fs <- cumsum(st <- table(match(s, su)))
            Ss <- length(s) - Fs + st
            c(max(su[Fs <= alpha * length(s)]),
              min(su[Ss <= alpha * length(s)]))
        }
        ### cf PVAL computation!!!
        rs <- object$perm$permStat
        qu <- .pq(round(rs, 10), alpha = 1 - conf.level)
        att.level <- mean(rs > qu[1] & rs < qu[2])
        attr(CINT, "Attained level") <- att.level
    }
}
@}

The \code{confint} method starts with Wald intervals, which are either returned or used as starting
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
        ### we also could invert p-values, but the
        ### p-value function might be discrete for permutation
        ### tests, in contrast to the test statistic
        return(STATISTIC - quantile)
    }

    if (test == "Permutation") {
        @<permutation confint@>
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
### asymptotic permutation test
(ft <- free1way(x))
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
N <- 25
w <- gl(2, N)
y <- rlogis(length(w), location = c(0, 1)[w])

#### link = logit is default
ft <- free1way(y ~ w)

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


\section{Mantel-Haenszel Test}

<<mh>>=
example(mantelhaen.test, echo = FALSE)
mantelhaen.test(UCBAdmissions, correct = FALSE)
ft <- free1way(UCBAdmissions)
summary(ft, test = "Wald")
exp(coef(ft))
exp(confint(ft, test = "Wald"))
exp(sapply(dimnames(UCBAdmissions)[[3L]], function(dept)
       confint(free1way(UCBAdmissions[,,dept]), test = "Permutation")))
sapply(dimnames(UCBAdmissions)[[3L]], function(dept)
       fisher.test(UCBAdmissions[,,dept], conf.int = TRUE)$conf.int)
@@

\section{\code{prop.test}}

<<pt>>=
prop.test(UCBAdmissions[,,1], correct = FALSE)
summary(free1way(UCBAdmissions[,,1]), test = "Rao")
@@


\section{Kruskal-Wallis Test}

<<kw>>=
example(kruskal.test, echo = FALSE)
kruskal.test(x ~ g)
free1way(x ~ g)
@@

\section{Savage Test}

We start without censoring (Savage test) and add strata

<<sw>>=
library("survival")
N <- 10
nd <- expand.grid(g = gl(3, N), s = gl(3, N))
nd$tm <- rexp(nrow(nd))
nd$ev <- TRUE
survdiff(Surv(tm, ev) ~ g + strata(s), data = nd, rho = 0)$chisq
cm <- coxph(Surv(tm, ev) ~ g + strata(s), data = nd)

(ft <- free1way(tm ~ g | s, data = nd, link = "cloglog"))
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
                  trafo(..., numeric_trafo = logrank_trafo, block = nd$s), 
                  teststat = "quad")
@@

Wilcoxon against proportional odds

<<Peto>>=
survdiff(Surv(tm, ev) ~ g + strata(s), data = nd, rho = 1)$chisq
(ft <- free1way(tm ~ g | s, data = nd, link = "logit"))
summary(ft)
summary(ft, test = "Rao")
summary(ft, test = "LRT")
summary(ft, test = "Wald")
summary(ft, test = "Permutation")
@@

\section{Log-rank Test}

And now with censoring. We cannot expect this to be identical with what
\pkg{survival} reports, as this package is based on the partial likelihood
and we operate on the full likelihood.

<<sw>>=
library("survival")
data("GBSG2", package = "TH.data")
survdiff(Surv(time, cens) ~ horTh + strata(tgrade), data = GBSG2, rho = 0)$chisq
cm <- coxph(Surv(time, cens) ~ horTh + strata(tgrade), data = GBSG2)

### no formula interface yet
ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh | tgrade, 
                           link = "cloglog"))
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

### test with many small strata 
(ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh | pnodes, 
                            link = "cloglog")))
@@

Wilcoxon against proportional odds

<<Peto>>=
survdiff(Surv(time, cens) ~ horTh + strata(tgrade), data = GBSG2, rho = 1)$chisq
(ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh | tgrade, 
                            link = "logit")))
summary(ft)
summary(ft, test = "Rao")
summary(ft, test = "LRT")
summary(ft, test = "Wald")
summary(ft, test = "Permutation")
@@


\section{van der Waerden Test}

Normal scores test against a generalised Cohen's $d$:

<<normal>>=
nd$y <- rnorm(nrow(nd))
free1way(y ~ g | s, data = nd, link = "probit")
independence_test(y ~ g | s, data = nd, ytrafo = function(...)
                  trafo(..., numeric_trafo = normal_trafo, block = nd$s), 
                  teststat = "quad")
@@

\section{Friedman Test}

Each observation is a block

<<friedman>>=
example(friedman.test, echo = FALSE)
rt <- expand.grid(str = gl(22, 1),
                  trt = gl(3, 1, labels = c("Round Out", 
                                            "Narrow Angle", 
                                            "Wide Angle")))
rt$tm <- c(RoundingTimes)
friedman.test(RoundingTimes)
(ft <- free1way(tm ~ trt | str, data = rt))
summary(ft)
@@

\section{Contrast Tests}

\code{free1way} output can be used to define multiple contrast tests and
corresponding confidence intervals via the \pkg{multcomp} package. For
example, Tukey-style simultaneous all-pair comparisons can be implemented via
<<Tukey>>=
tk <- free1way(Ozone ~ Month, data = airquality)
library("multcomp")
confint(glht(tk, linfct = mcp(Month = "Tukey")))
@@

\chapter{Model Diagnostics}

The classical shift model $F_Y(y \mid T = 2) = F_Y(y - \mu \mid T = 1)$
can be criticised using confidence bands for QQ-plots in \code{qqplot},
because the parameter $\mu$ shows up as a vertical shift of the diagonal
if the model is appropriate.

Likewise, model~(\ref{model}) can be graphically assessed using the PP-plot.
We concentrate on the two-sample case. The shift parameter $\delta_2$ gives
rise to the model-based PP graph $(p, F(F^{-1}(p) - \delta_2))$ and a
confidence \emph{band} can be obtained from a confidence \emph{interval} for
$\delta_2$. The PP-plot is, up to rescalings, identical to the ROC curve.

@d ROC bands
@{
 if (!is.null(conf.level)) {
    prb <- seq_len(1000) / 1001
    res <- c(x, y)
    grp <- gl(2, 1, labels = c(xlab, ylab))
    grp <- grp[rep(1:2, c(length(x), length(y)))]
    args <- conf.args
    args$y <- res
    args$groups <- grp
    args$border <- args$col <- args$type <- NULL
    f1w <- do.call("free1way", args)

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
implementation of \code{qqplot}, allowing to plot the empirical 
\citep{WilkGnanadesikan1968} and
corresponding model-based \citep{SewakHothorn2023} probability-probability plot, 
the latter for a certain choice of link function:

@d ppplot
@{
ppplot <- function(x, y, plot.it = TRUE,
                   xlab = paste("Cumulative probabilities for", deparse1(substitute(x))),
                   ylab = paste("Cumulative probabilities for", deparse1(substitute(y))), 
                   main = "P-P plot",
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
    sy <- sort(unique(c(x, y)))
    py <- ecdf(y)(sy)
    px <- ex(sy)
    ret <- stepfun(px, c(0, py))
    if (!plot.it)
        return(ret)

    plot(ret, xlim = c(0, 1), ylim = c(0, 1), 
         xlab = xlab, ylab = ylab, main = main, 
         verticals = FALSE, ...)

    @<ROC bands@>

    plot(ret, add = TRUE, verticals = FALSE, ...)

    return(invisible(ret)) 
}
@}

A correct logistic model with log-odds ratio three is shown in Figure~\ref{fig:PO}
and an incorrect proportional hazards model for the same data in
Figure~\ref{fig:PH}.

\begin{figure}
<<ppplot, fig = TRUE>>=
y <- rlogis(50)
x <- rlogis(50, location = 3)
ppplot(y, x, conf.level = .95)
@@
\caption{Data sampled from a proportional-odds model with
probability-probability (PP)  curve and $95\%$ confidence band obtained from
a proportional-odds model. \label{fig:PO}}
\end{figure}



\begin{figure}
<<ppplot-savage, fig = TRUE>>=
ppplot(y, x, conf.args = list(link = "cloglog", type = "Wald", 
                              col = NA, border = NULL),
       conf.level = .95)
@@
\caption{Data sampled from a proportional-odds model with
probability-probability (PP)  curve and $95\%$ confidence band obtained from
a proportional-hazards model. \label{fig:PH}}
\end{figure}


\chapter{Random Number Generation} \label{ch:rng}

With~\ref{model} we know that
\begin{eqnarray*}
U = F_Y(Y \mid \rS = b, \rT = k) = F\left(F^{-1}\left(F_Y(Y \mid \rS = b, \rT = 1)\right) - \delta_k\right), \quad k = 2, \dots, K
\end{eqnarray*}
follows a standard uniform distribution on the unit interval. This means
that we can sample from the distribution of $Y$ using
\begin{eqnarray*}
F_Y^{-1}\left(F(F^{-1}(U) + \delta_k) \mid \rS = b, \rT = 1)\right).
\end{eqnarray*}
It is therefore enough to draw samples from $F(F^{-1}(U) + \delta_k)$, that
is, assuming a uniform distribution for $F_Y$ in each control group. Because
of the invariance with respect to monotone transformations, transforming all
observations by the same quantile function changes the outcome distributions
but not the shift effects.

@d design args
@{
K <- length(delta) + 1L
if (is.null(names(delta))) 
    names(delta) <- LETTERS[seq_len(K)[-1]]
if (length(alloc_ratio) == 1L) 
    alloc_ratio <- rep_len(alloc_ratio, K - 1)
if (length(alloc_ratio) != K - 1L)
    stop("Incorrect argument 'alloc_ratio'")
if (length(strata_ratio) == 1L) 
    strata_ratio <- rep_len(strata_ratio, B - 1)
if (length(strata_ratio) != B - 1L)
    stop("Incorrect argument 'strata_ratio'")
### sample size per group (columns) and stratum (rows)
N <- n * matrix(c(1, alloc_ratio), nrow = B, ncol = K, byrow = TRUE) * 
         matrix(c(1, strata_ratio), nrow = B, ncol = K)
@}

@d rfree1way
@{
.rfree1way <- function(n, delta = 0, link = c("logit", "probit", "cloglog", "loglog")) {

    logU <- log(ret <- runif(n))

    @<link2fun@>

    trt <- (abs(delta) > 0)
    ret[trt] <- .p(link, .q(link, logU[trt], log.p = TRUE) + delta[trt])

    return(ret)
}

rfree1way <- function(n, prob = NULL, alloc_ratio = 1, 
                      blocks = ifelse(is.null(prob), 1, NCOL(prob)), 
                      strata_ratio = 1, delta = 0, offset = 0, 
                      link = c("logit", "probit", "cloglog", "loglog"))
{
    B <- blocks

    @<design args@>

    rownames(N) <- paste0("block", seq_len(B))
    ctrl <- "Control"
    colnames(N) <- c(ctrl, names(delta))

    if (length(offset) != K)
        offset <- rep_len(offset, K)

    trt <- gl(K, 1, labels = colnames(N))
    blk <- gl(B, 1, labels = rownames(N))
    ret <- expand.grid(groups = trt, blocks = blk)
    if (B == 1L) ret$blocks <- NULL
    ret <- ret[rep(seq_len(nrow(ret)), times = N), , drop = FALSE]
    ret$y <- .rfree1way(nrow(ret), 
                        delta = offset[ret$groups] + c(0, delta)[ret$groups], 
                        link = link)
    if (is.null(prob)) return(ret)

    ### return discrete distribution
    if (!is.matrix(prob))
        prob <- matrix(prob, nrow = NROW(prob), ncol = B)
    if (ncol(prob) != B)
        stop(gettextf("Incorrect number of columns for 'prob' in %s",
                      "rfree1way"),
                      domain = NA)
    prob <- prop.table(prob, margin = 2L)
    ret <- do.call("rbind", lapply(1:ncol(prob), function(b) {
        if (B > 1)
            ret <- subset(ret, blocks == levels(blocks)[b])
        ret$y <- cut(ret$y, breaks = c(-Inf, cumsum(prob[,b])), 
                     labels = paste0("Y", 1:nrow(prob)), ordered_result = TRUE)
        ret
    }))
    return(ret)
}
@}

<<rfree1way>>=
(logOR <- c(log(1.5), log(2)))
nd <- rfree1way(150, delta = logOR)
coef(ft <- free1way(y ~ groups, data = nd))
sqrt(diag(vcov(ft)))
logLik(ft)
nd$y <- qchisq(nd$y, df = 3)
coef(ft <- free1way(y ~ groups, data = nd))
sqrt(diag(vcov(ft)))
logLik(ft)
N <- 25
pvals <- replicate(Nsim, 
{
  nd <- rfree1way(n = N, blocks = 2, delta = c(.25, .5), alloc_ratio = 2)
  summary(free1way(y ~ groups | blocks, data = nd), test = "Permutation")$p.value
})

power.free1way.test(n = N, blocks = 2, delta = c(.25, .5), alloc_ratio = 2)
mean(pvals < .05)
@@

This function can also be used to simulate survival times, for example, from
a proportional hazards model with a censoring probability of $.25$ in the
control arm and of $.5$ in the treated arm, under random censoring (that is,
event and censoring times are independent given treatment).

<<rfree1waysurv>>=
N <- 1000
nd <- rfree1way(N, delta = 1, link = "cloglog")
nd$C <- rfree1way(n = N, delta = 1, offset = -c(qlogis(.25), qlogis(.5)), 
                  link = "cloglog")$y
nd$y <- Surv(pmin(nd$y, nd$C), nd$y < nd$C)
### check censoring probability
1 - tapply(nd$y[,2], nd$groups, mean)
summary(free1way(y ~ groups, data = nd, link = "cloglog"))
summary(coxph(y ~ groups, data = nd))
@@

Next we start implementing a function for simulating $C \times K$ tables. We need
to specify the number of observations in each treatment group (\code{c}),
the discrete distribution of the control (\code{r}), a model (\code{link}), and a
treatment effect (\code{delta}, in line with \code{power.XYZ.test}). In
essence, we draw samples from the multinomial distribution after computing
the relevant discrete density.

@d r2dsim
@{
.r2dsim <- function(n, r, c, delta = 0,
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
    h0 <- .q(link, p0[-length(p0)]) ### last element of p0 is one

    h1 <- h0 - matrix(delta, nrow = length(prob) - 1L, ncol = K - 1, byrow = TRUE)
    p1 <- rbind(.p(link, h1), 1)
    p <- cbind(p0, p1)
    ret <- vector(mode = "list", length = n)

    for (i in seq_len(n)) {
        tab <- sapply(seq_len(K), function(k)
            unclass(table(cut(runif(colsums[k]), breaks = c(-Inf, p[,k])))))
        ret[[i]] <- as.table(array(unlist(tab), dim = c(length(prob), K), 
                          dimnames = list(names(prob), 
                                          names(colsums))))
    }
    return(ret)
}
@}




\chapter{Power and Sample Size}

The term ``distribution-free'' refers to the invariance of the reference
distribution with respect to the distribution of an absolutely continuous
outcome under control. Unfortunately, this is no longer true for
non-continuous outcomes (due to ties) and under the alternative. That means
that sample size assessments always take place under certain assumptions
regarding the outcome distribution.

With the infrastructure from Chapter~\ref{ch:rng}, 
we are now ready to put together a function for power evaluation and sample
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
    prob <- matrix(prob, nrow = NROW(prob), ncol = blocks)
prob <- prop.table(prob, margin = 2L)
C <- nrow(prob)
B <- ncol(prob)
if (is.null(colnames(prob))) 
    colnames(prob) <- paste0("stratum", seq_len(B))
p0 <- apply(prob, 2, cumsum)
h0 <- .q(link, p0[-nrow(p0),,drop = FALSE])

@<design args@>

rownames(N) <- colnames(prob)
ctrl <- "Control"
dn <- dimnames(prob)
if (!is.null(names(dn)[1L]))
    ctrl <- names(dn)[1L]
colnames(N) <- c(ctrl, names(delta))
@}

For estimating the Fisher information, we draw samples from the discrete
outcome distribution and evaluate the observed Fisher information for the,
here and now known true parameters. The average of these Fisher information
matrices is then used as an estimate for the expected Fisher information.
For small sample sizes less than $100$, we draw larger samples (at least
$1000$) and adjust the obtained Fisher information accordingly to reduce
sampling error.

@d estimate Fisher information
@{
he <- 0
deltamu <- delta - mu
Nboost <- ifelse(n < 100, ceiling(1000 / n), 1)
for (i in seq_len(nsim)) {
    parm <- deltamu
    x <- as.table(array(0, dim = c(C, K, B)))
    for (b in seq_len(B)) {
        x[,,b] <- .r2dsim(1L, r = prob[, b], c = Nboost * N[b,], 
                          delta = delta, link = link)[[1L]]
        rs <- which(.rowSums(x[,,b], m = dim(x)[1L], n = dim(x)[2L]) > 0)
        theta <- h0[pmin(nrow(h0), rs), b]
        parm <- c(parm, theta[-length(theta)])
    }
    ### evaluate observed hessian for true parameters parm and x data
    he <- he + .free1wayML(x, link = link, mu = mu, start = parm, 
                           fix = seq_along(parm))$hessian / Nboost
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
                    blocks = blocks,
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
if (K == 2L) ret[["Standard error"]] <- se
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
power.free1way.test <- function(n = NULL, prob = if (is.null(n)) NULL else rep.int(1 / n, n), 
                                alloc_ratio = 1, blocks = if (is.null(prob)) 1 else NCOL(prob), 
                                strata_ratio = 1, 
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

    @<r2dsim@>

    if (is.null(n)) 
        n <- ceiling(uniroot(function(n) {
                 @<power call@>
             }, interval = c(5, 1e+03), tol = tol, extendInt = "upX")$root)
    else if (is.null(delta)) {
        ### 2-sample only
        if (length(alloc_ratio) > 1L)
            stop(gettextf("Effect size can only computed for two sample problems in %s",
                          "power.free1way.test"),
                          domain = NA)
        delta <- uniroot(function(delta) {
                 @<power call@>
    ### <TH> interval depending on alternative, symmetry? </TH>
            }, interval = c(0, 10), tol = tol, extendInt = "upX")$root
        }
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) {
                @<power call@>
            }, interval = c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root

    ### n is available now
    if (is.null(prob)) prob <- rep(1 / n, n)
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
        if (alternative != "two.sided")
            stop(gettextf("%s only allows two-sided alternatives in the presence of more than two groups",
                          "power.free1way.test"),
                          domain = NA)
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

### approximate formula in Hmisc::popower
library("Hmisc")
popower(p = rep(1 / N, N), odds.ratio = exp(delta), n = 2 * N)
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

We next use the \code{rfree1way} function to sample from $4 \times 3$ tables with odds ratios $2$ and $3$
and compare the resulting power with result obtained from the approximated
Fisher information. By default, the continuous control distribution is
uniform on the unit interval, thus \code{cut} with breaks defined by the
target control discrete probability distribution generates the outcome.
The plot shows the distribution of the parameter
estimates and the corresponding population values as red dots
(Figure~\ref{fig:POsim}).	

\begin{figure}
<<table, fig = TRUE>>=
prb <- rep.int(1, 4) / 4
pw <- numeric(Nsim)
cf <- matrix(0, nrow = Nsim, ncol = length(delta))
colnames(cf) <- names(delta)
for (i in seq_along(pw)) {
    nd <- rfree1way(n = N, prob = prb, delta = delta)
    ft <- free1way(y ~ groups, data = nd)
    cf[i,] <- coef(ft)
    pw[i] <- summary(ft, test = "Permutation")$p.value
}
mean(pw < .05)
boxplot(cf)
points(c(1:2), delta, pch = 19, col = "red")
power.free1way.test(n = N, prob = prb, delta = delta)
@@
\caption{Power simulation for proportional-odds model and corresponding
power approximation. \label{fig:POsim}}
\end{figure}

In the last example, we sample from $4 \times 3$ tables with odds ratios $2$ and $3$ for three
strata with different control distributions, see Figure~\ref{fig:POstrata}, and again compare the
simulation results to the power function.

\begin{figure}
<<stable, fig = TRUE>>=
prb <- cbind(S1 = rep(1, 4), 
             S2 = c(1, 2, 1, 2), 
             S3 = 1:4)
dimnames(prb) <- list(Ctrl = paste0("i", seq_len(nrow(prb))),
                      Strata = colnames(prb))

pw <- numeric(Nsim)
cf <- matrix(0, nrow = Nsim, ncol = length(delta))
colnames(cf) <- names(delta)
for (i in seq_along(pw)) {
    nd <- rfree1way(n = N, prob = prb, delta = delta)
    ft <- free1way(y ~ groups | blocks, data = nd)
    cf[i,] <- coef(ft)
    pw[i] <- summary(ft, test = "Permutation")$p.value
}
mean(pw < .05)
boxplot(cf)
points(c(1:2), delta, pch = 19, col = "red")
@@
\caption{Power simulation for stratified proportional-odds model and corresponding
power approximation. \label{fig:POstrata}}
\end{figure}

<<powertest>>=
power.free1way.test(n = N, prob = prb, delta = delta, seed = 3)
power.free1way.test(power = .8, prob = prb, delta = delta, seed = 3)
power.free1way.test(n = 19, prob = prb, delta = delta, seed = 3)
@@

\chapter{Firth Correction} \label{ch:Firth}

Sometimes, especially under complete separation, the maximum likelihood
estimator does not exist. We could think of offering the option to add a
penalty term to the log-likelihood, for example half of the log-determinant
of the Hessian as suggested by \cite{Firth1993}. Here is an example

<<Firth>>=
N <- 20
w <- gl(2, N)
y <- rnorm(length(w), mean = c(-2, 3)[w])

x <- free1way(y ~ w, link = "probit")
coef(x)
logLik(x)

pll <- function(cf) {

    start <- x$par
    start[1] <- cf
    x$profile(start, fix = 1)
}

### https://doi.org/10.1111/j.0006-341X.2001.00114.x
### https://doi.org/10.1111/j.1467-9876.2012.01057.x
### https://doi.org/10.1186/s12874-017-0313-9
### https://files.osf.io/v1/resources/fet4d_v3/providers/osfstorage/682fb176db88f967facacb5a?format=pdf&action=download&direct&version=1
### https://doi.org/10.1002/sim.6537
fun <- function(cf) {
    ret <- pll(cf)
    ret$value - .5 * determinant(ret$hessian, logarithm = TRUE)$modulus
}

ci <- confint(x, level = .99, test = "Wald")
grd <- seq(from = ci[1], to = ci[2], length.out = 50)

optim(coef(x), fn = fun, method = "Brent", 
      lower = min(grd), upper = max(grd))[c("par", "value")]
@@

The \code{correctFirth} argument can be used to request this type of bias
correction from \code{free1way} (this argument should be added to
\code{free1way.table} and documented)
<<correctFirth>>=
free1way(y ~ w, link = "probit", correctFirth = TRUE)
@@

\chapter*{Index}

\section*{Files}

@f

\section*{Fragments}

@m

\section*{Identifiers}

@u

\bibliographystyle{plainnat}
\bibliography{\Sexpr{gsub("\\.bib", "", system.file("REFERENCES.bib", package = "free1way"))}}

\end{document}
