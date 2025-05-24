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

\usepackage[%
backref,%
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

\usepackage[round]{natbib}

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
\newcommand{\Y}{\mathbf{Y}}
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


\author{Torsten Hothorn \\ Universit\"at Z\"urich}

\title{Stratified $K$ Sample Inference}

\begin{document}

\pagenumbering{roman}
\maketitle
\tableofcontents

\chapter{Introduction}
\pagenumbering{arabic}

Treatment group $X \in \{1, \dots, K\}$, outcome $Y \in \samY$ at least ordered,
stratum $S \in \{1, \dots, B\}$ with conditional distribution function
$F_Y(y \mid S = b, X = k) = \Prob(Y \le y \mid S = b, X = k)$.

Model

$F(y \mid S = b, X = k) = g(F(y \mid S = b, X = 1), \beta_k)$ for all $b = 1,
\dots, B$, $k = 2, \dots, K$, and $y \in \samY$. $\beta_1 := 0$.

$g(p, \beta) = p^{\exp(\beta)}$ or $g(p, \beta) = 1 - (1 - p)^{\exp(\beta)}$ or
$g(p, \beta) = \text{expit}(\text{logit}(p) - \beta)$ or $g(p, \beta) =
\Phi(\Phi^{-1}(p) - \beta)$.

For some absolute continuous cdf $F$ with log-concave density $f = F^\prime$, we write
$F_Y(y \mid S = b, X = k) = F(F^{-1}(F_Y(y \mid S = b, X = 0)) - \beta_k)$.

Hypthesis
$H_0: \beta_k = \mu_k, k = 2, \dots, K$

$\samY = \{y_1 < y_2 < \cdots < y_C\}$

$F_Y(y_c \mid S = b, X = k) = F(\vartheta_{c,b} - \beta_k)$ with parameters
$\vartheta_{0,b} = -\infty < \vartheta_{1,b} < \cdots < \vartheta_{C,b} = \infty$.

For observation $(y_i = y_c, x_i = k, s_i = b)$, the log-likelihood contribution is
$\log(F(\vartheta_{c,b} - \beta_k) - F(\vartheta_{c - 1,b} - \beta_k)$

Data as $C \times K \times B$ contingency table.

$\theta = (\beta_2, \dots, \beta_K, \vartheta_{1,b}, \vartheta_{2,b} -
\vartheta_{1,b}, \dots, \vartheta_{C-1,b} - \vartheta_{C-2,b})$ with

@o strKsam.R -cp
@{
@<cumsumrev@>
@<negative logLik@>
@<negative score@>
@}

@d parm to prob
@{
b <- seq_along(ncol(x) - 1L)
beta <- c(0, mu + parm[b])
theta <- c(-Inf, cumsum(parm[- b]), Inf)
tmb <- theta - matrix(beta, nrow = length(theta),  
                            ncol = ncol(x),
                            byrow = TRUE)
Ftmb <- F(tmb)
prb <- Ftmb[, - 1L, drop = FALSE] - 
       Ftmb[, - nrow(Ftmb), drop = FALSE]
@}

@d negative logLik
@{
.nll <- function(parm, x, mu = numeric(ncol(x))) {
    @<parm to prob@>
    - sum(x * log(prb))
}
@}

@d cumsumrev
@{
.rcr1 <- function(z)
    rev(cumsum(rev(z)))
@}

@d negative score
@{
.nsc <- function(parm, x, mu = numeric(ncol(x))) {
    @<parm to prob@>
    ftmb <- f(tmb)

    ret <- numeric(length(parm))
    zu <- x * ftmb[, - 1, drop = FALSE] / prb
    zl <- x * ftmb[, - nrow(ftmb), drop = FALSE] / prb
    ret[b] <- colSums(zl)[-1L] -
              colSums(zu[-nrow(zu),,drop = FALSE])[-1L]
    ret[-b] <- Reduce("+", 
                      lapply(1:ncol(x), 
                          function(j) {
                              .rcr1(zu[-nrow(zu),j]) - 
                              .rcr1(zl[-1,j])
                          })
                     )
    -ret
}
@}


<<>>=
source("strKsam.R")
w <- matrix(c(10, 5, 11, 8), ncol = 2)
d <- expand.grid(y = gl(2, 1), x = gl(2, 1))
d$w <- c(w)
m <- glm(y ~ x, data = d, weights = w, family = binomial())
(cf <- coef(m))
logLik(m)
F <- plogis
.nll(rev(cf), w)
@@

\chapter*{Index}

\section*{Files}

@f

\section*{Fragments}

@m

\section*{Identifiers}

@u

%\bibliographystyle{plainnat}
%\bibliography{libcoin}

\end{document}
