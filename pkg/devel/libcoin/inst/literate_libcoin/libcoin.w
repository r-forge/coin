
\documentclass[a4paper]{report}

%% packages
\usepackage{amsfonts,amstext,amsmath,amssymb,amsthm}

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
pdfpagemode=None,%
pdfstartview=FitBH,%
linkcolor={linkcolor},%
anchorcolor={linkcolor},%
citecolor={linkcolor},%
filecolor={linkcolor},%
menucolor={linkcolor},%
pagecolor={linkcolor},%
urlcolor={linkcolor}%
]{hyperref}

\setlength{\oddsidemargin}{0in}
\setlength{\evensidemargin}{0in}
\setlength{\topmargin}{0in}
\addtolength{\topmargin}{-\headheight}
\addtolength{\topmargin}{-\headsep}
\setlength{\textheight}{8.9in}
\setlength{\textwidth}{6.5in}
\setlength{\marginparwidth}{0.5in}

\newcommand{\pkg}[1]{\textbf{#1}}
\newcommand{\proglang}[1]{\textsf{#1}}

\newcommand{\R}{\mathbb{R} }
\newcommand{\Prob}{\mathbb{P} }
\newcommand{\N}{\mathbb{N} }
\newcommand{\C}{\mathbb{C} }
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
\newcommand{\ws}{\mathbf{w}_\cdot}
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


\author{Torsten Hothorn \\ Universit\"at Z\"urich}

\title{The \pkg{libcoin} Package}

\begin{document}

\pagenumbering{roman}
\maketitle
\tableofcontents

\chapter{Introduction}
\pagenumbering{arabic}

In the following we assume that we are provided with $n$ observations
\begin{eqnarray*}
(\Y_i, \X_i, w_i, b_i), \quad i = 1, \dots, N.
\end{eqnarray*}
The variables $\Y$ and $\X$ from sample spaces $\mathcal{Y}$ and
$\mathcal{X}$ may
be measured at arbitrary scales and may be multivariate as well. In addition
to those measurements, case weights $w_i \in \N$ and a factor $b_i \in \{1, \dots, B\}$
coding for $B$ independent blocks may
be available.
We are interested in testing the null hypothesis of independence of $\Y$ and
$\X$
\begin{eqnarray*}
H_0: D(\Y \mid \X) = D(\Y)
\end{eqnarray*}
against arbitrary alternatives. \cite{strasserweber1999} suggest to derive
scalar test statistics for testing $H_0$ from multivariate linear statistics
of a specific linear form. Let $\A \subseteq \{1, \dots, N\}$ denote some subset of the
observation numbers and consider the linear statistic
\begin{eqnarray} \label{linstat}
\T(A) = \vec\left(\sum_{i \in \A} w_i g(\X_i) h(\Y_i, \{\Y_i \mid i \in \A\})^\top\right)
\in \R^{pq}.
\end{eqnarray}
Here, $g: \mathcal{X} \rightarrow \R^P$ is a transformation of
$\X$ known as the \emph{regression function} and $h: \mathcal{Y} \times
\mathcal{Y}^n \rightarrow \R^Q$ is a transformation of $\Y$ known as the
\emph{influence function}, where the latter may depend on $\Y_i$ for $i \in \A$
in a permutation symmetric way.  We will give specific examples on how to choose
$g$ and $h$ later on.

With $\x_i = g(\X_i) \in \R^P$ and $\y_i = h(\Y_i, \{\Y_i, i \in \A\}) \in \R^Q$ we write
\begin{eqnarray} \label{linstat}
\T(A) = \vec\left(\sum_{i \in \A} w_i \x_i \y_i^\top\right)
\in \R^{PQ}.
\end{eqnarray}
The \pkg{libcoin} package doesn't handle neither $g$ nor $h$, this is the job
of \pkg{coin} and we therefore continue with $\x_i$ and $\y_i$.

The distribution of $\T$  depends on the joint
distribution of $\Y$ and $\X$, which is unknown under almost all practical
circumstances. At least under the null hypothesis one can dispose of this
dependency by fixing $\X_i, i \in \A$ and conditioning on all possible
permutations $S(A)$ of the responses $\Y_i, i \in \A$.
This principle leads to test procedures known
as \textit{permutation tests}.
The conditional expectation $\mu(A) \in \R^{PQ}$ and covariance
$\Sigma(A) \in \R^{PQ \times PQ}$
of $\T$ under $H_0$ given
all permutations $\sigma \in S(A)$ of the responses are derived by
\cite{strasserweber1999}:
\begin{eqnarray}
\mu(A) & = & \E(\T(A) \mid S(A)) = \vec \left( \left( \sum_{i \in \A} w_i \x_i \right) \E(h \mid S(A))^\top
\right), \nonumber \\
\Sigma(A) & = & \V(\T(A) \mid S(A)) \nonumber \\
& = &
    \frac{\ws}{\ws(A) - 1}  \V(h \mid S(A)) \otimes
        \left(\sum_{i \in \A} w_i  \x_i \otimes w_i \x_i^\top \right)
\label{expectcovar}
\\
& - & \frac{1}{\ws(A) - 1}  \V(h \mid S(A))  \otimes \left(
        \sum_{i \in \A} w_i \x_i \right)
\otimes \left( \sum_{i \in \A} w_i \x_i\right)^\top
\nonumber
\end{eqnarray}
where $\ws(A) = \sum_{i \in \A} w_i$ denotes the sum of the case weights,
and $\otimes$ is the Kronecker product. The conditional expectation of the
influence function is
\begin{eqnarray*}
\E(h \mid S(A)) = \ws(A)^{-1} \sum_{i \in \A} w_i \y_i \in
\R^Q
\end{eqnarray*}
with corresponding $Q \times Q$ covariance matrix
\begin{eqnarray*}
\V(h \mid S(A)) = \ws(A)^{-1} \sum_{i \in \A} w_i \left(\y_i - \E(h \mid S(A)) \right) \left(\y_i  - \E(h \mid S(A))\right)^\top.
\end{eqnarray*}

With $A_b = \{i \mid b_i = b\}$ we get $\T = \sum_{b = 1}^B T(A_b)$,
$\mu = \sum_{b = 1}^B \mu(A_b)$ and $\Sigma = \sum_{b = 1}^B \Sigma(A_b)$.

Having the conditional expectation and covariance at hand we are able to
standardize a linear statistic $\T \in \R^{PQ}$ of the form
(\ref{linstat}). Univariate test statistics~$c$ mapping an observed linear
statistic $\mathbf{t} \in
\R^{PQ}$ into the real line can be of arbitrary form.  An obvious choice is
the maximum of the absolute values of the standardized linear statistic
\begin{eqnarray*}
c_\text{max}(\mathbf{t}, \mu, \Sigma)  = \max \left| \frac{\mathbf{t} -
\mu}{\text{diag}(\Sigma)^{1/2}} \right|
\end{eqnarray*}
utilizing the conditional expectation $\mu$ and covariance matrix
$\Sigma$. The application of a quadratic form $c_\text{quad}(\mathbf{t}, \mu,
\Sigma)  =
(\mathbf{t} - \mu) \Sigma^+ (\mathbf{t} - \mu)^\top$ is one alternative, although
computationally more expensive because the Moore-Penrose
inverse $\Sigma^+$ of $\Sigma$ is involved.

The definition of one- and two-sided $p$-values used for the computations in
the \pkg{libcoin} package is
\begin{eqnarray*}
P(c(\T, \mu, \Sigma) &\le& c(\mathbf{t}, \mu, \Sigma)) \quad \text{(less)} \\
P(c(\T, \mu, \Sigma) &\ge& c(\mathbf{t}, \mu, \Sigma)) \quad \text{(greater)}\\
P(|c(\T, \mu, \Sigma)| &\le& |c(\mathbf{t}, \mu, \Sigma)|) \quad \text{(two-sided).}
\end{eqnarray*}
Note that for quadratic forms only two-sided $p$-values are available
and that in the one-sided case maximum type test statistics are replaced by
\begin{eqnarray*}
\min \left( \frac{\mathbf{t} - \mu}{\text{diag}(\Sigma)^{1/2}} \right)
    \quad \text{(less) and }
\max \left( \frac{\mathbf{t} - \mu}{\text{diag}(\Sigma)^{1/2}} \right)
    \quad \text{(greater).}
\end{eqnarray*}

\chapter{C Code}



@S

\section{Sums}

@s

The corresponding header file contains definitions of
functions to be used outside \verb|Sums.c|

@o Sums.h -cc
@{
@<Function Prototypes@>
@}

@d Function Prototypes
@{
@<RC\_Sums Prototype@>;
@<RC\_KronSums Prototype@>;
@<RC\_KronSums\_Permutation Prototype@>;
@<RC\_colSums Prototype@>;
@<RC\_OneTableSums Prototype@>;
@<RC\_TwoTableSums Prototype@>;
@<RC\_ThreeTableSums Prototype@>;
@<RC\_LinearStatistic Prototype@>;
@<RC\_ExpectationInfluence Prototype@>;
@<RC\_CovarianceInfluence Prototype@>;
@<RC\_ExpectationX Prototype@>;
@<RC\_CovarianceX Prototype@>;
@<RC\_order\_subset\_wrt\_block Prototype@>;
@}

The \proglang{C} file \verb|Sums.c| defines the \proglang{C}
functions and a corresponding \proglang{R} interface (via \verb|.C()|)

@o Sums.c -cc
@{
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#define S(i, j, n) ((i) >= (j) ? (n) * (j) + (i) - (j) * ((j) + 1) / 2 : (n) * (i) + (j) - (i) * ((i) + 1) / 2)
#define DoSymmetric 1
#define DoCenter 1
#define Power1 1
#define Power2 2
#define Offset0 0
@<Function Definitions@>
@}

@d Function Definitions
@{

int NCOL
(
    SEXP x
) {

    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    return(INTEGER(a)[1]);
}
@<C\_KronSums\_dweights\_dsubset@>
@<C\_KronSums\_iweights\_dsubset@>
@<C\_KronSums\_iweights\_isubset@>
@<C\_KronSums\_dweights\_isubset@>
@<C\_XfactorKronSums\_dweights\_dsubset@>
@<C\_XfactorKronSums\_iweights\_dsubset@>
@<C\_XfactorKronSums\_iweights\_isubset@>
@<C\_XfactorKronSums\_dweights\_isubset@>
@<RC\_KronSums@>
@<R\_KronSums@>
@<C\_colSums\_dweights\_dsubset@>
@<C\_colSums\_iweights\_dsubset@>
@<C\_colSums\_iweights\_isubset@>
@<C\_colSums\_dweights\_isubset@>
@<RC\_colSums@>
@<R\_colSums@>
@<C\_Sums\_dweights\_dsubset@>
@<C\_Sums\_iweights\_dsubset@>
@<C\_Sums\_iweights\_isubset@>
@<C\_Sums\_dweights\_isubset@>
@<RC\_Sums@>
@<R\_Sums@>
@<C\_KronSums\_Permutation\_isubset@>
@<C\_KronSums\_Permutation\_dsubset@>
@<C\_XfactorKronSums\_Permutation\_isubset@>
@<C\_XfactorKronSums\_Permutation\_dsubset@>
@<RC\_KronSums\_Permutation@>
@<R\_KronSums\_Permutation@>
@<C\_OneTableSums\_dweights\_dsubset@>
@<C\_OneTableSums\_iweights\_dsubset@>
@<C\_OneTableSums\_iweights\_isubset@>
@<C\_OneTableSums\_dweights\_isubset@>
@<RC\_OneTableSums@>
@<R\_OneTableSums@>
@<C\_TwoTableSums\_dweights\_dsubset@>
@<C\_TwoTableSums\_iweights\_dsubset@>
@<C\_TwoTableSums\_iweights\_isubset@>
@<C\_TwoTableSums\_dweights\_isubset@>
@<RC\_TwoTableSums@>
@<R\_TwoTableSums@>
@<C\_ThreeTableSums\_dweights\_dsubset@>
@<C\_ThreeTableSums\_iweights\_dsubset@>
@<C\_ThreeTableSums\_iweights\_isubset@>
@<C\_ThreeTableSums\_dweights\_isubset@>
@<RC\_ThreeTableSums@>
@<R\_ThreeTableSums@>
@<RC\_LinearStatistic@>
@<R\_LinearStatistic@>
@<RC\_ExpectationInfluence@>
@<R\_ExpectationInfluence@>
@<RC\_CovarianceInfluence@>
@<R\_CovarianceInfluence@>
@<RC\_ExpectationX@>
@<R\_ExpectationX@>
@<RC\_CovarianceX@>
@<R\_CovarianceX@>
@<C\_setup\_subset@>
@<C\_setup\_subset\_block@>
@<C\_order\_subset\_wrt\_block@>
@<RC\_order\_subset\_wrt\_block@>;
@<R\_order\_subset\_wrt\_block@>;
@}

The \proglang{R} interfaces are used to implement
regression tests to be called from within \proglang{R}

@d R N Input
@{
    SEXP N,
@}

@d C integer N Input
@{
    R_xlen_t N,
@}

@d R x Input
@{
    SEXP x,
@}

@d C integer x Input
@{
    int *x,
    @<C integer N Input@>
    int P,
@}

@d C real x Input
@{
    double *x,
    @<C integer N Input@>
    int P,
@}

@d R y Input
@{
    SEXP y,
@}

@d C integer y Input
@{
    int *y,
    int Q,
@}

@d C real y Input
@{
    double *y,
    int Q,
@}

@d R weights Input
@{
    SEXP weights,
@}

@d C integer weights Input
@{
    int *weights,
    int HAS_WEIGHTS,
@}

@d C real weights Input
@{
    double *weights,
    int HAS_WEIGHTS,
@}

@d R subset Input
@{
    SEXP subset
@}

@d C integer subset Input
@{
    int *subset,
    const R_xlen_t offset,
    const R_xlen_t Nsubset
@}

@d C real subset Input
@{
    double *subset,
    const R_xlen_t offset,
    const R_xlen_t Nsubset
@}

@d R block Input
@{
    SEXP block,
@}

@d R blockTable Input
@{
    SEXP blockTable
@}

@d C integer block Input
@{
    int *block,
    int L,
@}


@d init subset loop
@{
    R_xlen_t diff = 0;
    s = subset;
    w = weights;
    if (Nsubset > 0)
        diff = (R_xlen_t) s[offset];
@}

@d start subset loop
@{
    for (R_xlen_t i = offset; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
@}

@d continue subset loop
@{
    if (Nsubset > 0) {
        diff = (R_xlen_t) s[1] - s[0];
        s++;
    } else {
        diff = 1;
    }
@}

<<Sums-setup>>=
### replace with library("libcoin")
dyn.load("Sums.so")
set.seed(29)
N <- 20L
P <- 3L
Q <- 4L
L <- 5L
x <- matrix(runif(N * P), nrow = N)
y <- matrix(runif(N * Q), nrow = N)
ix <- sample(1:P, size = N, replace = TRUE)
iX <- model.matrix(~ as.factor(ix) - 1)
iy <- sample(1:Q, size = N, replace = TRUE)
weights <- sample(0:5, size = N, replace = TRUE)
block <- sample(gl(L, ceiling(N / L))[1:N])
subset <- sort(sample(1:N, floor(N * 1.5), replace = TRUE))
subsety <- sample(1:N, floor(N * 1.5), replace = TRUE)
r1 <- rep(1:ncol(x), ncol(y))
r2 <- rep(1:ncol(y), each = ncol(x))
@@

\subsection{Linear Statistics}

<<LinearStatistics>>=
a0 <- colSums(x[subset,r1] * y[subset,r2] * weights[subset])
a1 <- .Call("R_LinearStatistic", x, P, y, weights, subset - 1L, integer(0))
a2 <- .Call("R_LinearStatistic", x, P, y, as.double(weights), as.double(subset - 1L), integer(0))
a3 <- .Call("R_LinearStatistic", x, P, y, weights, as.double(subset - 1L), integer(0))
a4 <- .Call("R_LinearStatistic", x, P, y, as.double(weights), subset - 1L, integer(0))

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))


a0 <- as.vector(colSums(iX[subset,r1] * y[subset,r2] * weights[subset]))
a1 <- .Call("R_LinearStatistic", ix, P, y, weights, subset - 1L, integer(0))
a2 <- .Call("R_LinearStatistic", ix, P, y, as.double(weights), as.double(subset - 1L), integer(0))
a3 <- .Call("R_LinearStatistic", ix, P, y, weights, as.double(subset - 1L), integer(0))
a4 <- .Call("R_LinearStatistic", ix, P, y, as.double(weights), subset - 1L, integer(0))

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a0 <- colSums(x[subset,r1] * y[subsety, r2])
a1 <- .Call("R_LinearStatistic", x, P, y, integer(0), subset - 1L, subsety -1L)
a2 <- .Call("R_LinearStatistic", x, P, y, integer(0), as.double(subset - 1L), as.double(subsety -1L))
stopifnot(all.equal(a0, a1) && all.equal(a0, a1))

a0 <- as.vector(colSums(iX[subset,r1] * y[subsety, r2]))
a1 <- .Call("R_LinearStatistic", ix, P, y, integer(0), subset - 1L, subsety -1L)
a1 <- .Call("R_LinearStatistic", ix, P, y, integer(0), as.double(subset - 1L), as.double(subsety -1L))
stopifnot(all.equal(a0, a1))

@@

@d R\_LinearStatistic
@{
SEXP R_LinearStatistic
(
    @<R x Input@>
    SEXP P,
    @<R y Input@>
    @<R weights Input@>
    @<R subset Input@>,
    SEXP subsety
) 
{
    SEXP ans;
    int Q;
    R_xlen_t N, Nsubset;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * Q));
    RC_LinearStatistic(x, N, INTEGER(P)[0], REAL(y), Q, 
                       weights, subset, Offset0, Nsubset, subsety, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

@d RC\_LinearStatistic Prototype
@{
void RC_LinearStatistic
(
    @<R x Input@>
    R_xlen_t N,
    int P,
    @<C real y Input@>
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    SEXP subsety,
    @<C KronSums Answer@>
) 
@}

@d RC\_LinearStatistic
@{
@<RC\_LinearStatistic Prototype@>
{
    double center;

    if (XLENGTH(subsety) == 0) {
        RC_KronSums(x, N, P, y, Q, !DoSymmetric, &center, &center, !DoCenter, weights, 
                    subset, Offset0, Nsubset, PQ_ans);
    } else {
        if (XLENGTH(weights) > 0) 
            error("weights given for permutation");
        if (XLENGTH(subset) != XLENGTH(subsety))
            error("incorrect subsets");
        RC_KronSums_Permutation(x, N, P, y, Q, subset, Offset0, Nsubset, 
                                subsety, PQ_ans);
    }
}
@|RC_LinearStatistic
@}

\subsection{Expectation and Covariance Influence}

<<ExpectationCovarianceInfluence>>=
sw <- sum(weights[subset])
expecty <- a0 <- colSums(y[subset, ] * weights[subset]) / sw
a1 <- .Call("R_ExpectationInfluence", y, weights, subset - 1L);
a2 <- .Call("R_ExpectationInfluence", y, as.double(weights), as.double(subset - 1L));
a3 <- .Call("R_ExpectationInfluence", y, weights, as.double(subset - 1L));
a4 <- .Call("R_ExpectationInfluence", y, as.double(weights), subset - 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@

@d R\_ExpectationInfluence
@{
SEXP R_ExpectationInfluence
(
    @<R y Input@>
    @<R weights Input@>
    @<R subset Input@>
) 
{
    SEXP ans;
    int Q;
    R_xlen_t N, Nsubset;
    double sw;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    sw = RC_Sums(N, weights, subset, Offset0, Nsubset);

    PROTECT(ans = allocVector(REALSXP, Q));
    RC_ExpectationInfluence(N, REAL(y), Q, weights, subset, Offset0, Nsubset, sw, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

@d RC\_ExpectationInfluence Prototype
@{
void RC_ExpectationInfluence
(
    @<C integer N Input@>
    @<C real y Input@>
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    double sumweights,
    @<C colSums Answer@>
) 
@}

@d RC\_ExpectationInfluence
@{
@<RC\_ExpectationInfluence Prototype@>
{
    double center;

    RC_colSums(y, N, Q, DoSymmetric, &center, !DoCenter, weights, 
               subset, Offset0, Nsubset, P_ans);
    for (int q = 0; q < Q; q++) 
        P_ans[q] = P_ans[q] / sumweights;
}
@|RC_ExpectationInfluence
@}

<<CovarianceInfluence>>=
sw <- sum(weights[subset])
yc <- t(t(y) - expecty)
r1y <- rep(1:ncol(y), ncol(y))
r2y <- rep(1:ncol(y), each = ncol(y))
a0 <- colSums(yc[subset, r1y] * yc[subset, r2y] * weights[subset]) / sw
a0 <- matrix(a0, ncol = ncol(y))
vary <- diag(a0)
a0 <- a0[upper.tri(a0, diag = TRUE)]
a1 <- .Call("R_CovarianceInfluence", y, weights, subset - 1L, 0L);
a2 <- .Call("R_CovarianceInfluence", y, as.double(weights), as.double(subset - 1L), 0L);
a3 <- .Call("R_CovarianceInfluence", y, weights, as.double(subset - 1L), 0L);
a4 <- .Call("R_CovarianceInfluence", y, as.double(weights), subset - 1L, 0L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a1 <- .Call("R_CovarianceInfluence", y, weights, subset - 1L, 1L);
a2 <- .Call("R_CovarianceInfluence", y, as.double(weights), as.double(subset - 1L), 1L);
a3 <- .Call("R_CovarianceInfluence", y, weights, as.double(subset - 1L), 1L);
a4 <- .Call("R_CovarianceInfluence", y, as.double(weights), subset - 1L, 1L);

a0 <- vary

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

@@

@d R\_CovarianceInfluence
@{
SEXP R_CovarianceInfluence
(
    @<R y Input@>
    @<R weights Input@>
    @<R subset Input@>,
    SEXP varonly
) 
{
    SEXP ans;
    SEXP ExpInf;
    int Q;
    R_xlen_t N, Nsubset;
    double sw;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    PROTECT(ExpInf = R_ExpectationInfluence(y, weights, subset));

    sw = RC_Sums(N, weights, subset, Offset0, Nsubset);

    if (INTEGER(varonly)[0]) {
        PROTECT(ans = allocVector(REALSXP, Q));
    } else {
        PROTECT(ans = allocVector(REALSXP, Q * (Q + 1) / 2));
    }
    RC_CovarianceInfluence(N, y, Q, weights, subset, Offset0, Nsubset, REAL(ExpInf), sw, 
                           INTEGER(varonly)[0], REAL(ans));
    UNPROTECT(2);
    return(ans);
}
@}

@d RC\_CovarianceInfluence Prototype
@{
void RC_CovarianceInfluence
(
    @<C integer N Input@>
    @<R y Input@>
    int Q,
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    double *ExpInf,
    double sumweights,
    int VARONLY,
    @<C KronSums Answer@>
) 
@}

@d RC\_CovarianceInfluence
@{
@<RC\_CovarianceInfluence Prototype@>
{
    if (VARONLY) {
        RC_colSums(REAL(y), N, Q, Power2, ExpInf, DoCenter, weights, 
                   subset, offset, Nsubset, PQ_ans);
        for (int q = 0; q < Q; q++) 
            PQ_ans[q] = PQ_ans[q] / sumweights;
    } else {
        RC_KronSums(y, N, Q, REAL(y), Q, DoSymmetric, ExpInf, ExpInf, DoCenter, weights, 
                    subset, offset, Nsubset, PQ_ans);
        for (int q = 0; q < Q * (Q + 1) / 2; q++) 
            PQ_ans[q] = PQ_ans[q] / sumweights;
    }
}
@|RC_CovarianceInfluence
@}

\subsection{Expectation and Covariance X}

<<ExpectationCovarianceX>>=
expectx <- a0 <- colSums(x[subset, ] * weights[subset]) 
a0
a1 <- .Call("R_ExpectationX", x, P, weights, subset - 1L);
a2 <- .Call("R_ExpectationX", x, P, as.double(weights), as.double(subset - 1L));
a3 <- .Call("R_ExpectationX", x, P, weights, as.double(subset - 1L));
a4 <- .Call("R_ExpectationX", x, P, as.double(weights), subset - 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

expectx <- a0 <- colSums(x[subset, ]^2 * weights[subset]) 
a1 <- .Call("R_CovarianceX", x, P, weights, subset - 1L, 1L);
a2 <- .Call("R_CovarianceX", x, P, as.double(weights), as.double(subset - 1L), 1L);
a3 <- .Call("R_CovarianceX", x, P, weights, as.double(subset - 1L), 1L);
a4 <- .Call("R_CovarianceX", x, P, as.double(weights), subset - 1L, 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

expectx <- a0 <- as.vector(colSums(iX[subset, ] * weights[subset]))
a0
a1 <- .Call("R_ExpectationX", ix, P + 1L, weights, subset - 1L)[-1];
a2 <- .Call("R_ExpectationX", ix, P + 1L, as.double(weights), as.double(subset - 1L))[-1];
a3 <- .Call("R_ExpectationX", ix, P + 1L, weights, as.double(subset - 1L))[-1];
a4 <- .Call("R_ExpectationX", ix, P + 1L, as.double(weights), subset - 1L)[-1];

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a1 <- .Call("R_CovarianceX", ix, P + 1L, weights, subset - 1L, 1L)[-1];
a2 <- .Call("R_CovarianceX", ix, P + 1L, as.double(weights), as.double(subset - 1L), 1L)[-1];
a3 <- .Call("R_CovarianceX", ix, P + 1L, weights, as.double(subset - 1L), 1L)[-1];
a4 <- .Call("R_CovarianceX", ix, P + 1L, as.double(weights), subset - 1L, 1L)[-1];

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

r1x <- rep(1:ncol(x), ncol(x))
r2x <- rep(1:ncol(x), each = ncol(x))
a0 <- colSums(iX[subset, r1x] * iX[subset, r2x] * weights[subset])
a0 <- matrix(a0, ncol = ncol(x))
vary <- diag(a0)
a0 <- a0[upper.tri(a0, diag = TRUE)]

a0
a1 <- .Call("R_CovarianceX", ix, P + 1L, weights, subset - 1L, 0L)[-1];
a1
a2 <- .Call("R_CovarianceX", ix, P + 1L, as.double(weights), as.double(subset - 1L), 0L)[-1];
a3 <- .Call("R_CovarianceX", ix, P + 1L, weights, as.double(subset - 1L), 0L)[-1];
a4 <- .Call("R_CovarianceX", ix, P + 1L, as.double(weights), subset - 1L, 0L)[-1];

#stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
#          all.equal(a0, a3) && all.equal(a0, a4))


@@

@d R\_ExpectationX
@{
SEXP R_ExpectationX
(
    @<R x Input@>
    SEXP P,
    @<R weights Input@>
    @<R subset Input@>
) 
{
    SEXP ans;
    R_xlen_t N, Nsubset;

    N = XLENGTH(x) / INTEGER(P)[0];
    Nsubset = XLENGTH(subset);

    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0]));
    RC_ExpectationX(x, N, INTEGER(P)[0], weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

@d RC\_ExpectationX Prototype
@{
void RC_ExpectationX
(
    @<R x Input@>
    R_xlen_t N,
    int P,
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<C OneTableSums Answer@>
) 
@}

@d RC\_ExpectationX
@{
@<RC\_ExpectationX Prototype@>
{
    double center;

    if (TYPEOF(x) == INTSXP) {
        RC_OneTableSums(INTEGER(x), N, P, weights, subset, offset, Nsubset, P_ans);
    } else {
        RC_colSums(REAL(x), N, P, DoSymmetric, &center, !DoCenter, weights, subset, offset, Nsubset, P_ans);
    }
}
@|RC_ExpectationX
@}

@d R\_CovarianceX
@{
SEXP R_CovarianceX
(
    @<R x Input@>
    SEXP P,
    @<R weights Input@>
    @<R subset Input@>,
    SEXP varonly
) 
{
    SEXP ans;
    SEXP ExpX;
    R_xlen_t N, Nsubset;

    N = XLENGTH(x) / INTEGER(P)[0];
    Nsubset = XLENGTH(subset);

    PROTECT(ExpX = R_ExpectationX(x, P, weights, subset));

    if (INTEGER(varonly)[0]) {
        PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0]));
    } else {
        PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * (INTEGER(P)[0] + 1) / 2));
    }
    RC_CovarianceX(x, N, INTEGER(P)[0], weights, subset, Offset0, Nsubset, REAL(ExpX), 
                   INTEGER(varonly)[0], REAL(ans));
    UNPROTECT(2);
    return(ans);
}
@}

@d RC\_CovarianceX Prototype
@{
void RC_CovarianceX
(
    @<R x Input@>
    R_xlen_t N,
    int P,
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    double *ExpX,
    int VARONLY,
    @<C KronSums Answer@>
) 
@}

@d RC\_CovarianceX
@{
@<RC\_CovarianceX Prototype@>
{
    double center;

    if (TYPEOF(x) == INTSXP) {
        if (VARONLY) {
            for (int p = 0; p < P; p++) PQ_ans[p] = ExpX[p];
        } else {
            for (int p = 0; p < P; p++)
                PQ_ans[S(p, p, P)] = ExpX[p];
        }
    } else {
        if (VARONLY) {
            RC_colSums(REAL(x), N, P, Power2, &center, !DoCenter, weights, 
                       subset, offset, Nsubset, PQ_ans);
        } else {
            RC_KronSums(x, N, P, REAL(x), P, DoSymmetric, &center, &center, !DoCenter, weights, 
                        subset, offset, Nsubset, PQ_ans);
        }
    }
}
@|RC_CovarianceX
@}


\subsection{Sums}

<<Sums>>=
a0 <- sum(weights[subset])
a1 <- .Call("R_Sums", N, weights, subset - 1L)
a2 <- .Call("R_Sums", N, as.double(weights), as.double(subset - 1L))
a3 <- .Call("R_Sums", N, weights, as.double(subset - 1L))
a4 <- .Call("R_Sums", N, as.double(weights), subset - 1L)
stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@

@d R\_Sums
@{
SEXP R_Sums
(
    @<R N Input@>
    @<R weights Input@>
    @<R subset Input@>
) 
{
    SEXP ans;
    R_xlen_t Nsubset;

    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = RC_Sums(INTEGER(N)[0], weights, subset, Offset0, Nsubset);
    UNPROTECT(1);

    return(ans);
}
@}

@d RC\_Sums Prototype
@{
double RC_Sums
(
    @<C integer N Input@>
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset
) 
@}

@d RC\_Sums
@{
@<RC\_Sums Prototype@>
{
    if (XLENGTH(weights) == 0) {
        if (XLENGTH(subset) == 0) {
            return((double) N);
        } else {
            return((double) Nsubset - offset);
        }
    } 
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            return(C_Sums_iweights_isubset(N, INTEGER(weights), XLENGTH(weights),
                                           INTEGER(subset), offset, Nsubset));
        } else {
            return(C_Sums_iweights_dsubset(N, INTEGER(weights), XLENGTH(weights), 
                                           REAL(subset), offset, Nsubset));
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            return(C_Sums_dweights_isubset(N, REAL(weights), XLENGTH(weights),
                                           INTEGER(subset), offset, Nsubset));
        } else {
            return(C_Sums_dweights_dsubset(N, REAL(weights), XLENGTH(weights),
                                           REAL(subset), offset, Nsubset));
        }
    }
}
@|RC_Sums
@}


@d C\_Sums\_dweights\_dsubset
@{
double C_Sums_dweights_dsubset
(
    @<C integer N Input@>
    @<C real weights Input@>
    @<C real subset Input@>
) 
{
    double *s, *w; 
    @<Sums Body@>
}
@|C_Sums_dweights_dsubset
@}

@d C\_Sums\_iweights\_dsubset
@{
double C_Sums_iweights_dsubset
(
    @<C integer N Input@>
    @<C integer weights Input@>
    @<C real subset Input@>
) 
{
    double *s;
    int *w; 
    @<Sums Body@>
}
@|C_Sums_iweights_dsubset
@}

@d C\_Sums\_iweights\_isubset
@{
double C_Sums_iweights_isubset
(
    @<C integer N Input@>
    @<C integer weights Input@>
    @<C integer subset Input@>
) 
{
    int *s, *w;
    @<Sums Body@>
}
@|C_Sums_iweights_isubset
@}

@d C\_Sums\_dweights\_isubset
@{
double C_Sums_dweights_isubset
(
    @<C integer N Input@>
    @<C real weights Input@>
    @<C integer subset Input@>
) 
{
    int *s; 
    double *w;
    @<Sums Body@>
}
@|C_Sums_dweights_isubset
@}


@d Sums Body
@{

    double ans = 0.0;

    if (Nsubset > 0) {
        if (!HAS_WEIGHTS) return((double) Nsubset - offset);
    } else {
        if (!HAS_WEIGHTS) return((double) N);
    }

    @<init subset loop@>    

    @<start subset loop@>    
    {
        w = w + diff;
        ans += w[0];
        @<continue subset loop@>    
    }
    w = w + diff;
    ans += w[0];
    return(ans);
@}


\subsection{Kronecker Sums}

<<KronSums>>=
r1 <- rep(1:ncol(x), ncol(y))
r2 <- rep(1:ncol(y), each = ncol(x))

a0 <- colSums(x[subset,r1] * y[subset,r2] * weights[subset])
a1 <- .Call("R_KronSums", x, P, y, weights, subset - 1L)
a2 <- .Call("R_KronSums", x, P, y, as.double(weights), as.double(subset - 1L))
a3 <- .Call("R_KronSums", x, P, y, weights, as.double(subset - 1L))
a4 <- .Call("R_KronSums", x, P, y, as.double(weights), subset - 1L)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@

@d R\_KronSums
@{
SEXP R_KronSums
(
    @<R x Input@>
    SEXP P,
    @<R y Input@>
    @<R weights Input@>
    @<R subset Input@>
) 
{
    SEXP ans;
    int Q;
    R_xlen_t N, Nsubset;
    double center;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * Q));
    RC_KronSums(x, N, INTEGER(P)[0], REAL(y), Q, !DoSymmetric, &center, &center, 
                !DoCenter, weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

@d RC\_KronSums Prototype
@{
void RC_KronSums
(
    @<RC KronSums Input@>
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<C KronSums Answer@>
) 
@}

@d RC\_KronSums
@{
@<RC\_KronSums Prototype@>
{
    if (TYPEOF(x) == INTSXP) {
        if (SYMMETRIC) error("not implemented");
        if (CENTER) error("not implemented");
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_XfactorKronSums_iweights_isubset(INTEGER(x), N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_XfactorKronSums_iweights_dsubset(INTEGER(x), N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_XfactorKronSums_dweights_isubset(INTEGER(x), N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_XfactorKronSums_dweights_dsubset(INTEGER(x), N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    }
    } else {
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_KronSums_iweights_isubset(REAL(x), N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_KronSums_iweights_dsubset(REAL(x), N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_KronSums_dweights_isubset(REAL(x), N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_KronSums_dweights_dsubset(REAL(x), N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    }
    }
}
@|RC_KronSums
@}

@d RC KronSums Input
@{
    @<R x Input@>
    R_xlen_t N,
    int P,
    @<C real y Input@>
    const int SYMMETRIC,
    double *centerx,
    double *centery,
    const int CENTER,
@}

@d C KronSums Input
@{
    @<C real x Input@>
    @<C real y Input@>
    const int SYMMETRIC,
    double *centerx,
    double *centery,
    const int CENTER,
@}

@d C KronSums Answer
@{
    double *PQ_ans
@}

@d C\_KronSums\_dweights\_dsubset
@{
void C_KronSums_dweights_dsubset
(
    @<C KronSums Input@>
    @<C real weights Input@>
    @<C real subset Input@>,
    @<C KronSums Answer@>
)
{
    double *s, *w; 
    @<KronSums Body@>
}
@|C_KronSums_dweights_dsubset
@}

@d C\_KronSums\_iweights\_dsubset
@{
void C_KronSums_iweights_dsubset
(
    @<C KronSums Input@>
    @<C integer weights Input@>    
    @<C real subset Input@>,
    @<C KronSums Answer@>
) 
{
    double *s;
    int *w; 
    @<KronSums Body@>
}
@|C_KronSums_iweights_dsubset
@}

@d C\_KronSums\_iweights\_isubset
@{
void C_KronSums_iweights_isubset
(
    @<C KronSums Input@>
    @<C integer weights Input@>    
    @<C integer subset Input@>,
    @<C KronSums Answer@>
) 
{
    int *s, *w;
    @<KronSums Body@>
}
@|C_KronSums_iweights_isubset
@}

@d C\_KronSums\_dweights\_isubset
@{
void C_KronSums_dweights_isubset
(
    @<C KronSums Input@>
    @<C real weights Input@>    
    @<C integer subset Input@>,
    @<C KronSums Answer@>
) {
    int *s; 
    double *w;
    @<KronSums Body@>
}
@|C_KronSums_dweights_isubset
@}

@d KronSums Body
@{
    double *xx, *yy, cx = 0.0, cy = 0.0;

    for (int q = 0; q < Q; q++) {
        for (int p = 0; p < (SYMMETRIC ? q + 1 : P); p++) {
            PQ_ans[0] = 0.0;
            yy = y + N * q;
            xx = x + N * p;
            if (CENTER) {
                cx = centerx[p];
                cy = centery[q];
            }
            @<init subset loop@>
            @<start subset loop@>
            {
                xx = xx + diff;
                yy = yy + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    PQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                } else {
                    PQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                }
                @<continue subset loop@>
            }
            xx = xx + diff;
            yy = yy + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
            } else {
                PQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
            }
            PQ_ans++;
        }
    }
@}

\subsection{Xfactor Kronecker Sums}

<<XfactorKronSums>>=

a0 <- as.vector(colSums(iX[subset,r1] * 
                        y[subset,r2] * weights[subset]))
a1 <- .Call("R_KronSums", ix, P, y, weights, subset - 1L)
a2 <- .Call("R_KronSums", ix, P, y, as.double(weights), as.double(subset - 1L))
a3 <- .Call("R_KronSums", ix, P, y, weights, as.double(subset - 1L))
a4 <- .Call("R_KronSums", ix, P, y, as.double(weights), subset - 1L)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@

@d C XfactorKronSums Input
@{
    @<C integer x Input@>
    @<C real y Input@>
@}

@d C\_XfactorKronSums\_dweights\_dsubset
@{
void C_XfactorKronSums_dweights_dsubset
(
    @<C XfactorKronSums Input@>
    @<C real weights Input@>
    @<C real subset Input@>,
    @<C KronSums Answer@>
)
{
    double *s, *w; 
    @<XfactorKronSums Body@>
}
@|C_XfactorKronSums_dweights_dsubset
@}

@d C\_XfactorKronSums\_iweights\_dsubset
@{
void C_XfactorKronSums_iweights_dsubset
(
    @<C XfactorKronSums Input@>
    @<C integer weights Input@>    
    @<C real subset Input@>,
    @<C KronSums Answer@>
) 
{
    double *s;
    int *w; 
    @<XfactorKronSums Body@>
}
@|C_XfactorKronSums_iweights_dsubset
@}

@d C\_XfactorKronSums\_iweights\_isubset
@{
void C_XfactorKronSums_iweights_isubset
(
    @<C XfactorKronSums Input@>
    @<C integer weights Input@>    
    @<C integer subset Input@>,
    @<C KronSums Answer@>
) 
{
    int *s, *w;
    @<XfactorKronSums Body@>
}
@|C_XfactorKronSums_iweights_isubset
@}

@d C\_XfactorKronSums\_dweights\_isubset
@{
void C_XfactorKronSums_dweights_isubset
(
    @<C XfactorKronSums Input@>
    @<C real weights Input@>    
    @<C integer subset Input@>,
    @<C KronSums Answer@>
) {
    int *s; 
    double *w;
    @<XfactorKronSums Body@>
}
@|C_XfactorKronSums_dweights_isubset
@}

@d XfactorKronSums Body
@{
    int *xx, ixi;
    double *yy;
 
    for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

    for (int q = 0; q < Q; q++) {
        yy = y + N * q;
        xx = x;
        @<init subset loop@>
        @<start subset loop@>
        {
            xx = xx + diff;
            yy = yy + diff;
            ixi = xx[0] - 1;
            if (HAS_WEIGHTS) {
                w = w + diff;
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0] * w[0];
            } else {
                if (ixi >= 0)
                    PQ_ans[ixi + q * P] += yy[0];
            }
            @<continue subset loop@>
        }
        xx = xx + diff;
        yy = yy + diff;
        ixi = xx[0] - 1;
        if (HAS_WEIGHTS) {
            w = w + diff;
            if (ixi >= 0)
                PQ_ans[ixi + q * P] += yy[0] * w[0];
        } else {
            if (ixi >= 0)
                PQ_ans[ixi + q * P] += yy[0];
        }
    }
@}


\subsection{Permuted Kronecker Sums}

<<KronSums-Permutation>>=
a0 <- colSums(x[subset,r1] * y[subsety, r2])
a1 <- .Call("R_KronSums_Permutation", x, P, y, subset - 1L, subsety -1L)
a2 <- .Call("R_KronSums_Permutation", x, P, y, as.double(subset - 1L), as.double(subsety -1L))
stopifnot(all.equal(a0, a1) && all.equal(a0, a1))
@@


@d R\_KronSums\_Permutation
@{
SEXP R_KronSums_Permutation
(
    @<R x Input@>
    SEXP P,
    @<R y Input@>
    @<R subset Input@>,
    SEXP subsety
) {

    SEXP ans;
    int Q;
    R_xlen_t N, Nsubset;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * Q));
    RC_KronSums_Permutation(x, N, INTEGER(P)[0], REAL(y), Q, subset, Offset0, Nsubset, 
                            subsety, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

@d RC\_KronSums\_Permutation Prototype
@{
void RC_KronSums_Permutation
(
    @<R x Input@>
    R_xlen_t N,
    int P,
    @<C real y Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    SEXP subsety,
    @<C KronSums Answer@>
) 
@}

@d RC\_KronSums\_Permutation
@{
@<RC\_KronSums\_Permutation Prototype@>
{
    if (TYPEOF(x) == INTSXP) {
    if (TYPEOF(subset) == INTSXP) {
        C_XfactorKronSums_Permutation_isubset(INTEGER(x), N, P, y, Q, 
                                       INTEGER(subset), offset, Nsubset, 
                                       INTEGER(subsety), PQ_ans);
        } else {
        C_XfactorKronSums_Permutation_dsubset(INTEGER(x), N, P, y, Q, 
                                       REAL(subset), offset, Nsubset, 
                                       REAL(subsety), PQ_ans);
    }
    } else {
    if (TYPEOF(subset) == INTSXP) {
        C_KronSums_Permutation_isubset(REAL(x), N, P, y, Q, 
                                       INTEGER(subset), offset, Nsubset, 
                                       INTEGER(subsety), PQ_ans);
        } else {
        C_KronSums_Permutation_dsubset(REAL(x), N, P, y, Q, 
                                       REAL(subset), offset, Nsubset, 
                                       REAL(subsety), PQ_ans);
    }
    }
}
@|RC_KronSums_Permutation
@}

@d C\_KronSums\_Permutation\_dsubset
@{
void C_KronSums_Permutation_dsubset
(
    @<C real x Input@>
    @<C real y Input@>
    @<C real subset Input@>,
    double *subsety,
    @<C KronSums Answer@>
) 
{
    double *sx, *sy;
    @<KronSums Permutation Body@>
}
@|C_KronSums_Permutation_dsubset
@}

@d C\_KronSums\_Permutation\_isubset
@{
void C_KronSums_Permutation_isubset
(
    @<C real x Input@>
    @<C real y Input@>
    @<C integer subset Input@>,
    int *subsety,
    @<C KronSums Answer@>
) 
{
    int *sx, *sy;
    @<KronSums Permutation Body@>
}
@|C_KronSums_Permutation_isubset
@}

@d KronSums Permutation Body
@{
    R_xlen_t diff;
    double *xx;

    for (int q = 0; q < Q; q++) {
        for (int p = 0; p < P; p++) {
            PQ_ans[0] = 0.0;
            xx = x + N * p;
            sx = subset;
            sy = subsety;
            diff = (R_xlen_t) sx[offset];
            for (R_xlen_t i = offset; i < Nsubset - 1; i++) {
                xx = xx + diff;
                PQ_ans[0] += xx[0] * y[ (R_xlen_t) sy[0] + q * N];
                diff = (R_xlen_t) sx[1] - sx[0];
                sx++;
                sy++;
            }
            xx = xx + diff;
            PQ_ans[0] += xx[0] * y[ (R_xlen_t) sy[0] + q * N];
            PQ_ans++;
        }
    }
@}

\subsection{Xfactor Permuted Kronecker Sums}

<<XfactorKronSums-Permutation>>=
a0 <- as.vector(colSums(iX[subset,r1] * y[subsety,
r2]))
a1 <- .Call("R_KronSums_Permutation", ix, P, y, subset - 1L, subsety -1L)
a1 <- .Call("R_KronSums_Permutation", ix, P, y, as.double(subset - 1L), as.double(subsety -1L))
stopifnot(all.equal(a0, a1))
@@


@d C\_XfactorKronSums\_Permutation\_dsubset
@{
void C_XfactorKronSums_Permutation_dsubset
(
    @<C integer x Input@>
    @<C real y Input@>
    @<C real subset Input@>,
    double *subsety,
    @<C KronSums Answer@>
) 
{
    double *sx, *sy;
    @<XfactorKronSums Permutation Body@>
}
@|C_XfactorKronSums_Permutation_dsubset
@}

@d C\_XfactorKronSums\_Permutation\_isubset
@{
void C_XfactorKronSums_Permutation_isubset
(
    @<C integer x Input@>
    @<C real y Input@>
    @<C integer subset Input@>,
    int *subsety,
    @<C KronSums Answer@>
) 
{
    int *sx, *sy;
    @<XfactorKronSums Permutation Body@>
}
@|C_XfactorKronSums_Permutation_isubset
@}

@d XfactorKronSums Permutation Body
@{
    R_xlen_t diff;
    int *xx;

    for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

    for (int q = 0; q < Q; q++) {
        xx = x;
        sx = subset;
        sy = subsety;
        diff = (R_xlen_t) sx[offset];
        for (R_xlen_t i = offset; i < Nsubset - 1; i++) {
            xx = xx + diff;
            PQ_ans[(xx[0] - 1) + q * P] += y[ (R_xlen_t) sy[0] + q * N];
            diff = (R_xlen_t) sx[1] - sx[0];
            sx++;
            sy++;
        }
        xx = xx + diff;
        PQ_ans[(xx[0] - 1) + q * P] += y[ (R_xlen_t) sy[0] + q * N];
    }
@}



\subsection{Column Sums}

<<colSums>>=
a0 <- colSums(x[subset,] * weights[subset])
a1 <- .Call("R_colSums", x, weights, subset - 1L)
a2 <- .Call("R_colSums", x, as.double(weights), as.double(subset - 1L))
a3 <- .Call("R_colSums", x, weights, as.double(subset - 1L))
a4 <- .Call("R_colSums", x, as.double(weights), subset - 1L)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@

@d R\_colSums
@{
SEXP R_colSums
(
    @<R x Input@>
    @<R weights Input@>
    @<R subset Input@>
) {

    SEXP ans;
    int P;
    R_xlen_t N, Nsubset;
    double center;

    P = NCOL(x);
    N = XLENGTH(x) / P;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, P));
    RC_colSums(REAL(x), N, P, DoSymmetric, &center, !DoCenter, weights, subset, Offset0, 
               Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

@d RC\_colSums Prototype
@{
void RC_colSums
(
    @<C colSums Input@>
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<C colSums Answer@>
) 
@}

@d RC\_colSums
@{
@<RC\_colSums Prototype@>
{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_colSums_iweights_isubset(x, N, P, power, centerx, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, P_ans);
        } else {
            C_colSums_iweights_dsubset(x, N, P, power, centerx, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, P_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_colSums_dweights_isubset(x, N, P, power, centerx, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, P_ans);
        } else {
            C_colSums_dweights_dsubset(x, N, P, power, centerx, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, P_ans);
        }
    }
}
@|RC_colSums
@}

@d C colSums Input
@{
    @<C real x Input@>
    const int power,
    double *centerx,
    const int CENTER,
@}

@d C colSums Answer
@{
    double *P_ans
@}

@d C\_colSums\_dweights\_dsubset
@{
void C_colSums_dweights_dsubset
(
    @<C colSums Input@>
    @<C real weights Input@>    
    @<C real subset Input@>,
    @<C colSums Answer@>
) 
{
    double *s, *w; 
    @<colSums Body@>
}
@|C_colSums_dweights_dsubset
@}

@d C\_colSums\_iweights\_dsubset
@{
void C_colSums_iweights_dsubset
(
    @<C colSums Input@>
    @<C integer weights Input@>    
    @<C real subset Input@>,
    @<C colSums Answer@>
) 
{
    double *s;
    int *w; 
    @<colSums Body@>
}
@|C_colSums_iweights_dsubset
@}

@d C\_colSums\_iweights\_isubset
@{
void C_colSums_iweights_isubset
(
    @<C colSums Input@>
    @<C integer weights Input@>    
    @<C integer subset Input@>,
    @<C colSums Answer@>
) 
{
    int *s, *w;
    @<colSums Body@>
}
@|C_colSums_iweights_isubset
@}

@d C\_colSums\_dweights\_isubset
@{
void C_colSums_dweights_isubset
(
    @<C colSums Input@>
    @<C real weights Input@>    
    @<C integer subset Input@>,
    @<C colSums Answer@>
) 
{
    int *s; 
    double *w;
    @<colSums Body@>
}
@|C_colSums_dweights_isubset
@}

@d colSums Body
@{
    double *xx, cx = 0.0;

    for (int p = 0; p < P; p++) {
        P_ans[0] = 0.0;
        xx = x + N * p;
        if (CENTER) {
            cx = centerx[p];
        }
        @<init subset loop@>
        @<start subset loop@>
        {
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[0] += pow(xx[0] - cx, power) * w[0];
            } else {
                P_ans[0] += pow(xx[0] - cx, power);
            }
            @<continue subset loop@>
        }
        xx = xx + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            P_ans[0] += pow(xx[0] - cx, power) * w[0];
        } else {
            P_ans[0] += pow(xx[0] - cx, power);
        }
        P_ans++;
    }
@}

\subsection{OneTable Sums}

<<OneTableSum>>=

a0 <- as.vector(xtabs(weights ~ ix, subset = subset))
a1 <- .Call("R_OneTableSums", ix, P + 1L, weights, subset - 1L)[-1]
a2 <- .Call("R_OneTableSums", ix, P + 1L, 
            as.double(weights), as.double(subset - 1L))[-1]
a3 <- .Call("R_OneTableSums", ix, P + 1L, 
            weights, as.double(subset - 1L))[-1]
a4 <- .Call("R_OneTableSums", ix, P + 1L, 
            as.double(weights), subset - 1L)[-1]

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@



@d R\_OneTableSums
@{
SEXP R_OneTableSums
(
    @<R x Input@>
    SEXP P,
    @<R weights Input@>
    @<R subset Input@>
) {

    SEXP ans;
    R_xlen_t N, Nsubset;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0]));
    RC_OneTableSums(INTEGER(x), N, INTEGER(P)[0],  
                 weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}


@d RC\_OneTableSums Prototype
@{
void RC_OneTableSums
(
    @<C OneTableSums Input@>
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<C OneTableSums Answer@>
) 
@}

@d RC\_OneTableSums
@{
@<RC\_OneTableSums Prototype@>
{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_OneTableSums_iweights_isubset(x, N, P, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, P_ans);
        } else {
            C_OneTableSums_iweights_dsubset(x, N, P, 
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, P_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_OneTableSums_dweights_isubset(x, N, P, 
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, P_ans);
        } else {
            C_OneTableSums_dweights_dsubset(x, N, P, 
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, P_ans);
        }
    }
}
@|RC_OneTableSums
@}

@d C OneTableSums Input
@{
    @<C integer x Input@>
@}

@d C OneTableSums Answer
@{
    double *P_ans
@}

@d C\_OneTableSums\_dweights\_dsubset
@{
void C_OneTableSums_dweights_dsubset
(
    @<C OneTableSums Input@>
    @<C real weights Input@>
    @<C real subset Input@>,
    @<C OneTableSums Answer@>
)
{
    double *s, *w; 
    @<OneTableSums Body@>
}
@|C_OneTableSums_dweights_dsubset
@}

@d C\_OneTableSums\_iweights\_dsubset
@{
void C_OneTableSums_iweights_dsubset
(
    @<C OneTableSums Input@>
    @<C integer weights Input@>
    @<C real subset Input@>,
    @<C OneTableSums Answer@>
)
{
    double *s;
    int *w; 
    @<OneTableSums Body@>
}
@|C_OneTableSums_iweights_dsubset
@}

@d C\_OneTableSums\_iweights\_isubset
@{
void C_OneTableSums_iweights_isubset
(
    @<C OneTableSums Input@>
    @<C integer weights Input@>    
    @<C integer subset Input@>,
    @<C OneTableSums Answer@>
)
{
    int *s, *w;
    @<OneTableSums Body@>
}
@|C_OneTableSums_iweights_isubset
@}

@d C\_OneTableSums\_dweights\_isubset
@{
void C_OneTableSums_dweights_isubset
(
    @<C OneTableSums Input@>
    @<C real weights Input@>
    @<C integer subset Input@>,
    @<C OneTableSums Answer@>
)
{
    int *s; 
    double *w;
    @<OneTableSums Body@>
}
@|C_OneTableSums_dweights_isubset
@}

@d OneTableSums Body
@{
    int *xx;

    for (int p = 0; p < P; p++) P_ans[p] = 0.0;

    xx = x;
    @<init subset loop@>
    @<start subset loop@>
    {
        xx = xx + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            P_ans[xx[0]] += (double) w[0];
        } else {
            P_ans[xx[0]]++;
        }
        @<continue subset loop@>
    }
    xx = xx + diff;
    if (HAS_WEIGHTS) {
        w = w + diff;
        P_ans[xx[0]] += w[0];
    } else {
        P_ans[xx[0]]++;
    }
@}

\subsection{TwoTable Sums}

<<TwoTableSum>>=

a0 <- c(xtabs(weights ~ ix + iy, subset = subset))
a1 <- .Call("R_TwoTableSums", ix, P + 1L, iy, Q + 1L, 
            weights, subset - 1L)
a1 <- c(matrix(a1, nrow = P + 1, ncol = Q + 1)[-1,-1])
a2 <- .Call("R_TwoTableSums", ix, P + 1L, iy, Q + 1L, 
            as.double(weights), as.double(subset - 1L))
a2 <- c(matrix(a2, nrow = P + 1, ncol = Q + 1)[-1,-1])
a3 <- .Call("R_TwoTableSums", ix, P + 1L, iy, Q + 1L, 
            weights, as.double(subset - 1L))
a3 <- c(matrix(a3, nrow = P + 1, ncol = Q + 1)[-1,-1])
a4 <- .Call("R_TwoTableSums", ix, P + 1L, iy, Q + 1L, 
            as.double(weights), subset - 1L)
a4 <- c(matrix(a4, nrow = P + 1, ncol = Q + 1)[-1,-1])

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@



@d R\_TwoTableSums
@{
SEXP R_TwoTableSums
(
    @<R x Input@>
    SEXP P,
    @<R y Input@>
    SEXP Q,
    @<R weights Input@>
    @<R subset Input@>
) {

    SEXP ans;
    R_xlen_t N, Nsubset;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * INTEGER(Q)[0]));
    RC_TwoTableSums(INTEGER(x), N, INTEGER(P)[0], INTEGER(y), INTEGER(Q)[0], 
                 weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}


@d RC\_TwoTableSums Prototype
@{
void RC_TwoTableSums
(
    @<C TwoTableSums Input@>
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<C TwoTableSums Answer@>
) 
@}

@d RC\_TwoTableSums
@{
@<RC\_TwoTableSums Prototype@>
{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_TwoTableSums_iweights_isubset(x, N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_TwoTableSums_iweights_dsubset(x, N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_TwoTableSums_dweights_isubset(x, N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_TwoTableSums_dweights_dsubset(x, N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    }
}
@|RC_TwoTableSums
@}

@d C TwoTableSums Input
@{
    @<C integer x Input@>
    @<C integer y Input@>
@}

@d C TwoTableSums Answer
@{
    double *PQ_ans
@}

@d C\_TwoTableSums\_dweights\_dsubset
@{
void C_TwoTableSums_dweights_dsubset
(
    @<C TwoTableSums Input@>
    @<C real weights Input@>
    @<C real subset Input@>,
    @<C TwoTableSums Answer@>
)
{
    double *s, *w; 
    @<TwoTableSums Body@>
}
@|C_TwoTableSums_dweights_dsubset
@}

@d C\_TwoTableSums\_iweights\_dsubset
@{
void C_TwoTableSums_iweights_dsubset
(
    @<C TwoTableSums Input@>
    @<C integer weights Input@>
    @<C real subset Input@>,
    @<C TwoTableSums Answer@>
)
{
    double *s;
    int *w; 
    @<TwoTableSums Body@>
}
@|C_TwoTableSums_iweights_dsubset
@}

@d C\_TwoTableSums\_iweights\_isubset
@{
void C_TwoTableSums_iweights_isubset
(
    @<C TwoTableSums Input@>
    @<C integer weights Input@>    
    @<C integer subset Input@>,
    @<C TwoTableSums Answer@>
)
{
    int *s, *w;
    @<TwoTableSums Body@>
}
@|C_TwoTableSums_iweights_isubset
@}

@d C\_TwoTableSums\_dweights\_isubset
@{
void C_TwoTableSums_dweights_isubset
(
    @<C TwoTableSums Input@>
    @<C real weights Input@>
    @<C integer subset Input@>,
    @<C TwoTableSums Answer@>
)
{
    int *s; 
    double *w;
    @<TwoTableSums Body@>
}
@|C_TwoTableSums_dweights_isubset
@}

@d TwoTableSums Body
@{
    int *xx, *yy;

    for (int p = 0; p < Q * P; p++) PQ_ans[p] = 0.0;

    yy = y;
    xx = x;
    @<init subset loop@>
    @<start subset loop@>
    {
        xx = xx + diff;
        yy = yy + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQ_ans[yy[0] * P + xx[0]] += (double) w[0];
        } else {
            PQ_ans[yy[0] * P + xx[0]]++;
        }
        @<continue subset loop@>
    }
    xx = xx + diff;
    yy = yy + diff;
    if (HAS_WEIGHTS) {
        w = w + diff;
        PQ_ans[yy[0] * P + xx[0]] += w[0];
    } else {
        PQ_ans[yy[0] * P + xx[0]]++;
    }
@}

\subsection{ThreeTable Sums}

<<ThreeTableSum>>=

a0 <- c(xtabs(weights ~ ix + iy + block, subset = subset))
a1 <- .Call("R_ThreeTableSums", ix, P + 1L, iy, Q + 1L, 
            block, L, weights, subset - 1L)
a1 <- c(array(a1, dim = c(P + 1, Q + 1, L))[-1,-1,])
a2 <- .Call("R_ThreeTableSums", ix, P + 1L, iy, Q + 1L, 
            block, L,
            as.double(weights), as.double(subset - 1L))
a2 <- c(array(a2, dim = c(P + 1, Q + 1, L))[-1,-1,])
a3 <- .Call("R_ThreeTableSums", ix, P + 1L, iy, Q + 1L, 
            block, L,
            weights, as.double(subset - 1L))
a3 <- c(array(a3, dim = c(P + 1, Q + 1, L))[-1,-1,])
a4 <- .Call("R_ThreeTableSums", ix, P + 1L, iy, Q + 1L, 
            block, L,
            as.double(weights), subset - 1L)
a4 <- c(array(a4, dim = c(P + 1, Q + 1, L))[-1,-1,])

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@



@d R\_ThreeTableSums
@{
SEXP R_ThreeTableSums
(
    @<R x Input@>
    SEXP P,
    @<R y Input@>
    SEXP Q,
    @<R block Input@>
    SEXP L,
    @<R weights Input@>
    @<R subset Input@>
) {

    SEXP ans;
    R_xlen_t N, Nsubset;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * INTEGER(Q)[0] * INTEGER(L)[0]));
    RC_ThreeTableSums(INTEGER(x), N, INTEGER(P)[0], INTEGER(y), INTEGER(Q)[0], 
                      INTEGER(block), INTEGER(L)[0],
                      weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}


@d RC\_ThreeTableSums Prototype
@{
void RC_ThreeTableSums
(
    @<C ThreeTableSums Input@>
    @<R weights Input@>
    @<R subset Input@>,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<C ThreeTableSums Answer@>
) 
@}

@d RC\_ThreeTableSums
@{
@<RC\_ThreeTableSums Prototype@>
{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_ThreeTableSums_iweights_isubset(x, N, P, y, Q, block, L, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQL_ans);
        } else {
            C_ThreeTableSums_iweights_dsubset(x, N, P, y, Q, block, L,
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQL_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_ThreeTableSums_dweights_isubset(x, N, P, y, Q, block, L,
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQL_ans);
        } else {
            C_ThreeTableSums_dweights_dsubset(x, N, P, y, Q, block, L,
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQL_ans);
        }
    }
}
@|RC_ThreeTableSums
@}

@d C ThreeTableSums Input
@{
    @<C integer x Input@>
    @<C integer y Input@>
    @<C integer block Input@>
@}

@d C ThreeTableSums Answer
@{
    double *PQL_ans
@}

@d C\_ThreeTableSums\_dweights\_dsubset
@{
void C_ThreeTableSums_dweights_dsubset
(
    @<C ThreeTableSums Input@>
    @<C real weights Input@>
    @<C real subset Input@>,
    @<C ThreeTableSums Answer@>
)
{
    double *s, *w; 
    @<ThreeTableSums Body@>
}
@|C_ThreeTableSums_dweights_dsubset
@}

@d C\_ThreeTableSums\_iweights\_dsubset
@{
void C_ThreeTableSums_iweights_dsubset
(
    @<C ThreeTableSums Input@>
    @<C integer weights Input@>
    @<C real subset Input@>,
    @<C ThreeTableSums Answer@>
)
{
    double *s;
    int *w; 
    @<ThreeTableSums Body@>
}
@|C_ThreeTableSums_iweights_dsubset
@}

@d C\_ThreeTableSums\_iweights\_isubset
@{
void C_ThreeTableSums_iweights_isubset
(
    @<C ThreeTableSums Input@>
    @<C integer weights Input@>    
    @<C integer subset Input@>,
    @<C ThreeTableSums Answer@>
)
{
    int *s, *w;
    @<ThreeTableSums Body@>
}
@|C_ThreeTableSums_iweights_isubset
@}

@d C\_ThreeTableSums\_dweights\_isubset
@{
void C_ThreeTableSums_dweights_isubset
(
    @<C ThreeTableSums Input@>
    @<C real weights Input@>
    @<C integer subset Input@>,
    @<C ThreeTableSums Answer@>
)
{
    int *s; 
    double *w;
    @<ThreeTableSums Body@>
}
@|C_ThreeTableSums_dweights_isubset
@}

@d ThreeTableSums Body
@{
    int *xx, *yy, *bb, PQ = P * Q;

    for (int p = 0; p < PQ * L; p++) PQL_ans[p] = 0.0;

    yy = y;
    xx = x;
    bb = block;
    @<init subset loop@>
    @<start subset loop@>
    {
        xx = xx + diff;
        yy = yy + diff;
        bb = bb + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += (double) w[0];
        } else {
            PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
        }
        @<continue subset loop@>
    }
    xx = xx + diff;
    yy = yy + diff;
    bb = bb + diff;
    if (HAS_WEIGHTS) {
        w = w + diff;
        PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]] += w[0];
    } else {
        PQL_ans[(bb[0] - 1) * PQ + yy[0] * P + xx[0]]++;
    }
@}

\subsection{Helpers}

<<Helpers>>=
a0 <- as.vector(do.call("c", tapply(subset, block[subset], function(s) s - 1)))
a1 <- .Call("R_order_subset_wrt_block", y, subset - 1L, integer(0), block, L + 1L);
stopifnot(all.equal(a0, a1))
a2 <- .Call("R_order_subset_wrt_block", y, subset - 1L, integer(0), integer(0), 2L);
stopifnot(all.equal(subset - 1L, a2))
a0 <- as.vector(do.call("c", tapply(1:N, block, function(s) s - 1)))
a3 <- .Call("R_order_subset_wrt_block", y, integer(0), integer(0), block, L + 1L);
stopifnot(all.equal(a0, a3))
a4 <- .Call("R_order_subset_wrt_block", y, integer(0), integer(0), integer(0), 2L);
stopifnot(all.equal(0:(N - 1), a4))
@@

@d R\_order\_subset\_wrt\_block
@{
SEXP R_order_subset_wrt_block
(
    @<R y Input@>
    @<R subset Input@>,
    @<R weights Input@>
    @<R block Input@>
    SEXP Nlevels
)
{

    R_xlen_t N = XLENGTH(y) / NCOL(y);
    SEXP blockTable, ans;

    if (XLENGTH(weights) > 0)
        error("cannot deal with weights here");

    if (INTEGER(Nlevels)[0] > 2) {
        PROTECT(blockTable = R_OneTableSums(block, Nlevels, weights, subset));
    } else {
        PROTECT(blockTable = allocVector(REALSXP, 2));
        REAL(blockTable)[0] = 0.0;
        REAL(blockTable)[1] = RC_Sums(N, weights, subset, Offset0, XLENGTH(subset));
    }
    
    PROTECT(ans = RC_order_subset_wrt_block(N, subset, block, blockTable));

    UNPROTECT(2);
    return(ans);
}
@}


@d RC\_order\_subset\_wrt\_block Prototype
@{
SEXP RC_order_subset_wrt_block
(
    @<C integer N Input@>
    @<R subset Input@>,
    @<R block Input@>
    @<R blockTable Input@>
)
@}

@d RC\_order\_subset\_wrt\_block
@{
@<RC\_order\_subset\_wrt\_block Prototype@>
{
    SEXP ans;
    int NOBLOCK = (XLENGTH(block) == 0 || XLENGTH(blockTable) == 2);

    if (XLENGTH(subset) > 0) {
        if (NOBLOCK) {
            return(subset);
        } else {
            PROTECT(ans = allocVector(TYPEOF(subset), XLENGTH(subset)));
            C_order_subset_wrt_block(subset, block, blockTable, ans);
            UNPROTECT(1);
            return(ans);
        }
    } else {
        PROTECT(ans = allocVector(TYPEOF(subset), N));
        if (NOBLOCK) {
            C_setup_subset(N, ans);
        } else {
            C_setup_subset_block(N, block, blockTable, ans);
        }
        UNPROTECT(1);
        return(ans);
    }
}
@}

@d C\_setup\_subset
@{
void C_setup_subset
(
    @<C integer N Input@>
    SEXP ans
)
{
    for (R_xlen_t i = 0; i < N; i++) {
        if (TYPEOF(ans) == INTSXP) {
            INTEGER(ans)[i] = i;
        } else {
            REAL(ans)[i] = (double) i;
        }
    }
}
@}

@d C\_setup\_subset\_block
@{
void C_setup_subset_block
(
    @<C integer N Input@>
    @<R block Input@>
    @<R blockTable Input@>,
    SEXP ans
)
{
    double *cumtable;
    int Nlevels = LENGTH(blockTable);

    cumtable = Calloc(Nlevels, double);
    for (int k = 0; k < Nlevels; k++) cumtable[k] = 0.0;

    /* table[0] are missings, ie block == 0 ! */
    for (int k = 1; k < Nlevels; k++)
        cumtable[k] = cumtable[k - 1] + REAL(blockTable)[k - 1];

    for (R_xlen_t i = 0; i < N; i++) {
        if (TYPEOF(ans) == INTSXP) {
            INTEGER(ans)[(int) cumtable[INTEGER(block)[i]]++] = i;
        } else {
            REAL(ans)[(R_xlen_t) cumtable[INTEGER(block)[i]]++] = (double) i;
        }
    }

    Free(cumtable);
}
@}

@d C\_order\_subset\_wrt\_block
@{
void C_order_subset_wrt_block
(
    @<R subset Input@>,
    @<R block Input@>
    @<R blockTable Input@>,
    SEXP ans
)
{
    double *cumtable;    
    int Nlevels = LENGTH(blockTable);

    cumtable = Calloc(Nlevels, double);
    for (int k = 0; k < Nlevels; k++) cumtable[k] = 0.0;

    /* table[0] are missings, ie block == 0 ! */
    for (int k = 1; k < Nlevels; k++)
        cumtable[k] = cumtable[k - 1] + REAL(blockTable)[k - 1];

    if (TYPEOF(subset) == INTSXP) {
        for (R_xlen_t i = 0; i < XLENGTH(subset); i++)
            INTEGER(ans)[(int) cumtable[INTEGER(block)[INTEGER(subset)[i]]]++] = INTEGER(subset)[i];
    } else {
        for (R_xlen_t i = 0; i < XLENGTH(subset); i++)
            REAL(ans)[(R_xlen_t) cumtable[INTEGER(block)[(R_xlen_t) REAL(subset)[i]]]++] = REAL(subset)[i];
    }

    Free(cumtable); 
}
@}


\section{R Code}

@s


\section{Example}

<<ex>>=
summary(1)
@@

@S

\end{document}
