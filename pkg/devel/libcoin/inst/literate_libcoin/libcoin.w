\documentclass[a4paper]{report}

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

\section{Header and Source Files}

@o libcoin_internal.h -cc
@{
@<R Includes@>
@<C Macros@>
@<C Global Variables@>
@}

@d R Includes
@{
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/stats_package.h> /* for S_rcont2 */
#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Lapack.h> /* for dgesdd */
@}

@d C Macros
@{
#define S(i, j, n) ((i) >= (j) ? (n) * (j) + (i) - (j) * ((j) + 1) / 2 : (n) * (i) + (j) - (i) * ((i) + 1) / 2)
#define LE(x, y, tol)  ((x) < (y)) || (fabs((x) - (y)) < (tol))
#define GE(x, y, tol)  ((x) > (y)) || (fabs((x) - (y)) < (tol))
@| S LE GE
@}

@d C Global Variables
@{
#define ALTERNATIVE_twosided            1
#define ALTERNATIVE_less                2
#define ALTERNATIVE_greater             3

#define TESTSTAT_maximum                1
#define TESTSTAT_quadratic              2

#define LinearStatistic_SLOT            0
#define Expectation_SLOT                1
#define Covariance_SLOT                 2
#define Variance_SLOT                   3
#define MPinv_SLOT                      4
#define ExpectationX_SLOT               5
#define varonly_SLOT                    6
#define dim_SLOT                        7
#define ExpectationInfluence_SLOT       8
#define CovarianceInfluence_SLOT        9
#define VarianceInfluence_SLOT          10
#define Xfactor_SLOT                    11
#define Work_SLOT                       12
#define tol_SLOT                        13
#define PermutedLinearStatistic_SLOT    14
#define TableBlock_SLOT                 15
#define Sumweights_SLOT                 16
#define Table_SLOT                      17

#define DoSymmetric 			1
#define DoCenter 			1
#define DoVarOnly 			1
#define Power1 				1
#define Power2 				2
#define Offset0 			0
@| LinearStatistic_SLOT Expectation_SLOT Covariance_SLOT Variance_SLOT
MPinv_SLOT ExpectationX_SLOT varonly_SLOT dim_SLOT
ExpectationInfluence_SLOT CovarianceInfluence_SLOT VarianceInfluence_SLOT
Xfactor_SLOT Work_SLOT tol_SLOT PermutedLinearStatistic_SLOT
TableBlock_SLOT Sumweights_SLOT Table_SLOT DoSymmetric DoCenter DoVarOnly Power1
Power2 Offset0
@}


The corresponding header file contains definitions of
functions to be used outside \verb|Sums.c|

@o Sums.h -cc
@{
#include "libcoin_internal.h"
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
#include "Sums.h"
#include <R_ext/stats_stubs.h> /* for S_rcont2 */
@<Function Definitions@>
@}

@d Function Definitions
@{
@<MoreUtils@>
@<Memory@>
@<P-Values@>
@<KronSums@>
@<colSums@>
@<SimpleSums@>
@<Tables@>
@<Utils@>
@<LinearStatistics@>
@<Permutations@>
@<ExpectationCovariances@>
@<Test Statistics@>
@<User Interface@>
@<2d User Interface@>
@}

The \proglang{R} interfaces are used to implement
regression tests to be called from within \proglang{R}

\section{Variables}

$N$ is the number of observations

@d R N Input
@{
    SEXP N,
@|N
@}

which at \proglang{C} level is represented as \verb|R_xlen_t| to allow for
$N > $ \verb|INT_MAX|

@d C integer N Input
@{
    R_xlen_t N
@|N
@}

The regressors $\x_i, i = 1, \dots, N$ 

@d R x Input
@{
    SEXP x,
@|x
@}

are either represented as a real matrix with $N$ rows and $P$ columns

@d C integer P Input
@{
    int P
@|P
@}

@d C real x Input
@{
    double *x,
    @<C integer N Input@>,
    @<C integer P Input@>,
@|x
@}

or as a factor (an integer at \proglang{C} level) at $P$ levels

@d C integer x Input
@{
    int *x,
    @<C integer N Input@>,
    @<C integer P Input@>,
@|x
@}

The influence functions are also either a $N \times Q$ real matrix

@d R y Input
@{
    SEXP y,
@|y
@}

@d C integer Q Input
@{
    int Q
@|Q
@}

@d C real y Input
@{
    double *y,
    @<C integer Q Input@>,
@|y
@}

or a factor at $Q$ levels

@d C integer y Input
@{
    int *y,
    @<C integer Q Input@>,
@|y
@}

The weights $w_i, i = 1, \dots, N$ 

@d R weights Input
@{
    SEXP weights
@|weights
@}

can be constant one \verb|XLENGTH(weights) == 0| or integer-valued, with 
\verb|HAS_WEIGHTS == 0| in the former case

@d C integer weights Input
@{
    int *weights,
    int HAS_WEIGHTS,
@|weights, HAS_WEIGHTS
@}

Weights larger than \verb|INT_MAX| are stored as double

@d C real weights Input
@{
    double *weights,
    int HAS_WEIGHTS,
@|weights, HAS_WEIGHTS
@}

The sum of all weights is a double

@d C sumweights Input
@{
    double sumweights
@|sumweights
@}

Subsets $\A \subseteq \{1, \dots, N\}$ are \proglang{R} style indices

@d R subset Input
@{
    SEXP subset
@|subset
@}

are either not existant (\verb|XLENGTH(subset) == 0|) or of length

@d C integer Nsubset Input
@{
    R_xlen_t Nsubset
@|Nsubset
@}

Optionally, one can specify a subset of the subset via

@d C subset range Input
@{
    R_xlen_t offset,
    @<C integer Nsubset Input@>
@|offset
@}

Subsets are stored either as integer

@d C integer subset Input
@{
    int *subset,
    @<C subset range Input@>
@|subset
@}

or double (to allow for indices larger than \verb|INT_MAX|)

@d C real subset Input
@{
    double *subset,
    @<C subset range Input@>
@|subset
@}

Blocks $b_i, i = 1, \dots, N$

@d R block Input
@{
    SEXP block,
@|block
@}

at $B$ levels

@d C integer Lb Input
@{
    int Lb
@|Lb
@}

are stored as a factor

@d C integer block Input
@{
    int *block,
    @<C integer Lb Input@>,
@|block
@}

The tabulation of $b$ (potentially in subsets) is

@d R blockTable Input
@{
    SEXP blockTable
@|blockTable
@}

where the table is of length $B + 1$ and the first element
counts the number of missing values.


<<Sums-setup>>=
### replace with library("libcoin")
dyn.load("Sums.so")
set.seed(29)
N <- 20L
P <- 3L
Lx <- 10L
Ly <- 5L
Q <- 4L
B <- 2L
iX2d <- rbind(0, matrix(runif(Lx * P), nrow = Lx))
ix <- sample(1:Lx, size = N, replace = TRUE)
levels(ix) <- 1:Lx
ixf <- factor(ix, levels = 1:Lx, labels = 1:Lx)
x <- iX2d[ix + 1,]
Xfactor <- diag(Lx)[ix,]
iY2d <- rbind(0, matrix(runif(Ly * Q), nrow = Ly))
iy <- sample(1:Ly, size = N, replace = TRUE)
levels(iy) <- 1:Ly
iyf <- factor(iy, levels = 1:Ly, labels = 1:Ly)
y <- iY2d[iy + 1,]
weights <- sample(0:5, size = N, replace = TRUE)
block <- sample(gl(B, ceiling(N / B))[1:N])
subset <- sort(sample(1:N, floor(N * 1.5), replace = TRUE))
subsety <- sample(1:N, floor(N * 1.5), replace = TRUE)
r1 <- rep(1:ncol(x), ncol(y))
r1Xfactor <- rep(1:ncol(Xfactor), ncol(y))
r2 <- rep(1:ncol(y), each = ncol(x))
r2Xfactor <- rep(1:ncol(y), each = ncol(Xfactor))
@@

<<Rlibcoin>>=
LECV <- function(X, Y, weights = integer(0), subset = integer(0), block = integer(0)) {

    if (length(weights) == 0) weights <- rep(1, nrow(X))
    if (length(subset) == 0) subset <- 1:nrow(X)
    idx <- rep(subset, weights[subset])
    X <- X[idx,,drop = FALSE]
    Y <- Y[idx,,drop = FALSE]
    sumweights <- length(idx)

    if (length(block) == 0) {
        ExpX <- colSums(X)
        ExpY <- colSums(Y) / sumweights
        yc <- t(t(Y) - ExpY)
        CovY <- crossprod(yc) / sumweights
        CovX <- crossprod(X)
        Exp <- kronecker(ExpY, ExpX) 
        Cov <- sumweights / (sumweights - 1) * kronecker(CovY, CovX) - 
               1 / (sumweights - 1) * kronecker(CovY, tcrossprod(ExpX))
 
        ret <- list(LinearStatistic = as.vector(crossprod(X, Y)),
                    Expectation = as.vector(Exp), 
                    Covariance = Cov[lower.tri(Cov, diag = TRUE)],
                    Variance = diag(Cov))
   } else {
        block <- block[idx]
        ret <- list(LinearStatistic = 0, Expectation = 0, Covariance = 0, Variance = 0)
        for (b in levels(block)) {
            tmp <- LECV(X = X, Y = Y, subset = which(block == b))
            for (l in names(ret)) ret[[l]] <- ret[[l]] + tmp[[l]]
        }
   }
   return(ret)
}

cmpr <- function(ret1, ret2) {
    ret1 <- ret1[!sapply(ret1, is.null)]
    ret2 <- ret2[!sapply(ret2, is.null)]
    nm1 <- names(ret1)
    nm2 <- names(ret2)
    nm <- c(nm1, nm2)
    nm <- names(table(nm))[table(nm) == 2]
    all.equal(ret1[nm], ret2[nm])
}

a <- .Call("R_ExpectationCovarianceStatistic", x, y, integer(0), integer(0),
            integer(0), 0L, 0.00001)
b <- LECV(x, y)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
           iY2d, iy, integer(0), integer(0), 
           integer(0), 0L, 0.00001)

cmpr(b, d)

a <- .Call("R_ExpectationCovarianceStatistic", x, y, integer(0), integer(0),
            block, 0L, 0.00001)
b <- LECV(x, y, block = block)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
           iY2d, iy, integer(0), integer(0), 
           block, 0L, 0.00001)

cmpr(b, d)

a <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, integer(0),
            integer(0), 0L, 0.00001)
b <- LECV(x, y, weights = weights)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
           iY2d, iy, weights, integer(0), 
           integer(0), 0L, 0.00001)

cmpr(b, d)

a <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, integer(0),
            block, 0L, 0.00001)
b <- LECV(x, y, weights = weights, block = block)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
           iY2d, iy, weights, integer(0), 
           block, 0L, 0.00001)

cmpr(b, d)


a <- .Call("R_ExpectationCovarianceStatistic", x, y, integer(0), subset,
            integer(0), 0L, 0.00001)
b <- LECV(x, y, subset = subset)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
           iY2d, iy, integer(0),  subset,
           integer(0), 0L, 0.00001)

cmpr(b, d)

a <- .Call("R_ExpectationCovarianceStatistic", x, y, integer(0), subset,
            block, 0L, 0.00001)
b <- LECV(x, y, subset = subset, block = block)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
           iY2d, iy, integer(0),  subset,
           block, 0L, 0.00001)

cmpr(b, d)



a <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset,
            integer(0), 0L, 0.00001)
b <- LECV(x, y, weights = weights, subset = subset)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
           iY2d, iy, weights, subset, 
           integer(0), 0L, 0.00001)

cmpr(b, d)

a <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset,
           block, 0L, 0.00001)
b <- LECV(x, y, weights = weights, subset = subset, block = block)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
           iY2d, iy, weights, subset, 
           block, 0L, 0.00001)

cmpr(b, d)

a <- .Call("R_ExpectationCovarianceStatistic", ix, y, integer(0), integer(0),
            integer(0), 0L, 0.00001)
b <- LECV(Xfactor, y)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
           iY2d, iy, integer(0), integer(0),
           integer(0), 0L, 0.00001)

cmpr(b, d)

a <- .Call("R_ExpectationCovarianceStatistic", ix, y, integer(0), integer(0),
           block, 0L, 0.00001)
b <- LECV(Xfactor, y, block = block)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
           iY2d, iy, integer(0), integer(0),
           block, 0L, 0.00001)

cmpr(b, d)

a <- .Call("R_ExpectationCovarianceStatistic", ix, y, weights, integer(0),
            integer(0), 0L, 0.00001)
b <- LECV(Xfactor, y, weights = weights)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
           iY2d, iy, weights, integer(0),
           integer(0), 0L, 0.00001)

cmpr(b, d)
a <- .Call("R_ExpectationCovarianceStatistic", ix, y, weights, integer(0),
           block, 0L, 0.00001)
b <- LECV(Xfactor, y, weights = weights, block = block)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
           iY2d, iy, weights, integer(0),
           block, 0L, 0.00001)

cmpr(b, d)


a <- .Call("R_ExpectationCovarianceStatistic", ix, y, integer(0), subset,
            integer(0), 0L, 0.00001)
b <- LECV(Xfactor, y, subset = subset)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
           iY2d, iy, integer(0), subset,
           integer(0), 0L, 0.00001)

cmpr(b, d)
a <- .Call("R_ExpectationCovarianceStatistic", ix, y, integer(0), subset,
           block, 0L, 0.00001)
b <- LECV(Xfactor, y, subset = subset, block = block)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
           iY2d, iy, integer(0), subset,
           block, 0L, 0.00001)

cmpr(b, d)


a <- .Call("R_ExpectationCovarianceStatistic", ix, y, weights, subset,
            integer(0), 0L, 0.00001)
b <- LECV(Xfactor, y, weights = weights, subset = subset)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
           iY2d, iy, weights, subset,
           integer(0), 0L, 0.00001)

cmpr(b, d)


a <- .Call("R_ExpectationCovarianceStatistic", ix, y, weights, subset,
           block, 0L, 0.00001)
b <- LECV(Xfactor, y, weights = weights, subset = subset, block = block)
cmpr(a, b)

d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
           iY2d, iy, weights, subset,
           block, 0L, 0.00001)

cmpr(b, d)

@@

\section{User Interface}

\subsection{1d Case}

<<1d>>=
library("libcoin")
LECV <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset, 
               integer(0), 0L, 0.00001)
lcv <- LinStatExpCov(X = x, Y = y, weights = weights, subset = subset)
all.equal(LECV, lcv)

iLECV <- .Call("R_ExpectationCovarianceStatistic", ix, y, weights, subset, 
               integer(0), 0L, 0.00001)
ilcv <- LinStatExpCov(X = ix, Y = y, weights = weights, subset = subset)
all.equal(iLECV, ilcv)

iLECVXfactor <- .Call("R_ExpectationCovarianceStatistic", Xfactor, y, weights, subset, 
               integer(0), 0L, 0.00001)
all.equal(iLECVXfactor, ilcv)
all.equal(iLECVXfactor, iLECV)

V <- matrix(0, nrow = P * Q, ncol = P * Q)
V[lower.tri(V, diag = TRUE)] <- LECV$Covariance
LSvar <- diag(V)
V <- matrix(0, nrow = Q, ncol = Q)
V[lower.tri(V, diag = TRUE)] <- LECV$CovarianceInfluence
Ivar <- diag(V)
(LEV <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset, 
              integer(0), 1L, 0.00001))
lv <- LinStatExpCov(X = x, Y = y, weights = weights, subset = subset, varonly = TRUE)
all.equal(LEV, lv)
stopifnot(all.equal(LECV$LinearStatistic, LEV$LinearStatistic) &&
          all.equal(LECV$Expectation, LEV$Expectation) &&
          all.equal(LECV$ExpectationInfluence, LEV$ExpectationInfluence) &&
          all.equal(Ivar, LEV$VarianceInfluence) &&
          all.equal(LSvar, LEV$Variance))

library("libcoin")
subset[block[subset] == 1]
subset[block[subset] == 2]
LECVb <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset, 
               block, 0L, 0.00001)
lcvb <- LinStatExpCov(X = x, Y = y, weights = weights, subset = subset, block = block)
all.equal(LECVb, lcvb)

LEVb <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset, 
               block, 1L, 0.00001)
lvb <- LinStatExpCov(X = x, Y = y, weights = weights, subset = subset, block = block, varonly =
TRUE)
all.equal(LEVb, lvb)
@@

@d User Interface
@{
@<RC\_ExpectationCovarianceStatistic@>
@<R\_ExpectationCovarianceStatistic@>
@<R\_PermutedLinearStatistic@>
@}

@d User Interface Inputs
@{
@<R x Input@>
@<R y Input@>
@<R weights Input@>,
@<R subset Input@>,
@<R block Input@>
@}

@d R\_ExpectationCovarianceStatistic
@{
SEXP R_ExpectationCovarianceStatistic
(
@<User Interface Inputs@>
SEXP varonly,
SEXP tol
) {
    SEXP ans;

    @<Setup Dimensions@>

    PROTECT(ans = RC_init_LECV_1d(P, Q, INTEGER(varonly)[0], Lb, isInteger(x), REAL(tol)[0]));

    RC_ExpectationCovarianceStatistic(x, y, weights, subset, block, ans);

    UNPROTECT(1);
    return(ans);
}
@|R_ExpectationCovarianceStatistic
@}


@d Setup Dimensions
@{
    int P, Q, Lb;

    if (isInteger(x)) {
        P = NLEVELS(x);
    } else {
        P = NCOL(x);
    }
    Q = NCOL(y);

    Lb = 1;
    if (LENGTH(block) > 0)
        Lb = NLEVELS(block);
@}


@d RC\_ExpectationCovarianceStatistic
@{
void RC_ExpectationCovarianceStatistic
(
@<User Interface Inputs@>
SEXP ans
) {

    @<C integer N Input@>;
    @<C integer P Input@>;
    @<C integer Q Input@>;
    @<C integer Lb Input@>;
    double *sumweights, *table;
    double *ExpInf, *VarInf, *CovInf, *ExpX, *ExpXtotal, *VarX, *CovX;
    SEXP nullvec, subset_block;

    /* note: x being an integer (Xfactor) with some 0 elements is not
             handled correctly (as sumweights doesnt't take this information
             into account; use subset to exclude these missings (as done
             in libcoin::LinStatExpCov) */

    @<Extract Dimensions@>

    @<Compute Linear Statistic@>

    @<Setup Memory and Subsets in Blocks@>

    /* start with subset[0] */
    R_xlen_t offset = (R_xlen_t) table[0];

    for (int b = 0; b < Lb; b++) {

        /* compute sum of weights in block b of subset */
        sumweights[b] = RC_Sums(N, weights, subset_block, 
                                offset, (R_xlen_t) table[b + 1]);

        /* an empty block level; don't do anything */
        if (sumweights[b] == 0) continue;

        @<Compute Expectation Linear Statistic@>

        if (C_get_varonly(ans)) {
            @<Compute Variance Linear Statistic@>
        } else {
            @<Compute Covariance Linear Statistic@>
        }

        /* next iteration starts with subset[table[b + 1]] */
        offset += (R_xlen_t) table[b + 1];
    }

    Free(ExpX); Free(VarX); Free(CovX);
    UNPROTECT(2);
}
@|RC_ExpectationCovarianceStatistic
@}


@d Extract Dimensions
@{
P = C_get_P(ans);
Q = C_get_Q(ans);
N = NROW(x);
Lb = C_get_Lb(ans);
@}

@d Compute Linear Statistic
@{
PROTECT(nullvec = allocVector(INTSXP, 0));
RC_LinearStatistic(x, N, P, REAL(y), Q, weights, subset, 
                   Offset0, XLENGTH(subset), nullvec,
                   C_get_LinearStatistic(ans));
@}

@d Setup Memory and Subsets in Blocks
@{
ExpInf = C_get_ExpectationInfluence(ans);
VarInf = C_get_VarianceInfluence(ans);
CovInf = C_get_CovarianceInfluence(ans);
ExpXtotal = C_get_ExpectationX(ans);
for (int p = 0; p < P; p++) ExpXtotal[p] = 0.0;
ExpX = Calloc(P, double);
VarX = Calloc(P, double);
CovX = Calloc(P * (P + 1) / 2, double);
table = C_get_TableBlock(ans);
sumweights = C_get_Sumweights(ans);

if (Lb == 1) {
    table[0] = 0.0;
    table[1] = RC_Sums(N, nullvec, subset, Offset0, XLENGTH(subset));
} else {
    RC_OneTableSums(INTEGER(block), N, Lb + 1, nullvec, subset, Offset0, 
                    XLENGTH(subset), table);
}
if (table[0] > 0)
    error("No missing values allowed in block");
PROTECT(subset_block = RC_order_subset_wrt_block(N, subset, block, 
                                                 VECTOR_ELT(ans, TableBlock_SLOT)));
@}

@d Compute Expectation Linear Statistic
@{
RC_ExpectationInfluence(N, y, Q, weights, subset_block, offset, 
                        (R_xlen_t) table[b + 1], sumweights[b], ExpInf + b * Q);
RC_ExpectationX(x, N, P, weights, subset_block, offset, 
                (R_xlen_t) table[b + 1], ExpX);
for (int p = 0; p < P; p++) ExpXtotal[p] += ExpX[p];
C_ExpectationLinearStatistic(P, Q, ExpInf + b * Q, ExpX, b, 
                             C_get_Expectation(ans));
@}

@d Compute Variance Linear Statistic
@{
RC_CovarianceInfluence(N, y, Q, weights, subset_block, offset, 
                       (R_xlen_t) table[b + 1], ExpInf + b * Q, sumweights[b], 
                       DoVarOnly, VarInf + b * Q);
RC_CovarianceX(x, N, P, weights, subset_block, offset, 
               (R_xlen_t) table[b + 1], ExpX, DoVarOnly, VarX);
C_VarianceLinearStatistic(P, Q, VarInf + b * Q, ExpX, VarX, sumweights[b], 
                          b, C_get_Variance(ans));
@}

@d Compute Covariance Linear Statistic
@{
RC_CovarianceInfluence(N, y, Q, weights, subset_block, (R_xlen_t) table[b], 
                       (R_xlen_t) table[b + 1], ExpInf + b * Q, sumweights[b], 
                       !DoVarOnly, CovInf + b * Q * (Q + 1) / 2);
RC_CovarianceX(x, N, P, weights, subset_block, (R_xlen_t) table[b], 
               (R_xlen_t) table[b + 1], ExpX, !DoVarOnly, CovX);
C_CovarianceLinearStatistic(P, Q, CovInf + b * Q * (Q + 1) / 2,
                            ExpX, CovX, sumweights[b], b,
                            C_get_Covariance(ans));
@}

<<permutations>>=
a <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset,
           integer(0), 0L, 0.00001)
.Call("R_PermutedLinearStatistic", x, y, weights, subset, integer(0), 10, list())
.Call("R_PermutedLinearStatistic", x, y, weights, subset, integer(0), 10, a)
.Call("R_PermutedLinearStatistic", x, y, weights, subset, block, 10, a)
@@

@d R\_PermutedLinearStatistic Prototype
@{
SEXP R_PermutedLinearStatistic
(
    @<User Interface Inputs@>
    SEXP nperm,
    @<R LECV Input@>
)
@}

@d R\_PermutedLinearStatistic
@{
@<R\_PermutedLinearStatistic Prototype@>
{
    SEXP ans, expand_subset, block_subset, perm, tmp, blockTable;
    double *linstat;
    int PQ;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;
    R_xlen_t inperm;

    @<Setup Dimensions@>
    PQ = P * Q;
    N = NROW(y);
    inperm = (R_xlen_t) REAL(nperm)[0];

    PROTECT(ans = allocMatrix(REALSXP, PQ, inperm));
    PROTECT(expand_subset = RC_setup_subset(N, weights, subset));
    Nsubset = XLENGTH(expand_subset);
    PROTECT(tmp = allocVector(REALSXP, Nsubset));
    PROTECT(perm = allocVector(REALSXP, Nsubset));

    GetRNGstate();
    if (Lb == 1) {
        for (R_xlen_t np = 0; np < inperm; np++) {
            @<Setup Linear Statistic@>
            C_doPermute(REAL(expand_subset), Nsubset, REAL(tmp), REAL(perm));
            @<Compute Permuted Linear Statistic@>
        }
    } else {
        PROTECT(blockTable = allocVector(REALSXP, Lb + 1));
        /* same as RC_OneTableSums(block, noweights, expand_subset) */
        RC_OneTableSums(INTEGER(block), XLENGTH(block), Lb + 1, weights, subset, Offset0,
                        XLENGTH(subset), REAL(blockTable));
        PROTECT(block_subset = RC_order_subset_wrt_block(XLENGTH(block), expand_subset, 
                                                         block, blockTable));
        for (R_xlen_t np = 0; np < inperm; np++) {
            @<Setup Linear Statistic@>
            C_doPermuteBlock(REAL(block_subset), Nsubset, REAL(blockTable), 
                             Lb + 1, REAL(tmp), REAL(perm));
            @<Compute Permuted Linear Statistic@>
        }
        UNPROTECT(2);
    }
    PutRNGstate();

    @<Standardise Linear Statistics@>

    UNPROTECT(4);
    return(ans);
}
@|R_PermutedLinearStatistic
@}

@d Setup Linear Statistic
@{
if (np % 256 == 0) R_CheckUserInterrupt();
linstat = REAL(ans) + PQ * np;
for (int p = 0; p < PQ; p++)
    linstat[p] = 0.0;
@}

@d Compute Permuted Linear Statistic
@{
/* does not require weights (as RC_LinearStatistic) */
RC_KronSums_Permutation(x, NROW(x), P, REAL(y), Q,
                        expand_subset, Offset0, Nsubset,
                        perm, linstat);
@}

@d Standardise Linear Statistics
@{
if (LENGTH(LECV) > 0) {
    for (R_xlen_t np = 0; np < inperm; np++) {
        if (C_get_varonly(LECV)) {
            C_standardise(PQ, REAL(ans) + PQ * np, C_get_Expectation(LECV),
                          C_get_Variance(LECV), 1, C_get_tol(LECV));
        } else {
            C_standardise(PQ, REAL(ans) + PQ * np, C_get_Expectation(LECV),
                          C_get_Covariance(LECV), 0, C_get_tol(LECV));
        }
    }
}
@}

\subsection{2d Case}

<<2d>>=
library("libcoin")
LECV2d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
              iY2d, iy, weights, subset, 
              integer(0), 0L, 0.00001)
lcv2d <- LinStatExpCov(X = iX2d, ix = ix, Y = iY2d, iy = iy, weights = weights, subset = subset)
all.equal(LECV2d, lcv2d)
all.equal(LECV2d, LECV)

LECV2d <- .Call("R_ExpectationCovarianceStatistic_2d", integer(0), ix, 
              iY2d, iy, weights, subset, 
              integer(0), 0L, 0.00001)
lcv2d <- LinStatExpCov(X = integer(0), ix = ix, Y = iY2d, iy = iy, weights = weights, subset = subset)
all.equal(LECV2d, lcv2d)

LECVXfactor <- .Call("R_ExpectationCovarianceStatistic", Xfactor, y,
              weights, subset, 
              integer(0), 0L, 0.00001)
all.equal(LECV2d, LECVXfactor)

@@

@d 2d User Interface
@{
@<RC\_ExpectationCovarianceStatistic\_2d@>
@<R\_ExpectationCovarianceStatistic\_2d@>
@<R\_PermutedLinearStatistic\_2d@>
@}

@d 2d User Interface Inputs
@{
@<R x Input@>
SEXP ix,
@<R y Input@>
SEXP iy,
@<R weights Input@>,
@<R subset Input@>,
@<R block Input@>
@}

@d R\_ExpectationCovarianceStatistic\_2d
@{
SEXP R_ExpectationCovarianceStatistic_2d
(
@<2d User Interface Inputs@>
SEXP varonly,
SEXP tol
) {
    SEXP ans;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;
    int Xfactor;

    N = XLENGTH(ix);
    Nsubset = XLENGTH(subset);
    Xfactor = XLENGTH(x) == 0;

    @<Setup Dimensions 2d@>

    PROTECT(ans = RC_init_LECV_2d(P, Q, INTEGER(varonly)[0], 
                                  Lx, Ly, Lb, Xfactor, REAL(tol)[0]));

    if (Lb == 1) {
        RC_TwoTableSums(INTEGER(ix), N, Lx + 1, INTEGER(iy), Ly + 1, 
                        weights, subset, Offset0, Nsubset, 
                        C_get_Table(ans));
    } else {
        RC_ThreeTableSums(INTEGER(ix), N, Lx + 1, INTEGER(iy), Ly + 1, 
                          INTEGER(block), Lb, weights, subset, Offset0, Nsubset, 
                          C_get_Table(ans));
    }
    RC_ExpectationCovarianceStatistic_2d(x, ix, y, iy, weights,
                                         subset, block, ans);

    UNPROTECT(1);
    return(ans);
}
@|R_ExpectationCovarianceStatistic_2d
@}

@d Setup Dimensions 2d
@{
int P, Q, Lb, Lx, Ly;

if (XLENGTH(x) == 0) {
    P = NLEVELS(ix);
} else {
    P = NCOL(x);
}
Q = NCOL(y);

Lb = 1;
if (LENGTH(block) > 0)
    Lb = NLEVELS(block);

Lx = NLEVELS(ix);
Ly = NLEVELS(iy);
@}

@d Linear Statistic 2d
@{
if (Xfactor) {
    for (int j = 1; j < Lyp1; j++) { /* j = 0 means NA */
        for (int i = 1; i < Lxp1; i++) { /* i = 0 means NA */
            for (int q = 0; q < Q; q++)
                linstat[q * (Lxp1 - 1) + (i - 1)] +=
                    btab[j * Lxp1 + i] * REAL(y)[q * Lyp1 + j];
        }
    }
} else {
    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++) {
            int qPp = q * P + p;
            int qLy = q * Lyp1;
            for (int i = 0; i < Lxp1; i++) {
                int pLxi = p * Lxp1 + i;
                for (int j = 0; j < Lyp1; j++)
                    linstat[qPp] += REAL(y)[qLy + j] * REAL(x)[pLxi] * btab[j * Lxp1 + i];
            }
        }
    }
}
@}

@d 2d Total Table
@{
for (int i = 0; i < Lxp1 * Lyp1; i++)
    table2d[i] = 0.0;
for (int b = 0; b < Lb; b++) {
    for (int i = 0; i < Lxp1; i++) {
        for (int j = 0; j < Lyp1; j++)
            table2d[j * Lxp1 + i] += table[b * Lxp1 * Lyp1 + j * Lxp1 + i];
    }
}
@}

@d Col Row Total Sums
@{
/* column sums */
for (int q = 0; q < Lyp1; q++) {
    csum[q] = 0;
    for (int p = 0; p < Lxp1; p++)
        csum[q] += btab[q * Lxp1 + p];
}
csum[0] = 0; /* NA */
/* row sums */

for (int p = 0; p < Lxp1; p++)  {
    rsum[p] = 0;
    for (int q = 0; q < Lyp1; q++)
        rsum[p] += btab[q * Lxp1 + p];
}
rsum[0] = 0; /* NA */
/* total sum */
sumweights[b] = 0;
for (int i = 1; i < Lxp1; i++) sumweights[b] += rsum[i];
@}

@d 2d Expectation
@{
RC_ExpectationInfluence(NROW(y), y, Q, Rcsum, subset, Offset0, 0, sumweights[b], ExpInf);

if (LENGTH(x) == 0) {
    for (int p = 0; p < P; p++)
        ExpX[p] = rsum[p + 1];
    } else {
        RC_ExpectationX(x, NROW(x), P, Rrsum, subset, Offset0, 0, ExpX);
}

C_ExpectationLinearStatistic(P, Q, ExpInf, ExpX, b, C_get_Expectation(ans));
@}

@d 2d Covariance
@{
if (C_get_varonly(ans)) {
    RC_CovarianceInfluence(NROW(y), y, Q, Rcsum, subset, Offset0, 0, ExpInf, sumweights[b],
                           DoVarOnly, C_get_VarianceInfluence(ans));
} else {
    RC_CovarianceInfluence(NROW(y), y, Q, Rcsum, subset, Offset0, 0, ExpInf, sumweights[b],
                           !DoVarOnly, C_get_CovarianceInfluence(ans));
}

if (C_get_varonly(ans)) {
    if (LENGTH(x) == 0) {
        for (int p = 0; p < P; p++) CovX[p] = ExpX[p];
    } else {
        RC_CovarianceX(x, NROW(x), P, Rrsum, subset, Offset0, 0, ExpX, DoVarOnly, CovX);
    }
    C_VarianceLinearStatistic(P, Q, C_get_VarianceInfluence(ans),
                              ExpX, CovX, sumweights[b], b,
                              C_get_Variance(ans));
} else {
    if (LENGTH(x) == 0) {
        for (int p = 0; p < P * (P + 1) / 2; p++) CovX[p] = 0.0;
        for (int p = 0; p < P; p++) CovX[S(p, p, P)] = ExpX[p];
    } else {
        RC_CovarianceX(x, NROW(x), P, Rrsum, subset, Offset0, 0, ExpX, !DoVarOnly, CovX);
    }
    C_CovarianceLinearStatistic(P, Q, C_get_CovarianceInfluence(ans),
                                ExpX, CovX, sumweights[b], b,
                                C_get_Covariance(ans));
}
@}

@d RC\_ExpectationCovarianceStatistic\_2d
@{
void RC_ExpectationCovarianceStatistic_2d
(
@<2d User Interface Inputs@>
SEXP ans
) {

    SEXP Rcsum, Rrsum;
    int P, Q, Lxp1, Lyp1, Lb, Xfactor;
    double *ExpInf, *ExpX, *CovX;
    double *table, *table2d, *csum, *rsum, *sumweights, *btab, *linstat;

    P = C_get_P(ans);
    Q = C_get_Q(ans);

    ExpInf = C_get_ExpectationInfluence(ans);
    ExpX = C_get_ExpectationX(ans);
    table = C_get_Table(ans);
    sumweights = C_get_Sumweights(ans);

    Lxp1 = C_get_dimTable(ans)[0];
    Lyp1 = C_get_dimTable(ans)[1];
    Lb = C_get_Lb(ans);
    Xfactor = C_get_Xfactor(ans);

    if (C_get_varonly(ans)) {
        CovX = Calloc(P, double);
    } else {
        CovX = Calloc(P * (P + 1) / 2, double);
    }

    table2d = Calloc(Lxp1 * Lyp1, double);
    PROTECT(Rcsum = allocVector(REALSXP, Lyp1));
    csum = REAL(Rcsum);
    PROTECT(Rrsum = allocVector(REALSXP, Lxp1));
    rsum = REAL(Rrsum);

    @<2d Total Table@>

    linstat = C_get_LinearStatistic(ans);
    for (int p = 0; p < P * Q; p++)
        linstat[p] = 0.0;

    for (int b = 0; b < Lb; b++) {
        btab = table + Lxp1 * Lyp1 * b;

        @<Linear Statistic 2d@>

        @<Col Row Total Sums@>

        @<2d Expectation@>

        @<2d Covariance@>

    }
    Free(table2d); 
    UNPROTECT(2);
}
@|RC_ExpectationCovarianceStatistic
@}

<<permutations-2d>>=
LECV2d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
              iY2d, iy, weights, subset, 
              integer(0), 0L, 0.00001)
LECV2d$Table
.Call("R_PermutedLinearStatistic_2d", iX2d, ix, 
              iY2d, iy, weights, subset, 
              integer(0), 10, LECV2d$Table, list())
LECV2d <- .Call("R_ExpectationCovarianceStatistic_2d", iX2d, ix, 
              iY2d, iy, weights, subset, 
              block, 0L, 0.00001)
LECV2d$Table
.Call("R_PermutedLinearStatistic_2d", iX2d, ix, 
              iY2d, iy, weights, subset, 
              block, 10, LECV2d$Table, list())
@@

@d R\_PermutedLinearStatistic\_2d Prototype
@{
SEXP R_PermutedLinearStatistic_2d
(
    @<2d User Interface Inputs@>
    SEXP nperm,
    SEXP itable,
    @<R LECV Input@>
)
@}

@d R\_PermutedLinearStatistic\_2d
@{
@<R\_PermutedLinearStatistic\_2d Prototype@>
{
    SEXP ans, Ritable;
    int *csum, *rsum, *sumweights, *jwork, *table, *rtable2, maxn = 0, Lxp1, Lyp1, *btab, PQ, Xfactor;
    R_xlen_t inperm;
    double *fact, *linstat;

    @<Setup Dimensions 2d@>

    PQ = P * Q;
    Xfactor = XLENGTH(x) == 0;
    Lxp1 = Lx + 1;
    Lyp1 = Ly + 1;
    inperm = (R_xlen_t) REAL(nperm)[0];

    PROTECT(ans = allocMatrix(REALSXP, PQ, inperm));

    @<Setup Working Memory@>

    @<Convert Table to Integer@>

    for (int b = 0; b < Lb; b++) {
        btab = INTEGER(Ritable) + Lxp1 * Lyp1 * b;
        @<Col Row Total Sums@>
        if (sumweights[b] > maxn) maxn = sumweights[b];
    }

    @<Setup Log-Factorials@>

    GetRNGstate();

    for (R_xlen_t np = 0; np < inperm; np++) {

        @<Setup Linear Statistic@>

        for (int p = 0; p < Lxp1 * Lyp1; p++)
            table[p] = 0;

        for (int b = 0; b < Lb; b++) {
            @<Compute Permuted Linear Statistic 2d@>
        }
    }

    PutRNGstate();

    @<Standardise Linear Statistics@>

    Free(csum); Free(rsum); Free(sumweights); Free(rtable2);
    Free(jwork); Free(fact);
    UNPROTECT(2);
    return(ans);
}
@|R_PermutedLinearStatistic_2d
@}

@d Convert Table to Integer
@{
PROTECT(Ritable = allocVector(INTSXP, LENGTH(itable)));
for (int i = 0; i < LENGTH(itable); i++) {
    if (REAL(itable)[i] > INT_MAX)
        error("cannot deal with weights larger INT_MAX in R_PermutedLinearStatistic_2d");
    INTEGER(Ritable)[i] = (int) REAL(itable)[i];
}
@}

@d Setup Working Memory
@{
csum = Calloc(Lyp1 * Lb, int);
rsum = Calloc(Lxp1 * Lb, int);
sumweights = Calloc(Lb, int);
table = Calloc(Lxp1 * Lyp1, int);
rtable2 = Calloc(Lx * Ly , int);
jwork = Calloc(Lyp1, int);
@}

@d Setup Log-Factorials
@{
fact = Calloc(maxn + 1, double);
/* Calculate log-factorials.  fact[i] = lgamma(i+1) */
fact[0] = fact[1] = 0.;
for(int j = 2; j <= maxn; j++)
    fact[j] = fact[j - 1] + log(j);
@}

@d Compute Permuted Linear Statistic 2d
@{
S_rcont2(&Lx, &Ly, rsum + Lxp1 * b + 1,
         csum + Lyp1 *b + 1, sumweights + b, fact, jwork, rtable2);

for (int j1 = 1; j1 <= Lx; j1++) {
    for (int j2 = 1; j2 <= Ly; j2++)
        table[j2 * Lxp1 + j1] = rtable2[(j2 - 1) * Lx + (j1 - 1)];
}
btab = table;
@<Linear Statistic 2d@>
@}


\section{Test Statistics}

@d Test Statistics
@{
@<C\_maxstand\_Covariance@>
@<C\_maxstand\_Variance@>
@<C\_minstand\_Covariance@>
@<C\_minstand\_Variance@>
@<C\_maxabsstand\_Covariance@>
@<C\_maxabsstand\_Variance@>
@<C\_quadform@>
@<C\_maxtype@>
@<C\_standardise@>
@<C\_ordered\_Xfactor@>
@<C\_unordered\_Xfactor@>
@}

@d C\_maxstand\_Covariance
@{
double C_maxstand_Covariance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar_sym,
    const double tol
) {

    double ans = R_NegInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(covar_sym[S(p, p, PQ)]);
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}
@|C_maxstand_Covariance
@}

@d C\_maxstand\_Variance
@{
double C_maxstand_Variance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *var,
    const double tol
) {

    double ans = R_NegInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(var[p]);
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}
@|C_maxstand_Variance
@}

@d C\_minstand\_Covariance
@{
double C_minstand_Covariance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar_sym,
    const double tol
) {

    double ans = R_PosInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(covar_sym[S(p, p, PQ)]);
        if (tmp < ans) ans = tmp;
    }
    return(ans);
}
@|C_minstand_Covariance
@}

@d C\_minstand\_Variance
@{
double C_minstand_Variance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *var,
    const double tol
) {

    double ans = R_PosInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(var[p]);
        if (tmp < ans) ans = tmp;
    }
    return(ans);
}
@|C_minstand_Variance
@}

@d C\_maxabsstand\_Covariance
@{
double C_maxabsstand_Covariance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar_sym,
    const double tol
) {

    double ans = R_NegInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = fabs((linstat[p] - expect[p]) /
                  sqrt(covar_sym[S(p, p, PQ)]));
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}
@|C_maxabsstand_Covariance
@}

@d C\_maxabsstand\_Variance
@{
double C_maxabsstand_Variance
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *var,
    const double tol
) {

    double ans = R_NegInf, tmp = 0.0;

    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = fabs((linstat[p] - expect[p]) / sqrt(var[p]));
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}
@|C_maxabsstand_Variance
@}

@d C\_quadform
@{
double C_quadform
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *MPinv_sym
) {

    double ans = 0.0, tmp = 0.0;

    for (int q = 0; q < PQ; q++) {
        tmp = 0.0;
        for (int p = 0; p < PQ; p++)
            tmp += (linstat[p] - expect[p]) * MPinv_sym[S(p, q, PQ)];
        ans += tmp * (linstat[q] - expect[q]);
    }
    return(ans);
}
@|C_quadform
@}

@d C\_maxtype
@{
double C_maxtype
(
    const int PQ,
    const double *linstat,
    const double *expect,
    const double *covar,
    const int varonly,
    const double tol,
    const int alternative
) {

    double ret = 0.0;

    if (varonly) {
        if (alternative ==  ALTERNATIVE_twosided) {
            ret = C_maxabsstand_Variance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_less) {
            ret = C_minstand_Variance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_greater) {
            ret = C_maxstand_Variance(PQ, linstat, expect, covar, tol);
        }
    } else {
        if (alternative ==  ALTERNATIVE_twosided) {
            ret = C_maxabsstand_Covariance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_less) {
            ret = C_minstand_Covariance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_greater) {
            ret = C_maxstand_Covariance(PQ, linstat, expect, covar, tol);
        }
    }
    return(ret);
}
@|C_maxtype
@}

@d C\_standardise
@{
void C_standardise
(
    const int PQ,
    double *linstat,            /* in place standardisation */
    const double *expect,
    const double *covar,
    const int varonly,
    const double tol
) {

    double var;

    for (int p = 0; p < PQ; p++) {
        if (varonly) {
            var = covar[p];
        } else {
            var = covar[S(p, p, PQ)];
        }
        if (var > tol) {
            linstat[p] = (linstat[p] - expect[p]) / sqrt(var);
        } else {
            linstat[p] = NAN;
        }
    }
}
@|C_standardise
@}

@d P-Values
@{
@<C\_chisq\_pvalue@>
@<C\_perm\_pvalue@>
@<C\_norm\_pvalue@>
@}

@d C\_chisq\_pvalue
@{
/* lower = 1 means p-value, lower = 0 means 1 - p-value */
double C_chisq_pvalue
(
    const double stat,
    const int df,
    const int lower,
    const int give_log
) {
    return(pchisq(stat, (double) df, lower, give_log));
}
@|C_chisq_pvalue
@}

@d C\_perm\_pvalue
@{
double C_perm_pvalue
(
    const int greater,
    const double nperm,
    const int lower,
    const int give_log
) {

    double ret;

    if (give_log) {
         if (lower) {
             ret = log1p(- (double) greater / nperm);
         } else {
             ret = log(greater) - log(nperm);
         }
    } else {
        if (lower) {
            ret = 1.0 - (double) greater / nperm;
        } else {
            ret = (double) greater / nperm;
        }
    }
    return(ret);
}
@|C_perm_pvalue
@}

@d C\_norm\_pvalue
@{
double C_norm_pvalue
(
    const double stat,
    const int alternative,
    const int lower,
    const int give_log
) {

    double ret;

    if (alternative == ALTERNATIVE_less) {
        return(pnorm(stat, 0.0, 1.0, 1 - lower, give_log));
    } else if (alternative == ALTERNATIVE_greater) {
        return(pnorm(stat, 0.0, 1.0, lower, give_log));
    } else if (alternative == ALTERNATIVE_twosided) {
        if (lower) {
            ret = pnorm(fabs(stat)*-1.0, 0.0, 1.0, 1, 0);
            if (give_log) {
                return(log1p(- 2 * ret));
            } else {
                return(1 - 2 * ret);
            }
        } else {
            ret = pnorm(fabs(stat)*-1.0, 0.0, 1.0, 1, give_log);
            if (give_log) {
                return(ret + log(2));
            } else {
                return(2 * ret);
            }
        }
    }
    return(NA_REAL);
}
@}

@d C\_maxtype\_pvalue
@{
double C_maxtype_pvalue
(
    const double stat,
    const double *Covariance,
    const int n,
    const int alternative,
    const int lower,
    const int give_log,
    int maxpts, /* const? */
    double releps,
    double abseps,
    double tol
) {

    int nu = 0, inform, i, j, sub, nonzero, *infin, *index, rnd = 0;
    double ans, myerror, *lowerbnd, *upperbnd, *delta, *corr, *sd;

    /* univariate problem */
    if (n == 1)
        return(C_norm_pvalue(stat, alternative, lower, give_log));

    if (n == 2)
         corr = Calloc(1, double);
    else
         corr = Calloc(n + ((n - 2) * (n - 1))/2, double);

    sd = Calloc(n, double);
    lowerbnd = Calloc(n, double);
    upperbnd = Calloc(n, double);
    infin = Calloc(n, int);
    delta = Calloc(n, double);
    index = Calloc(n, int);

    /* determine elements with non-zero variance */

    nonzero = 0;
    for (i = 0; i < n; i++) {
        if (Covariance[S(i, i, n)] > tol) {
            index[nonzero] = i;
            nonzero++;
        }
    }

    /* mvtdst assumes the unique elements of the triangular
       covariance matrix to be passes as argument CORREL
    */

    for (int nz = 0; nz < nonzero; nz++) {

        /* handle elements with non-zero variance only */
        i = index[nz];

        /* standard deviations */
        sd[i] = sqrt(Covariance[S(i, i, n)]);

        if (alternative == ALTERNATIVE_less) {
            lowerbnd[nz] = stat;
            upperbnd[nz] = R_PosInf;
            infin[nz] = 1;
        } else if (alternative == ALTERNATIVE_greater) {
            lowerbnd[nz] = R_NegInf;
            upperbnd[nz] = stat;
            infin[nz] = 0;
        } else if (alternative == ALTERNATIVE_twosided) {
            lowerbnd[nz] = fabs(stat) * -1.0;
            upperbnd[nz] = fabs(stat);
            infin[nz] = 2;
        }

        delta[nz] = 0.0;

        /* set up vector of correlations, i.e., the upper
           triangular part of the covariance matrix) */
        for (int jz = 0; jz < nz; jz++) {
            j = index[jz];
            sub = (int) (jz + 1) + (double) ((nz - 1) * nz) / 2 - 1;
            if (sd[i] == 0.0 || sd[j] == 0.0)
                corr[sub] = 0.0;
            else
                corr[sub] = Covariance[S(i, j, n)] / (sd[i] * sd[j]);
        }
    }

    /* call mvtnorm's mvtdst C function defined in mvtnorm/include/mvtnormAPI.h */
    mvtnorm_C_mvtdst(&nonzero, &nu, lowerbnd, upperbnd, infin, corr, delta,
                     &maxpts, &abseps, &releps, &myerror, &ans, &inform, &rnd);

    /* inform == 0 means: everything is OK */
    switch (inform) {
        case 0: break;
        case 1: warning("cmvnorm: completion with ERROR > EPS"); break;
        case 2: warning("cmvnorm: N > 1000 or N < 1");
                ans = 0.0;
                break;
        case 3: warning("cmvnorm: correlation matrix not positive semi-definite");
                ans = 0.0;
                break;
        default: warning("cmvnorm: unknown problem in MVTDST");
                ans = 0.0;
    }
    Free(corr); Free(sd); Free(lowerbnd); Free(upperbnd);
    Free(infin); Free(delta);

    /* ans = 1 - p-value */
    if (lower) {
        if (give_log)
            return(log(ans)); /* log(1 - p-value) */
        return(ans); /* 1 - p-value */
    } else {
        if (give_log)
            return(log1p(ans)); /* log(p-value) */
        return(1 - ans); /* p-value */
    }
}
@|C_maxtype_pvalue
@}

@d maxstat Xfactor Variables
@{
SEXP LECV,
const int minbucket,
const int teststat,
int *wmax,
double *maxstat,
double *pval,
const int lower,
const int give_log
@}

@d C\_ordered\_Xfactor
@{
void C_ordered_Xfactor
(
@<maxstat Xfactor Variables@>
) {

    @<Setup maxstat Variables@> 

    @<Setup maxstat Memory@> 

    wmax[0] = NA_INTEGER;

    for (int p = 0; p < P; p++) {
        sumleft += ExpX[p];
        sumright -= ExpX[p];

        for (int q = 0; q < Q; q++) {
            mlinstat[q] += linstat[q * P + p];
            for (int np = 0; np < nperm; np++)
                mblinstat[q + np * Q] += blinstat[q * P + p + np * PQ];
            mexpect[q] += expect[q * P + p];
            if (Lb == 1) {
                @<Compute maxstat Variance / Covariance Directly@>
            } else {
                @<Compute maxstat Variance / Covariance from Total Covariance@>
            }
        }

        if ((sumleft >= minbucket) && (sumright >= minbucket) && (ExpX[p] > 0)) {

            ls = mlinstat; 
            @<Compute maxstat Test Statistic@>
            if (tmp > maxstat[0]) {
                wmax[0] = p;
                maxstat[0] = tmp;
            }

            for (int np = 0; np < nperm; np++) {
                ls = mblinstat + np * Q;
                @<Compute maxstat Test Statistic@>
                if (tmp > bmaxstat[np])
                    bmaxstat[np] = tmp;
            }
        }
    }
    @<Compute maxstat Permutation P-Value@>
    Free(mlinstat); Free(mexpect); Free(mblinstat); Free(bmaxstat);
    Free(mvar); Free(mcovar); Free(mMPinv);
}
@}

@d Setup maxstat Variables
@{
double *linstat, *expect, *covar, *varinf, *covinf, *ExpX, *blinstat, tol, *ls;
int P, Q, Lb;
R_xlen_t nperm;

double *mlinstat, *mblinstat, *mexpect, *mvar, *mcovar, *mMPinv, *bmaxstat,
       tmp, sumleft, sumright, sumweights;
int rank, PQ, greater;

Q = C_get_Q(LECV);
P = C_get_Q(LECV);
PQ = P * Q;
Lb = C_get_Lb(LECV);
if (Lb > 1 && C_get_varonly(LECV)) 
    error("need covarinance for maximally statistics with blocks");
linstat = C_get_LinearStatistic(LECV);
expect = C_get_Expectation(LECV);
covar = C_get_Covariance(LECV);
ExpX = C_get_ExpectationX(LECV);
varinf = C_get_VarianceInfluence(LECV);
covinf = C_get_CovarianceInfluence(LECV);
nperm = C_get_nperm(LECV);
if (nperm > 0)
    blinstat = C_get_PermutedLinearStatistic(LECV);
tol = C_get_tol(LECV);
@}

@d Setup maxstat Memory
@{
mlinstat = Calloc(Q, double);
mexpect = Calloc(Q, double);
if (teststat == TESTSTAT_maximum) {
   mvar = Calloc(Q, double);
   /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
   mcovar = Calloc(1, double);
   mMPinv = Calloc(1, double);
} else {
   mcovar = Calloc(Q * (Q + 1) / 2, double);
   mMPinv = Calloc(Q * (Q + 1) / 2, double);
   /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
   mvar = Calloc(1, double);
}
if (nperm > 0) {
    mblinstat = Calloc(Q * nperm, double);
    bmaxstat = Calloc(nperm, double);
} else { /* not needed, but allocate anyway to make -Wmaybe-uninitialized happy */
    mblinstat = Calloc(1, double);
    bmaxstat = Calloc(1, double);
}

maxstat[0] = 0.0;

for (int q = 0; q < Q; q++) {
    mlinstat[q] = 0.0;
    mexpect[q] = 0.0;
    if (teststat == TESTSTAT_maximum)
        mvar[q] = 0.0;
    for (int np = 0; np < nperm; np++) mblinstat[q + np * Q] = 0.0;
}
if (teststat == TESTSTAT_quadratic) {
    for (int q = 0; q < Q * (Q + 1) / 2; q++)
        mcovar[q] = 0.0;
}

sumleft = 0.0;
sumright = 0.0;
for (int p = 0; p < P; p++)
    sumright += ExpX[p];
sumweights = sumright;
@}

@d Compute maxstat Variance / Covariance from Total Covariance
@{
if (teststat == TESTSTAT_maximum) {
    for (int pp = 0; pp < p; pp++)
        mvar[q] += 2 * covar[S(pp + q * P, p + P * q, P * Q)];
     mvar[q] += covar[S(p + q * P, p + P * q, P * Q)];
} else {
     for (int qq = 0; qq <= q; qq++) {
         for (int pp = 0; pp < p; pp++)
             mcovar[S(q, qq, Q)] += 2 * covar[S(pp + q * P, p + P * qq, P * Q)];
         mcovar[S(q, qq, Q)] += covar[S(p + q * P, p + P * qq, P * Q)];
     }
}
@}

@d Compute maxstat Variance / Covariance Directly
@{
/* does not work with blocks! */
if (teststat == TESTSTAT_maximum) {
    C_VarianceLinearStatistic(1, Q, varinf, &sumleft, &sumleft,
                              sumweights, 0, mvar);
} else {
    C_CovarianceLinearStatistic(1, Q, covinf, &sumleft, &sumleft,
                                sumweights, 0, mcovar);
}
@}

@d Compute maxstat Test Statistic
@{
if (teststat == TESTSTAT_maximum) {
    tmp = C_maxtype(Q, ls, mexpect, mvar, 1, tol,
                    ALTERNATIVE_twosided);
} else {
    C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
    tmp = C_quadform(Q, ls, mexpect, mMPinv);
}
@}

@d Compute maxstat Permutation P-Value
@{
if (nperm > 0) {
    greater = 0;
    for (int np = 0; np < nperm; np++) {
        if (bmaxstat[np] > maxstat[0]) greater++;
    }
    pval[0] = C_perm_pvalue(greater, nperm, lower, give_log);
}
@}

@d C\_unordered\_Xfactor
@{
void C_unordered_Xfactor
(
@<maxstat Xfactor Variables@>
) {

    double *mtmp;
    int qPp, nc, *levels, Pnonzero, *indl, *contrast;

    @<Setup maxstat Variables@>

    @<Setup maxstat Memory@>
    mtmp = Calloc(P, double);

    for (int p = 0; p < P; p++) wmax[p] = NA_INTEGER;

    @<Count Levels@>

    for (int j = 1; j < mi; j++) { /* go though all splits */
 
        @<Setup unordered maxstat Contrasts@>

        @<Compute unordered maxstat Linear Statistic and Expectation@>

        if (Lb == 1) {
            @<Compute unordered maxstat Variance / Covariance from Total Covariance@>
        } else {
            @<Compute unordered maxstat Variance / Covariance Directly@>
        }

        if ((sumleft >= minbucket) && (sumright >= minbucket)) {

            ls = mlinstat;
            @<Compute maxstat Test Statistic@>
            if (tmp > maxstat[0]) {
                for (int p = 0; p < Pnonzero; p++)
                    wmax[levels[p]] = contrast[levels[p]];
                maxstat[0] = tmp;
            }

            for (int np = 0; np < nperm; np++) {
                ls = mblinstat + np * Q;
                @<Compute maxstat Test Statistic@>
                if (tmp > bmaxstat[np])
                    bmaxstat[np] = tmp;
            }
        }
    }

    @<Compute maxstat Permutation P-Value@>

    Free(mlinstat); Free(mexpect); Free(levels); Free(contrast); Free(indl); Free(mtmp);
    Free(mblinstat); Free(bmaxstat); Free(mvar); Free(mcovar); Free(mMPinv);
}


@}

@d Count Levels
@{
contrast = Calloc(P, int);
Pnonzero = 0;
for (int p = 0; p < P; p++) {
    if (ExpX[p] > 0) Pnonzero++;
}
levels = Calloc(Pnonzero, int);
nc = 0;
for (int p = 0; p < P; p++) {
    if (ExpX[p] > 0) {
        levels[nc] = p;
        nc++;
    }
}

if (Pnonzero >= 31)
    error("cannot search for unordered splits in >= 31 levels");

int mi = 1;
for (int l = 1; l < Pnonzero; l++) mi *= 2;
indl = Calloc(Pnonzero, int);
for (int p = 0; p < Pnonzero; p++) indl[p] = 0;
@}

@d Setup unordered maxstat Contrasts
@{
/* indl determines if level p is left or right */
int jj = j;
for (int l = 1; l < Pnonzero; l++) {
    indl[l] = (jj%2);
    jj /= 2;
}

sumleft = 0.0;
sumright = 0.0;
for (int p = 0; p < P; p++) contrast[p] = 0;
for (int p = 0; p < Pnonzero; p++) {
    sumleft += indl[p] * ExpX[levels[p]];
    sumright += (1 - indl[p]) * ExpX[levels[p]];
    contrast[levels[p]] = indl[p];
}
@}

@d Compute unordered maxstat Linear Statistic and Expectation
@{
for (int q = 0; q < Q; q++) {
    mlinstat[q] = 0.0;
    mexpect[q] = 0.0;
    for (int np = 0; np < nperm; np++)
        mblinstat[q + np * Q] = 0.0;
    for (int p = 0; p < P; p++) {
        qPp = q * P + p;
        mlinstat[q] += contrast[p] * linstat[qPp];
        mexpect[q] += contrast[p] * expect[qPp];
        for (int np = 0; np < nperm; np++)
            mblinstat[q + np * Q] += contrast[p] * blinstat[q * P + p + np * PQ];
    }
}
@}

@d Compute unordered maxstat Variance / Covariance from Total Covariance
@{
if (teststat == TESTSTAT_maximum) {
    for (int q = 0; q < Q; q++) {
        mvar[q] = 0.0;
        for (int p = 0; p < P; p++) {
            qPp = q * P + p;
            mtmp[p] = 0.0;
            for (int pp = 0; pp < P; pp++)
                mtmp[p] += contrast[pp] * covar[S(pp + q * P, qPp, PQ)];
	}
        for (int p = 0; p < P; p++)
            mvar[q] += contrast[p] * mtmp[p];
    }
} else {
    for (int q = 0; q < Q; q++) {
        for (int qq = 0; qq <= q; qq++)
            mcovar[S(q, qq, Q)] = 0.0;
        for (int qq = 0; qq <= q; qq++) {
            for (int p = 0; p < P; p++) {
                mtmp[p] = 0.0;
                for (int pp = 0; pp < P; pp++)
                    mtmp[p] += contrast[pp] * covar[S(pp + q * P, p + P * qq, P * Q)];
            }
            for (int p = 0; p < P; p++)
                mcovar[S(q, qq, Q)] += contrast[p] * mtmp[p];
        }
    }
}
@}

@d Compute unordered maxstat Variance / Covariance Directly
@{
if (teststat == TESTSTAT_maximum) {
    C_VarianceLinearStatistic(1, Q, varinf, &sumleft, &sumleft,
                              sumweights, 0, mvar);
} else {
    C_CovarianceLinearStatistic(1, Q, covinf, &sumleft, &sumleft,
                                sumweights, 0, mcovar);
}
@}


\section{Linear Statistics}

<<LinearStatistics>>=
a0 <- colSums(x[subset,r1] * y[subset,r2] * weights[subset])
a1 <- .Call("R_LinearStatistic", x, P, y, weights, subset, integer(0))
a2 <- .Call("R_LinearStatistic", x, P, y, as.double(weights), as.double(subset), integer(0))
a3 <- .Call("R_LinearStatistic", x, P, y, weights, as.double(subset), integer(0))
a4 <- .Call("R_LinearStatistic", x, P, y, as.double(weights), subset, integer(0))

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECV$LinearStatistic))


a0 <- as.vector(colSums(Xfactor[subset,r1Xfactor] * y[subset,r2Xfactor] * weights[subset]))
a1 <- .Call("R_LinearStatistic", ix, Lx, y, weights, subset, integer(0))
a2 <- .Call("R_LinearStatistic", ix, Lx, y, as.double(weights), as.double(subset), integer(0))
a3 <- .Call("R_LinearStatistic", ix, Lx, y, weights, as.double(subset), integer(0))
a4 <- .Call("R_LinearStatistic", ix, Lx, y, as.double(weights), subset, integer(0))

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a0 <- colSums(x[subset,r1] * y[subsety, r2])
a1 <- .Call("R_LinearStatistic", x, P, y, integer(0), subset, subsety)
a2 <- .Call("R_LinearStatistic", x, P, y, integer(0), as.double(subset), as.double(subsety))
stopifnot(all.equal(a0, a1) && all.equal(a0, a1))

a0 <- as.vector(colSums(Xfactor[subset,r1Xfactor] * y[subsety, r2Xfactor]))
a1 <- .Call("R_LinearStatistic", ix, Lx, y, integer(0), subset, subsety)
a1 <- .Call("R_LinearStatistic", ix, Lx, y, integer(0), as.double(subset), as.double(subsety))
stopifnot(all.equal(a0, a1))

@@

@d LinearStatistics
@{
@<RC\_LinearStatistic@>
@<R\_LinearStatistic@>
@}

@d R\_LinearStatistic
@{
SEXP R_LinearStatistic
(
    @<R x Input@>
    SEXP P,
    @<R y Input@>
    @<R weights Input@>,
    @<R subset Input@>,
    SEXP subsety
) 
{
    SEXP ans;
    @<C integer Q Input@>;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * Q));
    RC_LinearStatistic(x, N, INTEGER(P)[0], REAL(y), Q, 
                       weights, subset, Offset0, Nsubset, subsety, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@|R_LinearStatistic
@}

@d RC\_LinearStatistic Prototype
@{
void RC_LinearStatistic
(
    @<R x Input@>
    @<C integer N Input@>,
    @<C integer P Input@>,
    @<C real y Input@>
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
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
                    subset, offset, Nsubset, PQ_ans);
    } else {
        if (XLENGTH(weights) > 0) 
            error("weights given for permutation");
        if (XLENGTH(subset) != XLENGTH(subsety))
            error("incorrect subsets");
        RC_KronSums_Permutation(x, N, P, y, Q, subset, offset, Nsubset, 
                                subsety, PQ_ans);
    }
}
@|RC_LinearStatistic
@}

\section{Expectation and Covariance}

@d ExpectationCovariances
@{
@<RC\_ExpectationInfluence@>
@<R\_ExpectationInfluence@>
@<RC\_CovarianceInfluence@>
@<R\_CovarianceInfluence@>
@<RC\_ExpectationX@>
@<R\_ExpectationX@>
@<RC\_CovarianceX@>
@<R\_CovarianceX@>
@<C\_ExpectationLinearStatistic@>
@<C\_CovarianceLinearStatistic@>
@<C\_VarianceLinearStatistic@>
@}

\subsection{Linear Statistic}

@d C\_ExpectationLinearStatistic
@{
void C_ExpectationLinearStatistic
(
    @<C integer P Input@>,
    @<C integer Q Input@>,
    double *ExpInf,
    double *ExpX,
    const int add,
    double *PQ_ans
) {

    if (!add)
        for (int p = 0; p < P * Q; p++) PQ_ans[p] = 0.0;

    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++)
            PQ_ans[q * P + p] += ExpX[p] * ExpInf[q];
    }
}
@|C_ExpectationLinearStatistic
@}

@d C\_CovarianceLinearStatistic 
@{
void C_CovarianceLinearStatistic
(
    @<C integer P Input@>,
    @<C integer Q Input@>,
    double *CovInf,
    double *ExpX,
    double *CovX,
    @<C sumweights Input@>,
    const int add,
    double *PQPQ_sym_ans
) {

    double f1 = sumweights / (sumweights - 1);
    double f2 = 1.0 / (sumweights - 1);
    double tmp, *PP_sym_tmp;


    if (P * Q == 1) {
        tmp = f1 * CovInf[0] * CovX[0];
        tmp -= f2 * CovInf[0] * ExpX[0] * ExpX[0];
        if (add) {
            PQPQ_sym_ans[0] += tmp;
        } else {
            PQPQ_sym_ans[0] = tmp;
        }
    } else {
        PP_sym_tmp = Calloc(P * (P + 1) / 2, double);
        C_KronSums_sym_(ExpX, 1, P,
                        PP_sym_tmp);
        for (int p = 0; p < P * (P + 1) / 2; p++)
            PP_sym_tmp[p] = f1 * CovX[p] - f2 * PP_sym_tmp[p];
        C_kronecker_sym(CovInf, Q, PP_sym_tmp, P, 1 - (add >= 1),
                        PQPQ_sym_ans);
        Free(PP_sym_tmp);
    }
}
@|C_CovarianceLinearStatistic
@}

@d C\_VarianceLinearStatistic 
@{
void C_VarianceLinearStatistic
(
    @<C integer P Input@>,
    @<C integer Q Input@>,
    double *VarInf,
    double *ExpX,
    double *VarX,
    @<C sumweights Input@>,    
    const int add,
    double *PQ_ans
) {


    if (P * Q == 1) {
        C_CovarianceLinearStatistic(P, Q, VarInf, ExpX, VarX,
                                    sumweights, (add >= 1),
                                    PQ_ans);
    } else {
        double *P_tmp;
        P_tmp = Calloc(P, double);
        double f1 = sumweights / (sumweights - 1);
        double f2 = 1.0 / (sumweights - 1);
        for (int p = 0; p < P; p++)
            P_tmp[p] = f1 * VarX[p] - f2 * ExpX[p] * ExpX[p];
        C_kronecker(VarInf, 1, Q, P_tmp, 1, P, 1 - (add >= 1),
                    PQ_ans);
        Free(P_tmp);
    }
}
@|C_VarianceLinearStatistic
@}


\subsection{Influence}

<<ExpectationCovarianceInfluence>>=
sumweights <- sum(weights[subset])
expecty <- a0 <- colSums(y[subset, ] * weights[subset]) / sumweights
a1 <- .Call("R_ExpectationInfluence", y, weights, subset);
a2 <- .Call("R_ExpectationInfluence", y, as.double(weights), as.double(subset));
a3 <- .Call("R_ExpectationInfluence", y, weights, as.double(subset));
a4 <- .Call("R_ExpectationInfluence", y, as.double(weights), subset);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECV$ExpectationInfluence))
@@

@d R\_ExpectationInfluence
@{
SEXP R_ExpectationInfluence
(
    @<R y Input@>
    @<R weights Input@>,
    @<R subset Input@>
) 
{
    SEXP ans;
    @<C integer Q Input@>;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;
    double sumweights;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    sumweights = RC_Sums(N, weights, subset, Offset0, Nsubset);

    PROTECT(ans = allocVector(REALSXP, Q));
    RC_ExpectationInfluence(N, y, Q, weights, subset, Offset0, Nsubset, sumweights, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@|R_ExpectationInfluence
@}

@d RC\_ExpectationInfluence Prototype
@{
void RC_ExpectationInfluence
(
    @<C integer N Input@>,
    @<R y Input@>
    @<C integer Q Input@>,
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
    @<C sumweights Input@>,    
    @<C colSums Answer@>
) 
@}

@d RC\_ExpectationInfluence
@{
@<RC\_ExpectationInfluence Prototype@>
{
    double center;

    RC_colSums(REAL(y), N, Q, Power1, &center, !DoCenter, weights, 
               subset, offset, Nsubset, P_ans);
    for (int q = 0; q < Q; q++) 
        P_ans[q] = P_ans[q] / sumweights;
}
@|RC_ExpectationInfluence
@}

<<CovarianceInfluence>>=
sumweights <- sum(weights[subset])
yc <- t(t(y) - expecty)
r1y <- rep(1:ncol(y), ncol(y))
r2y <- rep(1:ncol(y), each = ncol(y))
a0 <- colSums(yc[subset, r1y] * yc[subset, r2y] * weights[subset]) / sumweights
a0 <- matrix(a0, ncol = ncol(y))
vary <- diag(a0)
a0 <- a0[lower.tri(a0, diag = TRUE)]
a1 <- .Call("R_CovarianceInfluence", y, weights, subset, 0L);
a2 <- .Call("R_CovarianceInfluence", y, as.double(weights), as.double(subset), 0L);
a3 <- .Call("R_CovarianceInfluence", y, weights, as.double(subset), 0L);
a4 <- .Call("R_CovarianceInfluence", y, as.double(weights), subset, 0L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECV$CovarianceInfluence) &&
          all.equal(vary, LEV$VarianceInfluence))

a1 <- .Call("R_CovarianceInfluence", y, weights, subset, 1L);
a2 <- .Call("R_CovarianceInfluence", y, as.double(weights), as.double(subset), 1L);
a3 <- .Call("R_CovarianceInfluence", y, weights, as.double(subset), 1L);
a4 <- .Call("R_CovarianceInfluence", y, as.double(weights), subset, 1L);

a0 <- vary

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LEV$VarianceInfluence))
@@

@d R\_CovarianceInfluence
@{
SEXP R_CovarianceInfluence
(
    @<R y Input@>
    @<R weights Input@>,
    @<R subset Input@>,
    SEXP varonly
) 
{
    SEXP ans;
    SEXP ExpInf;
    @<C integer Q Input@>;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;
    double sumweights;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    PROTECT(ExpInf = R_ExpectationInfluence(y, weights, subset));

    sumweights = RC_Sums(N, weights, subset, Offset0, Nsubset);

    if (INTEGER(varonly)[0]) {
        PROTECT(ans = allocVector(REALSXP, Q));
    } else {
        PROTECT(ans = allocVector(REALSXP, Q * (Q + 1) / 2));
    }
    RC_CovarianceInfluence(N, y, Q, weights, subset, Offset0, Nsubset, REAL(ExpInf), sumweights, 
                           INTEGER(varonly)[0], REAL(ans));
    UNPROTECT(2);
    return(ans);
}
@|R_CovarianceInfluence
@}

@d RC\_CovarianceInfluence Prototype
@{
void RC_CovarianceInfluence
(
    @<C integer N Input@>,
    @<R y Input@>
    @<C integer Q Input@>,
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
    double *ExpInf,
    @<C sumweights Input@>,
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

\subsection{X}

<<ExpectationCovarianceX>>=
a0 <- colSums(x[subset, ] * weights[subset]) 
a0
a1 <- .Call("R_ExpectationX", x, P, weights, subset);
a2 <- .Call("R_ExpectationX", x, P, as.double(weights), as.double(subset));
a3 <- .Call("R_ExpectationX", x, P, weights, as.double(subset));
a4 <- .Call("R_ExpectationX", x, P, as.double(weights), subset);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECV$ExpectationX))

a0 <- colSums(x[subset, ]^2 * weights[subset]) 
a1 <- .Call("R_CovarianceX", x, P, weights, subset, 1L);
a2 <- .Call("R_CovarianceX", x, P, as.double(weights), as.double(subset), 1L);
a3 <- .Call("R_CovarianceX", x, P, weights, as.double(subset), 1L);
a4 <- .Call("R_CovarianceX", x, P, as.double(weights), subset, 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a0 <- as.vector(colSums(Xfactor[subset, ] * weights[subset]))
a0
a1 <- .Call("R_ExpectationX", ix, Lx, weights, subset);
a2 <- .Call("R_ExpectationX", ix, Lx, as.double(weights), as.double(subset));
a3 <- .Call("R_ExpectationX", ix, Lx, weights, as.double(subset));
a4 <- .Call("R_ExpectationX", ix, Lx, as.double(weights), subset);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a1 <- .Call("R_CovarianceX", ix, Lx, weights, subset, 1L);
a2 <- .Call("R_CovarianceX", ix, Lx, as.double(weights), as.double(subset), 1L);
a3 <- .Call("R_CovarianceX", ix, Lx, weights, as.double(subset), 1L);
a4 <- .Call("R_CovarianceX", ix, Lx, as.double(weights), subset, 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

r1x <- rep(1:ncol(Xfactor), ncol(Xfactor))
r2x <- rep(1:ncol(Xfactor), each = ncol(Xfactor))
a0 <- colSums(Xfactor[subset, r1x] * Xfactor[subset, r2x] * weights[subset])
a0 <- matrix(a0, ncol = ncol(Xfactor))
vary <- diag(a0)
a0 <- a0[upper.tri(a0, diag = TRUE)]

a0
a1 <- .Call("R_CovarianceX", ix, Lx, weights, subset, 0L);
a1
a2 <- .Call("R_CovarianceX", ix, Lx, as.double(weights), as.double(subset), 0L);
a3 <- .Call("R_CovarianceX", ix, Lx, weights, as.double(subset), 0L);
a4 <- .Call("R_CovarianceX", ix, Lx, as.double(weights), subset, 0L);

#stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
#          all.equal(a0, a3) && all.equal(a0, a4))


@@

@d R\_ExpectationX
@{
SEXP R_ExpectationX
(
    @<R x Input@>
    SEXP P,
    @<R weights Input@>,
    @<R subset Input@>
) 
{
    SEXP ans;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;

    N = XLENGTH(x) / INTEGER(P)[0];
    Nsubset = XLENGTH(subset);

    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0]));
    RC_ExpectationX(x, N, INTEGER(P)[0], weights, subset, 
                    Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@|R_ExpectationX
@}

@d RC\_ExpectationX Prototype
@{
void RC_ExpectationX
(
    @<R x Input@>
    @<C integer N Input@>,
    @<C integer P Input@>,
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
    @<C OneTableSums Answer@>
) 
@}

@d RC\_ExpectationX
@{
@<RC\_ExpectationX Prototype@>
{
    double center;

    if (TYPEOF(x) == INTSXP) {
        double* Pp1tmp = Calloc(P + 1, double);
        RC_OneTableSums(INTEGER(x), N, P + 1, weights, subset, offset, Nsubset, Pp1tmp);
        for (int p = 0; p < P; p++) P_ans[p] = Pp1tmp[p + 1];
        Free(Pp1tmp);
    } else {
        RC_colSums(REAL(x), N, P, Power1, &center, !DoCenter, weights, subset, offset, Nsubset, P_ans);
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
    @<R weights Input@>,
    @<R subset Input@>,
    SEXP varonly
) 
{
    SEXP ans;
    SEXP ExpX;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;

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
@|R_CovarianceX
@}

@d RC\_CovarianceX Prototype
@{
void RC_CovarianceX
(
    @<R x Input@>
    @<C integer N Input@>,
    @<C integer P Input@>,
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
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

\section{Computing Sums}

The core concept of all functions in the section is the computation of
various sums over observations, weights, or blocks. We start with an
initialisation of the loop over all observations

@d init subset loop
@{
    R_xlen_t diff = 0;
    s = subset + offset;
    w = weights;
    /* subset is R-style index in 1:N */
    if (Nsubset > 0)
        diff = (R_xlen_t) s[0] - 1;
@}

and loop over $i = 1, \dots, N$ when no subset was specified or
over the subset of the subset given by \verb|offset| and \verb|Nsubset|,
allowing for number of observations larger than \verb|INT_MAX|

@d start subset loop
@{
    for (R_xlen_t i = 0; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) 
@}

After computions in the loop, we compute the next element

@d continue subset loop
@{
    if (Nsubset > 0) {
        /* NB: diff also works with R style index */
        diff = (R_xlen_t) s[1] - s[0];
        if (diff < 0)
            error("subset not sorted");
        s++;
    } else {
        diff = 1;
    }
@}


\subsection{Simple Sums}

<<SimpleSums>>=
a0 <- sum(weights[subset])
a1 <- .Call("R_Sums", N, weights, subset)
a2 <- .Call("R_Sums", N, as.double(weights), as.double(subset))
a3 <- .Call("R_Sums", N, weights, as.double(subset))
a4 <- .Call("R_Sums", N, as.double(weights), subset)
stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@


@d SimpleSums
@{
@<C\_Sums\_dweights\_dsubset@>
@<C\_Sums\_iweights\_dsubset@>
@<C\_Sums\_iweights\_isubset@>
@<C\_Sums\_dweights\_isubset@>
@<RC\_Sums@>
@<R\_Sums@>
@}


@d R\_Sums
@{
SEXP R_Sums
(
    @<R N Input@>
    @<R weights Input@>,
    @<R subset Input@>
) 
{
    SEXP ans;
    @<C integer Nsubset Input@>;

    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = RC_Sums(INTEGER(N)[0], weights, subset, Offset0, Nsubset);
    UNPROTECT(1);

    return(ans);
}
@|R_Sums
@}

@d RC\_Sums Prototype
@{
double RC_Sums
(
    @<C integer N Input@>,
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>
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
            return((double) Nsubset);
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
    @<C integer N Input@>,
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
    @<C integer N Input@>,
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
    @<C integer N Input@>,
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
    @<C integer N Input@>,
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
a1 <- .Call("R_KronSums", x, P, y, weights, subset)
a2 <- .Call("R_KronSums", x, P, y, as.double(weights), as.double(subset))
a3 <- .Call("R_KronSums", x, P, y, weights, as.double(subset))
a4 <- .Call("R_KronSums", x, P, y, as.double(weights), subset)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@

@d KronSums
@{
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
@<C\_KronSums\_Permutation\_isubset@>
@<C\_KronSums\_Permutation\_dsubset@>
@<C\_XfactorKronSums\_Permutation\_isubset@>
@<C\_XfactorKronSums\_Permutation\_dsubset@>
@<RC\_KronSums\_Permutation@>
@<R\_KronSums\_Permutation@>
@}


@d R\_KronSums
@{
SEXP R_KronSums
(
    @<R x Input@>
    SEXP P,
    @<R y Input@>
    @<R weights Input@>,
    @<R subset Input@>
) 
{
    SEXP ans;
    @<C integer Q Input@>;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;

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
@|R_KronSums
@}

@d RC\_KronSums Prototype
@{
void RC_KronSums
(
    @<RC KronSums Input@>
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
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
    @<C integer N Input@>,
    @<C integer P Input@>,
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
    int idx;

    for (int p = 0; p < P; p++) {
        for (int q = (SYMMETRIC ? p : 0); q < Q; q++) {
            /* SYMMETRIC is column-wise, default
               is row-wise (maybe need to change this) */
            if (SYMMETRIC) {
                idx = S(p, q, P); 
            } else {
                idx = q * P + p;  
            }
            PQ_ans[idx] = 0.0;
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
                    PQ_ans[idx] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                } else {
                    PQ_ans[idx] += (xx[0] - cx) * (yy[0] - cy);
                }
                @<continue subset loop@>
            }
            xx = xx + diff;
            yy = yy + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                PQ_ans[idx] += (xx[0] - cx) * (yy[0] - cy) * w[0];
            } else {
                PQ_ans[idx] += (xx[0] - cx) * (yy[0] - cy);
            }
        }
    }
@}

\subsubsection{Xfactor Kronecker Sums}

<<XfactorKronSums>>=

a0 <- as.vector(colSums(Xfactor[subset,r1Xfactor] * 
                        y[subset,r2Xfactor] * weights[subset]))
a1 <- .Call("R_KronSums", ix, Lx, y, weights, subset)
a2 <- .Call("R_KronSums", ix, Lx, y, as.double(weights), as.double(subset))
a3 <- .Call("R_KronSums", ix, Lx, y, weights, as.double(subset))
a4 <- .Call("R_KronSums", ix, Lx, y, as.double(weights), subset)

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


\subsubsection{Permuted Kronecker Sums}

<<KronSums-Permutation>>=
a0 <- colSums(x[subset,r1] * y[subsety, r2])
a1 <- .Call("R_KronSums_Permutation", x, P, y, subset, subsety)
a2 <- .Call("R_KronSums_Permutation", x, P, y, as.double(subset), as.double(subsety))
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
    @<C integer Q Input@>;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * Q));
    RC_KronSums_Permutation(x, N, INTEGER(P)[0], REAL(y), Q, subset, Offset0, Nsubset, 
                            subsety, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@|R_KronSums_Permutation
@}

@d RC\_KronSums\_Permutation Prototype
@{
void RC_KronSums_Permutation
(
    @<R x Input@>
    @<C integer N Input@>,
    @<C integer P Input@>,
    @<C real y Input@>
    @<R subset Input@>,
    @<C subset range Input@>,
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
            /* subset is R-style index in 1:N */
            diff = (R_xlen_t) sx[offset] - 1;
            for (R_xlen_t i = offset; i < Nsubset - 1; i++) {
                xx = xx + diff;
                /* subsety is R-style index in 1:N; - 1 needed here */
                PQ_ans[0] += xx[0] * y[ (R_xlen_t) sy[0] - 1 + q * N];
                diff = (R_xlen_t) sx[1] - sx[0];
                sx++;
                sy++;
            }
            xx = xx + diff;
            PQ_ans[0] += xx[0] * y[ (R_xlen_t) sy[0] - 1 + q * N];
            PQ_ans++;
        }
    }
@}

\subsubsection{Xfactor Permuted Kronecker Sums}

<<XfactorKronSums-Permutation>>=
a0 <- as.vector(colSums(Xfactor[subset,r1Xfactor] * y[subsety, r2Xfactor]))
a1 <- .Call("R_KronSums_Permutation", ix, Lx, y, subset, subsety)
a1 <- .Call("R_KronSums_Permutation", ix, Lx, y, as.double(subset), as.double(subsety))
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
        /* subset is R-style index in 1:N */
        diff = (R_xlen_t) sx[offset] - 1;
        for (R_xlen_t i = offset; i < Nsubset - 1; i++) {
            xx = xx + diff;
            /* subsety is R-style index in 1:N; - 1 needed here */
            PQ_ans[(xx[0] - 1) + q * P] += y[ (R_xlen_t) sy[0] - 1 + q * N];
            diff = (R_xlen_t) sx[1] - sx[0];
            sx++;
            sy++;
        }
        xx = xx + diff;
        PQ_ans[(xx[0] - 1) + q * P] += y[ (R_xlen_t) sy[0] - 1 + q * N];
    }
@}



\subsection{Column Sums}

<<colSums>>=
a0 <- colSums(x[subset,] * weights[subset])
a1 <- .Call("R_colSums", x, weights, subset)
a2 <- .Call("R_colSums", x, as.double(weights), as.double(subset))
a3 <- .Call("R_colSums", x, weights, as.double(subset))
a4 <- .Call("R_colSums", x, as.double(weights), subset)

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@


@d colSums
@{
@<C\_colSums\_dweights\_dsubset@>
@<C\_colSums\_iweights\_dsubset@>
@<C\_colSums\_iweights\_isubset@>
@<C\_colSums\_dweights\_isubset@>
@<RC\_colSums@>
@<R\_colSums@>
@}


@d R\_colSums
@{
SEXP R_colSums
(
    @<R x Input@>
    @<R weights Input@>,
    @<R subset Input@>
) {

    SEXP ans;
    int P;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;
    double center;

    P = NCOL(x);
    N = XLENGTH(x) / P;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, P));
    RC_colSums(REAL(x), N, P, Power1, &center, !DoCenter, weights, subset, Offset0, 
               Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@|R_colSums
@}

@d RC\_colSums Prototype
@{
void RC_colSums
(
    @<C colSums Input@>
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
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

\subsection{Tables}

\subsubsection{OneTable Sums}

<<OneTableSum>>=

a0 <- as.vector(xtabs(weights ~ ixf, subset = subset))
a1 <- .Call("R_OneTableSums", ix, Lx + 1L, weights, subset)[-1]
a2 <- .Call("R_OneTableSums", ix, Lx + 1L, 
            as.double(weights), as.double(subset))[-1]
a3 <- .Call("R_OneTableSums", ix, Lx + 1L, 
            weights, as.double(subset))[-1]
a4 <- .Call("R_OneTableSums", ix, Lx + 1L, 
            as.double(weights), subset)[-1]

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))
@@


@d Tables
@{
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
@}


@d R\_OneTableSums
@{
SEXP R_OneTableSums
(
    @<R x Input@>
    SEXP P,
    @<R weights Input@>,
    @<R subset Input@>
) {

    SEXP ans;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0]));
    RC_OneTableSums(INTEGER(x), N, INTEGER(P)[0],  
                 weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@|R_OneTableSums
@}


@d RC\_OneTableSums Prototype
@{
void RC_OneTableSums
(
    @<C OneTableSums Input@>
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
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

\subsubsection{TwoTable Sums}

<<TwoTableSum>>=

a0 <- c(xtabs(weights ~ ixf + iyf, subset = subset))
a1 <- .Call("R_TwoTableSums", ix, Lx + 1L, iy, Ly + 1L, 
            weights, subset)
a1 <- c(matrix(a1, nrow = Lx + 1, ncol = Ly + 1)[-1,-1])
a2 <- .Call("R_TwoTableSums", ix, Lx + 1L, iy, Ly + 1L, 
            as.double(weights), as.double(subset))
a2 <- c(matrix(a2, nrow = Lx + 1, ncol = Ly + 1)[-1,-1])
a3 <- .Call("R_TwoTableSums", ix, Lx + 1L, iy, Ly + 1L, 
            weights, as.double(subset))
a3 <- c(matrix(a3, nrow = Lx + 1, ncol = Ly + 1)[-1,-1])
a4 <- .Call("R_TwoTableSums", ix, Lx + 1L, iy, Ly + 1L, 
            as.double(weights), subset)
a4 <- c(matrix(a4, nrow = Lx + 1, ncol = Ly + 1)[-1,-1])

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
    @<R weights Input@>,
    @<R subset Input@>
) {

    SEXP ans;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * INTEGER(Q)[0]));
    RC_TwoTableSums(INTEGER(x), N, INTEGER(P)[0], INTEGER(y), INTEGER(Q)[0], 
                 weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@|R_TwoTableSums
@}


@d RC\_TwoTableSums Prototype
@{
void RC_TwoTableSums
(
    @<C TwoTableSums Input@>
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
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

\subsubsection{ThreeTable Sums}

<<ThreeTableSum>>=

a0 <- c(xtabs(weights ~ ixf + iyf + block, subset = subset))
a1 <- .Call("R_ThreeTableSums", ix, Lx + 1L, iy, Ly + 1L, 
            block, B, weights, subset)
a1 <- c(array(a1, dim = c(Lx + 1, Ly + 1, B))[-1,-1,])
a2 <- .Call("R_ThreeTableSums", ix, Lx + 1L, iy, Ly + 1L, 
            block, B,
            as.double(weights), as.double(subset))
a2 <- c(array(a2, dim = c(Lx + 1, Ly + 1, B))[-1,-1,])
a3 <- .Call("R_ThreeTableSums", ix, Lx + 1L, iy, Ly + 1L, 
            block, B,
            weights, as.double(subset))
a3 <- c(array(a3, dim = c(Lx + 1, Ly + 1, B))[-1,-1,])
a4 <- .Call("R_ThreeTableSums", ix, Lx + 1L, iy, Ly + 1L, 
            block, B,
            as.double(weights), subset)
a4 <- c(array(a4, dim = c(Lx + 1, Ly + 1, B))[-1,-1,])

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
    @<R weights Input@>,
    @<R subset Input@>
) {

    SEXP ans;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * INTEGER(Q)[0] * INTEGER(L)[0]));
    RC_ThreeTableSums(INTEGER(x), N, INTEGER(P)[0], INTEGER(y), INTEGER(Q)[0], 
                      INTEGER(block), INTEGER(L)[0],
                      weights, subset, Offset0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@|R_ThreeTableSums
@}


@d RC\_ThreeTableSums Prototype
@{
void RC_ThreeTableSums
(
    @<C ThreeTableSums Input@>
    @<R weights Input@>,
    @<R subset Input@>,
    @<C subset range Input@>,
    @<C ThreeTableSums Answer@>
) 
@}

@d RC\_ThreeTableSums
@{
@<RC\_ThreeTableSums Prototype@>
{
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_ThreeTableSums_iweights_isubset(x, N, P, y, Q, block, Lb, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQL_ans);
        } else {
            C_ThreeTableSums_iweights_dsubset(x, N, P, y, Q, block, Lb,
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQL_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_ThreeTableSums_dweights_isubset(x, N, P, y, Q, block, Lb,
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQL_ans);
        } else {
            C_ThreeTableSums_dweights_dsubset(x, N, P, y, Q, block, Lb,
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

    for (int p = 0; p < PQ * Lb; p++) PQL_ans[p] = 0.0;

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

\section{Utilities}



\subsection{Blocks}

<<Helpers>>=
a0 <- as.vector(do.call("c", tapply(subset, block[subset], function(s) s)))
a1 <- .Call("R_order_subset_wrt_block", y, subset, integer(0), block, B + 1L);
stopifnot(all.equal(a0, a1))
a2 <- .Call("R_order_subset_wrt_block", y, subset, integer(0), integer(0), 2L);
stopifnot(all.equal(subset, a2))
a0 <- as.vector(do.call("c", tapply(1:N, block, function(s) s)))
a3 <- .Call("R_order_subset_wrt_block", y, integer(0), integer(0), block, B + 1L);
stopifnot(all.equal(a0, a3))
a4 <- .Call("R_order_subset_wrt_block", y, integer(0), integer(0), integer(0), 2L);
stopifnot(all.equal(1:N, a4))
@@

@d Utils
@{
@<C\_setup\_subset@>
@<C\_setup\_subset\_block@>
@<C\_order\_subset\_wrt\_block@>
@<RC\_order\_subset\_wrt\_block@>
@<R\_order\_subset\_wrt\_block@>
@}

@d R\_order\_subset\_wrt\_block
@{
SEXP R_order_subset_wrt_block
(
    @<R y Input@>
    @<R subset Input@>,
    @<R weights Input@>,
    @<R block Input@>
    SEXP Nlevels
)
{

    @<C integer N Input@>;
    SEXP blockTable, ans;

    N = XLENGTH(y) / NCOL(y);

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
@|R_order_subset_wrt_block
@}


@d RC\_order\_subset\_wrt\_block Prototype
@{
SEXP RC_order_subset_wrt_block
(
    @<C integer N Input@>,
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
@|RC_order_subset_wrt_block
@}

@d C\_setup\_subset
@{
void C_setup_subset
(
    @<C integer N Input@>,
    SEXP ans
)
{
    for (R_xlen_t i = 0; i < N; i++) {
        /* ans is R style index in 1:N */
        if (TYPEOF(ans) == INTSXP) {
            INTEGER(ans)[i] = i + 1;
        } else {
            REAL(ans)[i] = (double) i + 1;
        }
    }
}
@|C_setup_subset
@}

@d C\_setup\_subset\_block
@{
void C_setup_subset_block
(
    @<C integer N Input@>,
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
        /* ans is R style index in 1:N */
        if (TYPEOF(ans) == INTSXP) {
            INTEGER(ans)[(int) cumtable[INTEGER(block)[i]]++] = i + 1;
        } else {
            REAL(ans)[(R_xlen_t) cumtable[INTEGER(block)[i]]++] = (double) i + 1;
        }
    }

    Free(cumtable);
}
@|C_setup_subset_block
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

    /* subset is R style index in 1:N */
    if (TYPEOF(subset) == INTSXP) {
        for (R_xlen_t i = 0; i < XLENGTH(subset); i++)
            INTEGER(ans)[(int) cumtable[INTEGER(block)[INTEGER(subset)[i] - 1]]++] = INTEGER(subset)[i];
    } else {
        for (R_xlen_t i = 0; i < XLENGTH(subset); i++)
            REAL(ans)[(R_xlen_t) cumtable[INTEGER(block)[(R_xlen_t) REAL(subset)[i] - 1]]++] = REAL(subset)[i];
    }

    Free(cumtable); 
}
@|C_order_subset_wrt_block
@}


@d R\_setup\_subset
@{
SEXP R_setup_subset
(
    @<R N Input@>
    @<R weights Input@>,
    @<R subset Input@>
)
{
    SEXP ans;

    PROTECT(ans = RC_setup_subset(INTEGER(N)[0], weights, subset));

    UNPROTECT(1);
    return(ans);
}
@|R_setup_subset
@}


@d RC\_setup\_subset Prototype
@{
SEXP RC_setup_subset
(
    @<C integer N Input@>,
    @<R weights Input@>,
    @<R subset Input@>
)
@}

Because this will only be used when really needed (in Permutations) we can
be a little bit more generous with memory here.

@d RC\_setup\_subset
@{
@<RC\_setup\_subset Prototype@>
{
    SEXP ans, mysubset;
    R_xlen_t sumweights;

    if (XLENGTH(weights) == 0 && XLENGTH(subset) > 0)
        return(subset);

    if (XLENGTH(subset) == 0) {
        PROTECT(mysubset = allocVector(REALSXP, N));
        C_setup_subset(N, mysubset);
    } else {
        PROTECT(mysubset = coerceVector(subset, REALSXP));
    }

    if (XLENGTH(weights) == 0) {
        UNPROTECT(1);
        return(mysubset);
    }
        
    sumweights = (R_xlen_t) RC_Sums(N, weights, mysubset, Offset0, XLENGTH(subset));
    PROTECT(ans = allocVector(REALSXP, sumweights));

    R_xlen_t itmp = 0;
    for (R_xlen_t i = 0; i < XLENGTH(mysubset); i++) {
        if (TYPEOF(weights) == REALSXP) {
            for (R_xlen_t j = 0; j < REAL(weights)[(R_xlen_t) REAL(mysubset)[i] - 1]; j++)
                REAL(ans)[itmp++] = REAL(mysubset)[i];
        } else {
            for (R_xlen_t j = 0; j < INTEGER(weights)[(R_xlen_t) REAL(mysubset)[i] - 1]; j++)
                REAL(ans)[itmp++] = REAL(mysubset)[i];
        }
    }
    UNPROTECT(2);
    return(ans);
}
@|RC_setup_subset
@}

\subsection{Permutation Helpers}

<<Permutations>>=
stopifnot(all(all.equal(1:N, .Call("R_setup_subset", N, integer(0), integer(0))),
              all.equal(rep(1:N, weights), .Call("R_setup_subset", N, weights,
integer(0))),
              all.equal(subset, .Call("R_setup_subset", N, integer(0),
subset)),
              all.equal(rep(subset, weights[subset]), .Call("R_setup_subset", N, weights,
              subset))))
@@

@d Permutations
@{
@<RC\_setup\_subset@>
@<R\_setup\_subset@>
@<C\_Permute@>
@<C\_doPermute@>
@<C\_PermuteBlock@>
@<C\_doPermuteBlock@>
@}

@d C\_Permute
@{
void C_Permute
(
    double *subset,
    @<C integer Nsubset Input@>,
    double *ans
) {

    R_xlen_t n = Nsubset, j;

    for (R_xlen_t i = 0; i < Nsubset; i++) {
        j = n * unif_rand();
        ans[i] = subset[j];
        subset[j] = subset[--n];
    }
}
@|C_Permute
@}

@d C\_doPermute
@{
void C_doPermute
(
    double *subset,
    @<C integer Nsubset Input@>,
    double *Nsubset_tmp,
    double *perm
) {
    Memcpy(Nsubset_tmp, subset, Nsubset);
    C_Permute(Nsubset_tmp, Nsubset, perm);
}
@|C_doPermute
@}

@d C\_PermuteBlock
@{
void C_PermuteBlock
(
    double *subset,
    double *table,
    int Nlevels,
    double *ans
) {

    double *px, *pans;

    px = subset;
    pans = ans;

    for (R_xlen_t j = 0; j < Nlevels; j++) {
        if (table[j] > 0) {
            C_Permute(px, (R_xlen_t) table[j], pans);
            px += (R_xlen_t) table[j];
            pans += (R_xlen_t) table[j];
        }
    }
}
@|C_PermuteBlock
@}

@d C\_doPermuteBlock
@{
void C_doPermuteBlock
(
    double *subset,
    @<C integer Nsubset Input@>,
    double *table,
    int Nlevels,
    double *Nsubset_tmp,
    double *perm
) {
    Memcpy(Nsubset_tmp, subset, Nsubset);
    C_PermuteBlock(Nsubset_tmp, table, Nlevels, perm);
}
@|C_doPermuteBlock
@}


\subsection{Other Utils}

@d MoreUtils
@{
@<NROW@>
@<NCOL@>
@<NLEVELS@>
@<C\_kronecker@>
@<C\_kronecker\_sym@>
@<C\_KronSums\_sym@>
@<C\_MPinv\_sym@>
@}


@d NROW
@{
int NROW
(
    SEXP x
) {

    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(XLENGTH(x));
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[0]);
    return(INTEGER(a)[0]);
}
@|NROW
@}

@d NCOL
@{
int NCOL
(
    SEXP x
) {

    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    if (TYPEOF(a) == REALSXP)
        return(REAL(a)[1]);
    return(INTEGER(a)[1]);
}
@|NCOL
@}

@d NLEVELS
@{
int NLEVELS	
(
    SEXP x
) {

   SEXP a;

    a = getAttrib(x, R_LevelsSymbol);
    if (a == R_NilValue)
        error("no levels attribute found");
   return(NROW(a));
}
@|NLEVELS
@}

@d C\_kronecker
@{
void C_kronecker
(
    const double *A,
    const int m,
    const int n,
    const double *B,
    const int r,
    const int s,
    const int overwrite,
    double *ans
) {

    int i, j, k, l, mr, js, ir;
    double y;

    if (overwrite) {
        for (i = 0; i < m * r * n * s; i++) ans[i] = 0.0;
    }

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++)
                    ans[(js + l) * mr + ir + k] += y * B[l * r + k];
            }
        }
    }
}
@|C_kronecker
@}

@d C\_kronecker\_sym
@{
void C_kronecker_sym
(
    const double *A,
    const int m,
    const double *B,
    const int r,
    const int overwrite,
    double *ans
) {

    int i, j, k, l, mr, js, ir, s;
    double y;

    mr = m * r;
    s = r;

    if (overwrite) {
        for (i = 0; i < mr * (mr + 1) / 2; i++) ans[i] = 0.0;
    }

    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j <= i; j++) {
            js = j * s;
            y = A[S(i, j, m)];
            for (k = 0; k < r; k++) {
                for (l = 0; l < (j < i ? s : k + 1); l++) {
                    ans[S(ir + k, js + l, mr)] += y * B[S(k, l, r)];
                }
            }
        }
    }
}
@|C_kronecker_sym
@}

@d C\_KronSums\_sym
@{
/* sum_i (t(x[i,]) %*% x[i,]) */
void C_KronSums_sym_
(
    @<C real x Input@>
    double *PP_sym_ans
) {

    int pN, qN, SpqP;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < N; i++)
                 PP_sym_ans[SpqP] +=  x[qN + i] * x[pN + i];
        }
    }
}
@|C_KronSums_sym
@}

@d C\_MPinv\_sym
@{
void C_MPinv_sym
(
    const double *x,
    const int n,
    const double tol,
    double *dMP,
    int *rank
) {

    double *val, *vec, dtol, *rx, *work, valinv;
    int valzero = 0, info = 0, kn;

    if (n == 1) {
        if (x[0] > tol) {
            dMP[0] = 1 / x[0];
            rank[0] = 1;
        } else {
            dMP[0] = 0;
            rank[0] = 0;
        }
    } else {
        rx = Calloc(n * (n + 1) / 2, double);
        Memcpy(rx, x, n * (n + 1) / 2);
        work = Calloc(3 * n, double);
        val = Calloc(n, double);
        vec = Calloc(n * n, double);

/*
        F77_CALL(dspev)("V", "L", &n, rx, val, vec, &n, work,
                        &info);
*/
        dtol = val[n - 1] * tol;

        for (int k = 0; k < n; k++)
            valzero += (val[k] < dtol);
        rank[0] = n - valzero;

        for (int k = 0; k < n * (n + 1) / 2; k++) dMP[k] = 0.0;

        for (int k = valzero; k < n; k++) {
            valinv = 1 / val[k];
            kn = k * n;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= i; j++) {
                    /* MP is symmetric */
                    dMP[S(i, j, n)] += valinv * vec[kn + i] * vec[kn + j];
                }
            }
        }
        Free(rx); Free(work); Free(val); Free(vec);
    }
}
@}

\section{Memory}

@d Memory
@{
@<C\_get\_P@>
@<C\_get\_Q@>
@<C\_get\_varonly@>
@<C\_get\_Xfactor@>
@<C\_get\_LinearStatistic@>
@<C\_get\_Expectation@>
@<C\_get\_Variance@>
@<C\_get\_Covariance@>
@<C\_get\_MPinv@>
@<C\_get\_ExpectationX@>
@<C\_get\_ExpectationInfluence@>
@<C\_get\_CovarianceInfluence@>
@<C\_get\_VarianceInfluence@>
@<C\_get\_Work@>
@<C\_get\_TableBlock@>
@<C\_get\_Sumweights@>
@<C\_get\_Table@>
@<C\_get\_dimTable@>
@<C\_get\_Lb@>
@<C\_get\_nperm@>
@<C\_get\_PermutedLinearStatistic@>
@<C\_get\_tol@>
@<RC\_init\_LECV\_1d@>
@<RC\_init\_LECV\_2d@>
@}

@d R LECV Input
@{
SEXP LECV
@|LECV
@}

@d C\_get\_P
@{
int C_get_P
(
@<R LECV Input@>
) {

    return(INTEGER(VECTOR_ELT(LECV, dim_SLOT))[0]);
}
@|C_get_P
@}

@d C\_get\_Q
@{
int C_get_Q
(
@<R LECV Input@>
) {

    return(INTEGER(VECTOR_ELT(LECV, dim_SLOT))[1]);
}
@|C_get_Q
@}

@d C\_get\_varonly
@{
int C_get_varonly
(
@<R LECV Input@>
) {

    return(INTEGER(VECTOR_ELT(LECV, varonly_SLOT))[0]);
}
@|C_get_varonly
@}

@d C\_get\_Xfactor
@{
int C_get_Xfactor
(
@<R LECV Input@>
) {

    return(INTEGER(VECTOR_ELT(LECV, Xfactor_SLOT))[0]);
}
@|C_get_Xfactor
@}

@d C\_get\_LinearStatistic
@{
double* C_get_LinearStatistic
(
@<R LECV Input@>
) {

    return(REAL(VECTOR_ELT(LECV, LinearStatistic_SLOT)));
}
@|C_get_LinearStatistic
@}

@d C\_get\_Expectation
@{
double* C_get_Expectation
(
@<R LECV Input@>
) {

    return(REAL(VECTOR_ELT(LECV, Expectation_SLOT)));
}
@|C_get_Expectation
@}

@d C\_get\_Variance
@{
double* C_get_Variance
(
@<R LECV Input@>
) {

    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    double *var, *covar;

    if (isNull(VECTOR_ELT(LECV, Variance_SLOT))) {
        SET_VECTOR_ELT(LECV, Variance_SLOT,
                       allocVector(REALSXP, PQ));
        if (!isNull(VECTOR_ELT(LECV, Covariance_SLOT))) {
            covar = REAL(VECTOR_ELT(LECV, Covariance_SLOT));
            var = REAL(VECTOR_ELT(LECV, Variance_SLOT));
            for (int p = 0; p < PQ; p++)
                var[p] = covar[S(p, p, PQ)];
        }
    }
    return(REAL(VECTOR_ELT(LECV, Variance_SLOT)));
}
@|C_get_Variance
@}

@d C\_get\_Covariance
@{
double* C_get_Covariance
(
@<R LECV Input@>
) {

    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract covariance from variance only object");
    if (C_get_varonly(LECV) && PQ == 1)
        return(C_get_Variance(LECV));
    return(REAL(VECTOR_ELT(LECV, Covariance_SLOT)));
}
@|C_get_Covariance
@}

@d C\_get\_MPinv
@{
double* C_get_MPinv
(
@<R LECV Input@>
) {

    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract MPinv from variance only object");
    /* allocate memory on as needed basis */
    if (isNull(VECTOR_ELT(LECV, MPinv_SLOT))) {
        SET_VECTOR_ELT(LECV, MPinv_SLOT,
                       allocVector(REALSXP,
                                   PQ * (PQ + 1) / 2));
    }
    return(REAL(VECTOR_ELT(LECV, MPinv_SLOT)));
}
@|C_get_MPinv
@}

@d C\_get\_ExpectationX
@{
double* C_get_ExpectationX
(
@<R LECV Input@>
) {
    return(REAL(VECTOR_ELT(LECV, ExpectationX_SLOT)));
}
@|C_get_ExpectationX
@}

@d C\_get\_ExpectationInfluence
@{
double* C_get_ExpectationInfluence
(
@<R LECV Input@>
) {

    return(REAL(VECTOR_ELT(LECV, ExpectationInfluence_SLOT)));
}
@|C_get_ExpectationInfluence
@}

@d C\_get\_CovarianceInfluence
@{
double* C_get_CovarianceInfluence
(
@<R LECV Input@>
) {

    return(REAL(VECTOR_ELT(LECV, CovarianceInfluence_SLOT)));
}
@|C_get_CovarianceInfluence
@}

@d C\_get\_VarianceInfluence
@{
double* C_get_VarianceInfluence
(
@<R LECV Input@>
) {

    return(REAL(VECTOR_ELT(LECV, VarianceInfluence_SLOT)));
}
@|C_get_VarianceInfluence
@}

@d C\_get\_Work
@{
double* C_get_Work
(
@<R LECV Input@>
) {

    return(REAL(VECTOR_ELT(LECV, Work_SLOT)));
}
@|C_get_Work
@}

@d C\_get\_TableBlock
@{
double* C_get_TableBlock
(
@<R LECV Input@>
) {

    if (VECTOR_ELT(LECV, TableBlock_SLOT) == R_NilValue)
        error("object does not contain table block slot");
    return(REAL(VECTOR_ELT(LECV, TableBlock_SLOT)));
}
@|C_get_TableBlock
@}

@d C\_get\_Sumweights
@{
double* C_get_Sumweights
(
@<R LECV Input@>
) {
    if (VECTOR_ELT(LECV, Sumweights_SLOT) == R_NilValue)
        error("object does not contain sumweights slot");
    return(REAL(VECTOR_ELT(LECV, Sumweights_SLOT)));
}
@|C_get_Sumweights
@}

@d C\_get\_Table
@{
double* C_get_Table
(
@<R LECV Input@>
) {

    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");
    return(REAL(VECTOR_ELT(LECV, Table_SLOT)));
}
@|C_get_Table
@}

@d C\_get\_dimTable
@{
int* C_get_dimTable
(
@<R LECV Input@>
) {

    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");
    return(INTEGER(getAttrib(VECTOR_ELT(LECV, Table_SLOT),
                             R_DimSymbol)));
}
@|C_get_dimTable
@}

@d C\_get\_Lb
@{
int C_get_Lb
(
@<R LECV Input@>
) {

    if (VECTOR_ELT(LECV, TableBlock_SLOT) != R_NilValue)
        return(LENGTH(VECTOR_ELT(LECV, Sumweights_SLOT)));
    return(C_get_dimTable(LECV)[2]);
}
@|C_get_Lb
@}

@d C\_get\_nperm
@{
R_xlen_t C_get_nperm
(
@<R LECV Input@>
) {
    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    return(XLENGTH(VECTOR_ELT(LECV, PermutedLinearStatistic_SLOT)) / PQ);
}
@|C_get_nperm
@}

@d C\_get\_PermutedLinearStatistic
@{
double* C_get_PermutedLinearStatistic
(
@<R LECV Input@>
) {

    return(REAL(VECTOR_ELT(LECV, PermutedLinearStatistic_SLOT)));
}
@|C_get_PermutedLinearStatistic
@}

@d C\_get\_tol
@{
double C_get_tol
(
@<R LECV Input@>
) {
    return(REAL(VECTOR_ELT(LECV, tol_SLOT))[0]);
}
@|C_get_tol
@}

@d Memory Input Checks
@{
if (P <= 0)
    error("P is not positive");

if (Q <= 0)
    error("Q is not positive");

if (Lb <= 0)
    error("Lb is not positive");

if (varonly < 0 || varonly > 1)
    error("varonly is not 0 or 1");

if (Xfactor < 0 || Xfactor > 1)
    error("Xfactor is not 0 or 1");

if (tol <= DBL_MIN)
    error("tol is not positive");
@}

@d Memory Names
@{
PROTECT(names = allocVector(STRSXP, Table_SLOT + 1));
SET_STRING_ELT(names, LinearStatistic_SLOT, mkChar("LinearStatistic"));
SET_STRING_ELT(names, Expectation_SLOT, mkChar("Expectation"));
SET_STRING_ELT(names, varonly_SLOT, mkChar("varonly"));
SET_STRING_ELT(names, Variance_SLOT, mkChar("Variance"));
SET_STRING_ELT(names, Covariance_SLOT, mkChar("Covariance"));
SET_STRING_ELT(names, MPinv_SLOT, mkChar("MPinv"));
SET_STRING_ELT(names, Work_SLOT, mkChar("Work"));
SET_STRING_ELT(names, ExpectationX_SLOT, mkChar("ExpectationX"));
SET_STRING_ELT(names, dim_SLOT, mkChar("dimension"));
SET_STRING_ELT(names, ExpectationInfluence_SLOT,
               mkChar("ExpectationInfluence"));
SET_STRING_ELT(names, Xfactor_SLOT, mkChar("Xfactor"));
SET_STRING_ELT(names, CovarianceInfluence_SLOT,
               mkChar("CovarianceInfluence"));
SET_STRING_ELT(names, VarianceInfluence_SLOT,
               mkChar("VarianceInfluence"));
SET_STRING_ELT(names, TableBlock_SLOT, mkChar("TableBlock"));
SET_STRING_ELT(names, Sumweights_SLOT, mkChar("Sumweights"));
SET_STRING_ELT(names, PermutedLinearStatistic_SLOT,
               mkChar("PermutedLinearStatistic"));
SET_STRING_ELT(names, tol_SLOT, mkChar("tol"));
SET_STRING_ELT(names, Table_SLOT, mkChar("Table"));
@}

@d R\_init\_LECV
@{
    SEXP vo, d, names, tolerance;
    int PQ; 

    @<Memory Input Checks@>

    PQ = P * Q;

    @<Memory Names@>

    /* Table_SLOT is always last and only used in 2d case
       ie omitted here */
    PROTECT(ans = allocVector(VECSXP, Table_SLOT + 1));
    SET_VECTOR_ELT(ans, LinearStatistic_SLOT, allocVector(REALSXP, PQ));
    SET_VECTOR_ELT(ans, Expectation_SLOT, allocVector(REALSXP, PQ));
    SET_VECTOR_ELT(ans, varonly_SLOT, vo = allocVector(INTSXP, 1));
    INTEGER(vo)[0] = varonly;
    if (varonly) {
        SET_VECTOR_ELT(ans, Variance_SLOT, allocVector(REALSXP, PQ));
    } else  {
        SET_VECTOR_ELT(ans, Covariance_SLOT,
                       allocVector(REALSXP, PQ * (PQ + 1) / 2));
    }
    SET_VECTOR_ELT(ans, ExpectationX_SLOT, allocVector(REALSXP, P));
    SET_VECTOR_ELT(ans, dim_SLOT, d = allocVector(INTSXP, 2));
    INTEGER(d)[0] = P;
    INTEGER(d)[1] = Q;
    SET_VECTOR_ELT(ans, ExpectationInfluence_SLOT,
                   allocVector(REALSXP, Lb * Q));

    /* should always _both_ be there */
    SET_VECTOR_ELT(ans, VarianceInfluence_SLOT,
                   allocVector(REALSXP, Lb * Q));
    SET_VECTOR_ELT(ans, CovarianceInfluence_SLOT,
                   allocVector(REALSXP, Lb * Q * (Q + 1) / 2));

    SET_VECTOR_ELT(ans, Xfactor_SLOT, allocVector(INTSXP, 1));
    INTEGER(VECTOR_ELT(ans, Xfactor_SLOT))[0] = Xfactor;
    SET_VECTOR_ELT(ans, TableBlock_SLOT, allocVector(REALSXP, Lb + 1));
    SET_VECTOR_ELT(ans, Sumweights_SLOT, allocVector(REALSXP, Lb));
    SET_VECTOR_ELT(ans, PermutedLinearStatistic_SLOT,
                   allocMatrix(REALSXP, 0, 0));
    SET_VECTOR_ELT(ans, tol_SLOT, tolerance = allocVector(REALSXP, 1));
    REAL(tolerance)[0] = tol;
    namesgets(ans, names);
@}

@d RC\_init\_LECV\_1d
@{
SEXP RC_init_LECV_1d
(
    @<C integer P Input@>,
    @<C integer Q Input@>,
    int varonly,
    @<C integer Lb Input@>,
    int Xfactor,
    double tol
) {

    SEXP ans;

    @<R\_init\_LECV@>

    if (varonly) {
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP, 3 * P + 1));
    } else  {
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP,
                           P + 2 * P * (P + 1) / 2 + 1));
    }

    SET_VECTOR_ELT(ans, TableBlock_SLOT,
                   allocVector(REALSXP, Lb + 1));

    SET_VECTOR_ELT(ans, Sumweights_SLOT,
                   allocVector(REALSXP, Lb));

    UNPROTECT(2);
    return(ans);
}
@|RC_init_LECV_1d
@}

@d RC\_init\_LECV\_2d
@{
SEXP RC_init_LECV_2d
(
    @<C integer P Input@>,
    @<C integer Q Input@>,
    int varonly,
    int Lx,
    int Ly,
    @<C integer Lb Input@>,
    int Xfactor,
    double tol
) {
    SEXP ans, tabdim, tab;

    if (Lx <= 0)
        error("Lx is not positive");

    if (Ly <= 0)
        error("Ly is not positive");

    @<R\_init\_LECV@>

    if (varonly) {
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP, 2 * P));
    } else  {
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP,
                           2 * P * (P + 1) / 2));
    }

    PROTECT(tabdim = allocVector(INTSXP, 3));
    INTEGER(tabdim)[0] = Lx + 1;
    INTEGER(tabdim)[1] = Ly + 1;
    INTEGER(tabdim)[2] = Lb;
    SET_VECTOR_ELT(ans, Table_SLOT,
                   tab = allocVector(REALSXP,
                       INTEGER(tabdim)[0] *
                       INTEGER(tabdim)[1] *
                       INTEGER(tabdim)[2]));
    dimgets(tab, tabdim);

    UNPROTECT(3);
    return(ans);
}
@|RC_init_LECV_2d
@}

@S

\chapter{R Code}

@s


\section{Example}

<<ex>>=
summary(1)
@@

@S

\end{document}
