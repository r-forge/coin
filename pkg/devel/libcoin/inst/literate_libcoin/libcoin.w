
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

#define DoSymmetric 1
#define DoCenter 1
#define DoVarOnly 1
#define Power1 1
#define Power2 2
#define Offset0 0
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
@<Function Definitions@>
@}

@d Function Definitions
@{
@<MoreUtils@>
@<Memory@>
@<KronSums@>
@<colSums@>
@<SimpleSums@>
@<Tables@>
@<Utils@>
@<LinearStatistics@>
@<ExpectationCovariances@>
@<User Interface@>
@}

The \proglang{R} interfaces are used to implement
regression tests to be called from within \proglang{R}

\section{Variables}

$N$ is the number of observations

@d R N Input
@{
    SEXP N,
@}

which at \proglang{C} level is represented as \verb|R_xlen_t| to allow for
$N > $ \verb|INT_MAX|

@d C integer N Input
@{
    R_xlen_t N
@}

The regressors $\x_i, i = 1, \dots, N$ 

@d R x Input
@{
    SEXP x,
@}

are either represented as a real matrix with $N$ rows and $P$ columns

@d C integer P Input
@{
    int P
@}

@d C real x Input
@{
    double *x,
    @<C integer N Input@>,
    @<C integer P Input@>,
@}

or as a factor (an integer at \proglang{C} level) at $P$ levels

@d C integer x Input
@{
    int *x,
    @<C integer N Input@>,
    @<C integer P Input@>,
@}

The influence functions are also either a $N \times Q$ real matrix

@d R y Input
@{
    SEXP y,
@}

@d C integer Q Input
@{
    int Q
@}

@d C real y Input
@{
    double *y,
    @<C integer Q Input@>,
@}

or a factor at $Q$ levels

@d C integer y Input
@{
    int *y,
    @<C integer Q Input@>,
@}

The weights $w_i, i = 1, \dots, N$ 

@d R weights Input
@{
    SEXP weights,
@}

can be constant one \verb|XLENGTH(weights) == 0| or integer-valued, with 
\verb|HAS_WEIGHTS == 0| in the former case

@d C integer weights Input
@{
    int *weights,
    int HAS_WEIGHTS,
@}

Weights larger than \verb|INT_MAX| are stored as double

@d C real weights Input
@{
    double *weights,
    int HAS_WEIGHTS,
@}

The sum of all weights is a double

@d C sumweights Input
@{
    double sumweights
@}

Subsets $\A \subseteq \{1, \dots, N\}$ 

@d R subset Input
@{
    SEXP subset
@}

are either not existant (\verb|XLENGTH(subset) == 0|) or of length

@d C integer Nsubset Input
@{
    R_xlen_t Nsubset
@}

Optionally, one can specify a subset of the subset via

@d C subset range Input
@{
    R_xlen_t offset,
    @<C integer Nsubset Input@>
@}

Subsets are stored either as integer

@d C integer subset Input
@{
    int *subset,
    @<C subset range Input@>
@}

or double (to allow for indices larger than \verb|INT_MAX|)

@d C real subset Input
@{
    double *subset,
    @<C subset range Input@>
@}

Blocks $b_i, i = 1, \dots, N$

@d R block Input
@{
    SEXP block,
@}

at $B$ levels

@d C integer Lb Input
@{
    int Lb
@}

are stored as a factor

@d C integer block Input
@{
    int *block,
    @<C integer Lb Input@>,
@}

The tabulation of $b$ (potentially in subsets) is

@d R blockTable Input
@{
    SEXP blockTable
@}

where the table is of length $B + 1$ and the first element
counts the number of missing values.


<<Sums-setup>>=
### replace with library("libcoin")
dyn.load("Sums.so")
set.seed(29)
N <- 20L
P <- 3L
Q <- 4L
L <- 2L
x <- matrix(runif(N * P), nrow = N)
y <- matrix(runif(N * Q), nrow = N)
ix <- sample(1:P, size = N, replace = TRUE)
iX <- diag(P)[ix,] ### model.matrix
iy <- sample(1:Q, size = N, replace = TRUE)
weights <- sample(0:5, size = N, replace = TRUE)
block <- sample(gl(L, ceiling(N / L))[1:N])
subset <- sort(sample(1:N, floor(N * 1.5), replace = TRUE))
subsety <- sample(1:N, floor(N * 1.5), replace = TRUE)
r1 <- rep(1:ncol(x), ncol(y))
r2 <- rep(1:ncol(y), each = ncol(x))
@@

\section{User Interface}

<<together>>=
library("libcoin")
(LECV <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset - 1L, 
               integer(0), 0L, 0.00001))
lcv <- LinStatExpCov(X = x, Y = y, weights = weights, subset = subset)
all.equal(LECV, lcv)
LECV$CovarianceInfluence
lcv$CovarianceInfluence

sw <- sum(weights[subset])
expecty <- colSums(y[subset, ] * weights[subset]) / sw
yc <- t(t(y) - expecty)
r1y <- rep(1:ncol(y), ncol(y))
r2y <- rep(1:ncol(y), each = ncol(y))
matrix(colSums(yc[subset, r1y] * yc[subset, r2y] * weights[subset]) / sw,
       ncol = Q)



V <- matrix(0, nrow = P * Q, ncol = P * Q)
V[lower.tri(V, diag = TRUE)] <- LECV$Covariance
LSvar <- diag(V)
V <- matrix(0, nrow = Q, ncol = Q)
V[lower.tri(V, diag = TRUE)] <- LECV$CovarianceInfluence
Ivar <- diag(V)
(LEV <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset - 1L, 
              integer(0), 1L, 0.00001))
lv <- LinStatExpCov(X = x, Y = y, weights = weights, subset = subset, varonly = TRUE)
all.equal(LEV, lv)
stopifnot(all.equal(LECV$LinearStatistic, LEV$LinearStatistic) &&
          all.equal(LECV$Expectation, LEV$Expectation) &&
          all.equal(LECV$ExpectationInfluence, LEV$ExpectationInfluence) &&
          all.equal(Ivar, LEV$VarianceInfluence) &&
          all.equal(LSvar, LEV$Variance))

library("libcoin")
LECVb <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset - 1L, 
               block, 0L, 0.00001)
lcvb <- LinStatExpCov(X = x, Y = y, weights = weights, subset = subset, block = block)
all.equal(LECVb, lcvb)

LEVb <- .Call("R_ExpectationCovarianceStatistic", x, y, weights, subset - 1L, 
               block, 1L, 0.00001)
lvb <- LinStatExpCov(X = x, Y = y, weights = weights, subset = subset, block = block, varonly =
TRUE)
all.equal(LEVb, lvb)


@@

@d User Interface
@{
@<RC\_ExpectationCovarianceStatistic@>
@<R\_ExpectationCovarianceStatistic@>
@}

@d User Interface Inputs
@{
@<R x Input@>
@<R y Input@>
@<R weights Input@>
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

    PROTECT(ans = R_init_LECV_1d(P, Q, varonly, Lb, Xfactor, tol));

    RC_ExpectationCovarianceStatistic(x, y, weights, subset, block, ans);

    UNPROTECT(5);
    return(ans);
}
@|R_ExpectationCovarianceStatistic
@}


@d Setup Dimensions
@{
    SEXP P, Q, Lb, Xfactor;

    PROTECT(P = ScalarInteger(0));
    PROTECT(Q = ScalarInteger(0));
    PROTECT(Lb = ScalarInteger(0));
    PROTECT(Xfactor = ScalarInteger(0));

    if (isInteger(x)) {
        INTEGER(P)[0] = NLEVELS(x);
        INTEGER(Xfactor)[0] = 1;
    } else {
        INTEGER(P)[0] = NCOL(x);
        INTEGER(Xfactor)[0] = 0;
    }
    INTEGER(Q)[0] = NCOL(y);

    INTEGER(Lb)[0] = 1;
    if (LENGTH(block) > 0)
        INTEGER(Lb)[0] = NLEVELS(block);
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


\section{Linear Statistics}

<<LinearStatistics>>=
a0 <- colSums(x[subset,r1] * y[subset,r2] * weights[subset])
a1 <- .Call("R_LinearStatistic", x, P, y, weights, subset - 1L, integer(0))
a2 <- .Call("R_LinearStatistic", x, P, y, as.double(weights), as.double(subset - 1L), integer(0))
a3 <- .Call("R_LinearStatistic", x, P, y, weights, as.double(subset - 1L), integer(0))
a4 <- .Call("R_LinearStatistic", x, P, y, as.double(weights), subset - 1L, integer(0))

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECV$LinearStatistic))


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
    @<R weights Input@>
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
@}

@d RC\_LinearStatistic Prototype
@{
void RC_LinearStatistic
(
    @<R x Input@>
    @<C integer N Input@>,
    @<C integer P Input@>,
    @<C real y Input@>
    @<R weights Input@>
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
sw <- sum(weights[subset])
expecty <- a0 <- colSums(y[subset, ] * weights[subset]) / sw
a1 <- .Call("R_ExpectationInfluence", y, weights, subset - 1L);
a2 <- .Call("R_ExpectationInfluence", y, as.double(weights), as.double(subset - 1L));
a3 <- .Call("R_ExpectationInfluence", y, weights, as.double(subset - 1L));
a4 <- .Call("R_ExpectationInfluence", y, as.double(weights), subset - 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECV$ExpectationInfluence))
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
    @<C integer Q Input@>;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;
    double sw;

    Q = NCOL(y);
    N = XLENGTH(y) / Q;
    Nsubset = XLENGTH(subset);

    sw = RC_Sums(N, weights, subset, Offset0, Nsubset);

    PROTECT(ans = allocVector(REALSXP, Q));
    RC_ExpectationInfluence(N, y, Q, weights, subset, Offset0, Nsubset, sw, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

@d RC\_ExpectationInfluence Prototype
@{
void RC_ExpectationInfluence
(
    @<C integer N Input@>,
    @<R y Input@>
    @<C integer Q Input@>,
    @<R weights Input@>
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
sw <- sum(weights[subset])
yc <- t(t(y) - expecty)
r1y <- rep(1:ncol(y), ncol(y))
r2y <- rep(1:ncol(y), each = ncol(y))
a0 <- colSums(yc[subset, r1y] * yc[subset, r2y] * weights[subset]) / sw
a0 <- matrix(a0, ncol = ncol(y))
vary <- diag(a0)
a0 <- a0[lower.tri(a0, diag = TRUE)]
a1 <- .Call("R_CovarianceInfluence", y, weights, subset - 1L, 0L);
a2 <- .Call("R_CovarianceInfluence", y, as.double(weights), as.double(subset - 1L), 0L);
a3 <- .Call("R_CovarianceInfluence", y, weights, as.double(subset - 1L), 0L);
a4 <- .Call("R_CovarianceInfluence", y, as.double(weights), subset - 1L, 0L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECV$CovarianceInfluence) &&
          all.equal(vary, LEV$VarianceInfluence))

a1 <- .Call("R_CovarianceInfluence", y, weights, subset - 1L, 1L);
a2 <- .Call("R_CovarianceInfluence", y, as.double(weights), as.double(subset - 1L), 1L);
a3 <- .Call("R_CovarianceInfluence", y, weights, as.double(subset - 1L), 1L);
a4 <- .Call("R_CovarianceInfluence", y, as.double(weights), subset - 1L, 1L);

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
    @<R weights Input@>
    @<R subset Input@>,
    SEXP varonly
) 
{
    SEXP ans;
    SEXP ExpInf;
    @<C integer Q Input@>;
    @<C integer N Input@>;
    @<C integer Nsubset Input@>;
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
    @<C integer N Input@>,
    @<R y Input@>
    @<C integer Q Input@>,
    @<R weights Input@>
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
a1 <- .Call("R_ExpectationX", x, P, weights, subset - 1L);
a2 <- .Call("R_ExpectationX", x, P, as.double(weights), as.double(subset - 1L));
a3 <- .Call("R_ExpectationX", x, P, weights, as.double(subset - 1L));
a4 <- .Call("R_ExpectationX", x, P, as.double(weights), subset - 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4) &&
          all.equal(a0, LECV$ExpectationX))

a0 <- colSums(x[subset, ]^2 * weights[subset]) 
a1 <- .Call("R_CovarianceX", x, P, weights, subset - 1L, 1L);
a2 <- .Call("R_CovarianceX", x, P, as.double(weights), as.double(subset - 1L), 1L);
a3 <- .Call("R_CovarianceX", x, P, weights, as.double(subset - 1L), 1L);
a4 <- .Call("R_CovarianceX", x, P, as.double(weights), subset - 1L, 1L);

stopifnot(all.equal(a0, a1) && all.equal(a0, a2) &&
          all.equal(a0, a3) && all.equal(a0, a4))

a0 <- as.vector(colSums(iX[subset, ] * weights[subset]))
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
@}

@d RC\_ExpectationX Prototype
@{
void RC_ExpectationX
(
    @<R x Input@>
    @<C integer N Input@>,
    @<C integer P Input@>,
    @<R weights Input@>
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
        RC_OneTableSums(INTEGER(x), N, P, weights, subset, offset, Nsubset, P_ans);
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
    @<R weights Input@>
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
@}

@d RC\_CovarianceX Prototype
@{
void RC_CovarianceX
(
    @<R x Input@>
    @<C integer N Input@>,
    @<C integer P Input@>,
    @<R weights Input@>
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
    if (Nsubset > 0)
        diff = (R_xlen_t) s[0];
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
        diff = (R_xlen_t) s[1] - s[0];
        s++;
    } else {
        diff = 1;
    }
@}


\subsection{Simple Sums}

<<SimpleSums>>=
a0 <- sum(weights[subset])
a1 <- .Call("R_Sums", N, weights, subset - 1L)
a2 <- .Call("R_Sums", N, as.double(weights), as.double(subset - 1L))
a3 <- .Call("R_Sums", N, weights, as.double(subset - 1L))
a4 <- .Call("R_Sums", N, as.double(weights), subset - 1L)
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
    @<R weights Input@>
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
@}

@d RC\_Sums Prototype
@{
double RC_Sums
(
    @<C integer N Input@>,
    @<R weights Input@>
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
a1 <- .Call("R_KronSums", x, P, y, weights, subset - 1L)
a2 <- .Call("R_KronSums", x, P, y, as.double(weights), as.double(subset - 1L))
a3 <- .Call("R_KronSums", x, P, y, weights, as.double(subset - 1L))
a4 <- .Call("R_KronSums", x, P, y, as.double(weights), subset - 1L)

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
    @<R weights Input@>
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
@}

@d RC\_KronSums Prototype
@{
void RC_KronSums
(
    @<RC KronSums Input@>
    @<R weights Input@>
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
//    for (int q = 0; q < Q; q++) {
//        for (int p = 0; p < (SYMMETRIC ? q + 1 : P); p++) {

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
//            PQ_ans++;
        }
    }
@}

\subsubsection{Xfactor Kronecker Sums}

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


\subsubsection{Permuted Kronecker Sums}

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

\subsubsection{Xfactor Permuted Kronecker Sums}

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
    @<R weights Input@>
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
@}

@d RC\_colSums Prototype
@{
void RC_colSums
(
    @<C colSums Input@>
    @<R weights Input@>
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
    @<R weights Input@>
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
@}


@d RC\_OneTableSums Prototype
@{
void RC_OneTableSums
(
    @<C OneTableSums Input@>
    @<R weights Input@>
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
@}


@d RC\_TwoTableSums Prototype
@{
void RC_TwoTableSums
(
    @<C TwoTableSums Input@>
    @<R weights Input@>
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
@}


@d RC\_ThreeTableSums Prototype
@{
void RC_ThreeTableSums
(
    @<C ThreeTableSums Input@>
    @<R weights Input@>
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

@d Utils
@{
@<C\_setup\_subset@>
@<C\_setup\_subset\_block@>
@<C\_order\_subset\_wrt\_block@>
@<RC\_order\_subset\_wrt\_block@>;
@<R\_order\_subset\_wrt\_block@>;
@}

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

@d MoreUtils
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

/**
    Computes the Kronecker product of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

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
@}

\section{Memory}

@d Memory
@{
int C_get_P
(
    SEXP LECV
) {

    return(INTEGER(VECTOR_ELT(LECV, dim_SLOT))[0]);
}

int C_get_Q
(
    SEXP LECV
) {

    return(INTEGER(VECTOR_ELT(LECV, dim_SLOT))[1]);
}

int C_get_varonly
(
    SEXP LECV
) {

    return(INTEGER(VECTOR_ELT(LECV, varonly_SLOT))[0]);
}

int C_get_Xfactor
(
    SEXP LECV
) {

    return(INTEGER(VECTOR_ELT(LECV, Xfactor_SLOT))[0]);
}

double* C_get_LinearStatistic
(
    SEXP LECV
) {

    return(REAL(VECTOR_ELT(LECV, LinearStatistic_SLOT)));
}

double* C_get_Expectation
(
    SEXP LECV
) {

    return(REAL(VECTOR_ELT(LECV, Expectation_SLOT)));
}

double* C_get_Variance
(
    SEXP LECV
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

double* C_get_Covariance
(
    SEXP LECV
) {

    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract covariance from variance only object");
    if (C_get_varonly(LECV) && PQ == 1)
        return(C_get_Variance(LECV));
    return(REAL(VECTOR_ELT(LECV, Covariance_SLOT)));
}

double* C_get_MPinv
(
    SEXP LECV
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


double* C_get_ExpectationX
(
    SEXP LECV
) {
    return(REAL(VECTOR_ELT(LECV, ExpectationX_SLOT)));
}

double* C_get_ExpectationInfluence
(
    SEXP LECV
) {

    return(REAL(VECTOR_ELT(LECV, ExpectationInfluence_SLOT)));
}

double* C_get_CovarianceInfluence
(
    SEXP LECV
) {

    return(REAL(VECTOR_ELT(LECV, CovarianceInfluence_SLOT)));
}

double* C_get_VarianceInfluence
(
    SEXP LECV
) {

    return(REAL(VECTOR_ELT(LECV, VarianceInfluence_SLOT)));
}

double* C_get_Work
(
    SEXP LECV
) {

    return(REAL(VECTOR_ELT(LECV, Work_SLOT)));
}

double* C_get_TableBlock
(
    SEXP LECV
) {

    if (VECTOR_ELT(LECV, TableBlock_SLOT) == R_NilValue)
        error("object does not contain table block slot");
    return(REAL(VECTOR_ELT(LECV, TableBlock_SLOT)));
}

double* C_get_Sumweights
(
    SEXP LECV
) {
    if (VECTOR_ELT(LECV, Sumweights_SLOT) == R_NilValue)
        error("object does not contain sumweights slot");
    return(REAL(VECTOR_ELT(LECV, Sumweights_SLOT)));
}

int* C_get_Table
(
    SEXP LECV
) {

    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");
    return(INTEGER(VECTOR_ELT(LECV, Table_SLOT)));
}

int* C_get_dimTable
(
    SEXP LECV
) {

    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");
    return(INTEGER(getAttrib(VECTOR_ELT(LECV, Table_SLOT),
                             R_DimSymbol)));
}


int C_get_Lb
(
    SEXP LECV
) {

    if (VECTOR_ELT(LECV, TableBlock_SLOT) != R_NilValue)
        return(LENGTH(VECTOR_ELT(LECV, Sumweights_SLOT)));
    return(C_get_dimTable(LECV)[2]);
}

int C_get_B
(
    SEXP LECV
) {

    return(NCOL(VECTOR_ELT(LECV, PermutedLinearStatistic_SLOT)));
}

double* C_get_PermutedLinearStatistic
(
    SEXP LECV
) {

    return(REAL(VECTOR_ELT(LECV, PermutedLinearStatistic_SLOT)));
}

double C_get_tol
(
    SEXP LECV
) {

    return(REAL(VECTOR_ELT(LECV, tol_SLOT))[0]);
}

SEXP R_init_LECV
(
    SEXP P,
    SEXP Q,
    SEXP varonly,
    SEXP Lb,
    SEXP Xfactor,
    SEXP tol
) {

    SEXP ans, vo, d, names, tolerance;
    int p, q, pq, lb;

    if (!isInteger(P) || LENGTH(P) != 1)
        error("P is not a scalar integer");
    if (INTEGER(P)[0] <= 0)
        error("P is not positive");

    if (!isInteger(Q) || LENGTH(Q) != 1)
        error("Q is not a scalar integer");
    if (INTEGER(Q)[0] <= 0)
        error("Q is not positive");

    if (!isInteger(Lb) || LENGTH(Lb) != 1)
        error("Lb is not a scalar integer");
    if (INTEGER(Lb)[0] <= 0)
        error("Lb is not positive");

    if (!isInteger(varonly) || LENGTH(varonly) != 1)
        error("varonly is not a scalar integer");
    if (INTEGER(varonly)[0] < 0 || INTEGER(varonly)[0] > 1)
            error("varonly is not 0 or 1");

    if (!isReal(tol) || LENGTH(tol) != 1)
        error("tol is not a scalar double");
    if (REAL(tol)[0] <= DBL_MIN)
            error("tol is not positive");

    p = INTEGER(P)[0];
    q = INTEGER(Q)[0];
    lb = INTEGER(Lb)[0];
    pq = p * q;

    /* Table_SLOT is always last and only used in 2d case
       ie omitted here */
    PROTECT(ans = allocVector(VECSXP, Table_SLOT + 1));
    PROTECT(names = allocVector(STRSXP, Table_SLOT + 1));

    SET_VECTOR_ELT(ans, LinearStatistic_SLOT,
                   allocVector(REALSXP, pq));
    SET_STRING_ELT(names, LinearStatistic_SLOT,
                   mkChar("LinearStatistic"));

    SET_VECTOR_ELT(ans, Expectation_SLOT,
                   allocVector(REALSXP, pq));
    SET_STRING_ELT(names, Expectation_SLOT,
                   mkChar("Expectation"));

    SET_VECTOR_ELT(ans, varonly_SLOT,
                   vo = allocVector(INTSXP, 1));
    SET_STRING_ELT(names, varonly_SLOT,
                   mkChar("varonly"));

    INTEGER(vo)[0] = INTEGER(varonly)[0];
    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, Variance_SLOT,
                       allocVector(REALSXP, pq));
    } else  {
        SET_VECTOR_ELT(ans, Covariance_SLOT,
                       allocVector(REALSXP,
                                   pq * (pq + 1) / 2));
    }

    SET_STRING_ELT(names, Variance_SLOT,
                   mkChar("Variance"));
    SET_STRING_ELT(names, Covariance_SLOT,
                   mkChar("Covariance"));
    SET_STRING_ELT(names, MPinv_SLOT,
                   mkChar("MPinv"));
    SET_STRING_ELT(names, Work_SLOT,
                   mkChar("Work"));

    SET_VECTOR_ELT(ans, ExpectationX_SLOT,
                   allocVector(REALSXP, p));
    SET_STRING_ELT(names, ExpectationX_SLOT,
                   mkChar("ExpectationX"));

    SET_VECTOR_ELT(ans, dim_SLOT,
                   d = allocVector(INTSXP, 2));
    SET_STRING_ELT(names, dim_SLOT,
                   mkChar("dimension"));
    INTEGER(d)[0] = p;
    INTEGER(d)[1] = q;

    SET_VECTOR_ELT(ans, ExpectationInfluence_SLOT,
                   allocVector(REALSXP, lb * q));
    SET_STRING_ELT(names, ExpectationInfluence_SLOT,
                   mkChar("ExpectationInfluence"));

    /* should always _both_ be there */
    SET_VECTOR_ELT(ans, VarianceInfluence_SLOT,
                   allocVector(REALSXP, lb * q));
    SET_STRING_ELT(names, VarianceInfluence_SLOT,
                   mkChar("VarianceInfluence"));
    SET_VECTOR_ELT(ans, CovarianceInfluence_SLOT,
                   allocVector(REALSXP, lb * q * (q + 1) / 2));
    SET_STRING_ELT(names, CovarianceInfluence_SLOT,
                   mkChar("CovarianceInfluence"));

    SET_VECTOR_ELT(ans, Xfactor_SLOT,
                   allocVector(INTSXP, 1));
    SET_STRING_ELT(names, Xfactor_SLOT,
                   mkChar("Xfactor"));
    INTEGER(VECTOR_ELT(ans, Xfactor_SLOT))[0] = INTEGER(Xfactor)[0];

    SET_VECTOR_ELT(ans, TableBlock_SLOT,
                   allocVector(REALSXP, lb + 1));
    SET_STRING_ELT(names, TableBlock_SLOT,
                   mkChar("TableBlock"));

    SET_VECTOR_ELT(ans, Sumweights_SLOT,
                   allocVector(REALSXP, lb));
    SET_STRING_ELT(names, Sumweights_SLOT,
                   mkChar("Sumweights"));

    SET_VECTOR_ELT(ans, PermutedLinearStatistic_SLOT,
                   allocMatrix(REALSXP, 0, 0));
    SET_STRING_ELT(names, PermutedLinearStatistic_SLOT,
                   mkChar("PermutedLinearStatistic"));

    SET_VECTOR_ELT(ans, tol_SLOT,
                   tolerance = allocVector(REALSXP, 1));
    SET_STRING_ELT(names, tol_SLOT,
                   mkChar("tol"));
    REAL(tolerance)[0] = REAL(tol)[0];

    SET_STRING_ELT(names, Table_SLOT,
                   mkChar("Table"));

    namesgets(ans, names);

    UNPROTECT(2);
    return(ans);
}


SEXP R_init_LECV_1d
(
    SEXP P,
    SEXP Q,
    SEXP varonly,
    SEXP Lb,
    SEXP Xfactor,
    SEXP tol
) {

    SEXP ans;
    int p, lb;

    p = INTEGER(P)[0];
    lb = INTEGER(Lb)[0];

    PROTECT(ans = R_init_LECV(P, Q, varonly, Lb, Xfactor, tol));

    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP, 3 * p + 1));
    } else  {
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP,
                           p + 2 * p * (p + 1) / 2 + 1));
    }

    SET_VECTOR_ELT(ans, TableBlock_SLOT,
                   allocVector(REALSXP, lb + 1));

    SET_VECTOR_ELT(ans, Sumweights_SLOT,
                   allocVector(REALSXP, lb));

    UNPROTECT(1);
    return(ans);
}

SEXP R_init_LECV_2d
(
    SEXP P,
    SEXP Q,
    SEXP varonly,
    SEXP Lx,
    SEXP Ly,
    SEXP Lb,
    SEXP Xfactor,
    SEXP tol
) {

    SEXP ans, tabdim, tab;
    int p, lb;

    if (!isInteger(Lx) || LENGTH(Lx) != 1)
        error("Lx is not a scalar integer");
    if (INTEGER(Lx)[0] <= 0)
        error("Lx is not positive");

    if (!isInteger(Ly) || LENGTH(Ly) != 1)
        error("Ly is not a scalar integer");
    if (INTEGER(Ly)[0] <= 0)
        error("Ly is not positive");

    p = INTEGER(P)[0];
    lb = INTEGER(Lb)[0];

    PROTECT(ans = R_init_LECV(P, Q, varonly, Lb, Xfactor, tol));

    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP, 2 * p));
    } else  {
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP,
                           2 * p * (p + 1) / 2));
    }

    PROTECT(tabdim = allocVector(INTSXP, 3));
    INTEGER(tabdim)[0] = INTEGER(Lx)[0] + 1;
    INTEGER(tabdim)[1] = INTEGER(Ly)[0] + 1;
    INTEGER(tabdim)[2] = lb;
    SET_VECTOR_ELT(ans, Table_SLOT,
                   tab = allocVector(REALSXP,
                       INTEGER(tabdim)[0] *
                       INTEGER(tabdim)[1] *
                       INTEGER(tabdim)[2]));
    dimgets(tab, tabdim);

    UNPROTECT(2);
    return(ans);
}
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
