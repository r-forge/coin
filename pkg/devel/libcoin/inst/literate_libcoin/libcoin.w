
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

The \proglang{C} file \verb|Sums.c| defines the \proglang{C}
functions and a corresponding \proglang{R} interface (via \verb|.C()|)

@o Sums.c -cc
@{
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
@<Function prototypes@>
@}

@d Function prototypes
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
@<RC\_KronSums\_Permutation@>
@<R\_KronSums\_Permutation@>
@<C\_TableSums\_dweights\_dsubset@>
@<C\_TableSums\_iweights\_dsubset@>
@<C\_TableSums\_iweights\_isubset@>
@<C\_TableSums\_dweights\_isubset@>
@<RC\_TableSums@>
@<R\_TableSums@>
@}


The \proglang{R} interfaces are used to implement
regression tests to be called from within \proglang{R}


\subsection{Sums}

<<regression-test-Sums>>=
### replace with library("libcoin")
dyn.load("Sums.so")
set.seed(29)
N <- 20L
P <- 3L
x <- matrix(runif(N * P), nrow = N)
weights <- sample(0:5, size = N, replace = TRUE)
subset <- sort(sample(1:N, floor(N/2)))
subsety <- sort(sample(1:N, floor(N/2)))

a0 <- sum(weights[subset])

a1 <- .Call("R_Sums", x, weights, subset - 1L)

a2 <- .Call("R_Sums", x, as.double(weights), as.double(subset - 1L))

max(a0 - a1)
max(a1 - a2)
@@



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


@d Sums Body
@{
    R_xlen_t diff;
    double ans = 0.0;

    s = subset;
    w = weights;
    if (Nsubset > 0) {
        diff = (R_xlen_t) s[offset];
    } else {
        diff = 0;
    }
    for (R_xlen_t i = offset; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) {
        w = w + diff;
        ans += w[0];
        if (Nsubset > 0) {
            diff = (R_xlen_t) s[1] - s[0];
            s++;
        } else {
            diff = 1;
        }
    }
    w = w + diff;
    ans += w[0];
    return(ans);
@}

@d Sums Inputs
@{
    int N,
@}

@d C\_Sums\_dweights\_dsubset
@{
double C_Sums_dweights_dsubset
(
    @<Sums Inputs@>
    double *weights,
    @<C real subset Input@>
) {

    double *s, *w; 

  @<Sums Body@>
}
@}

@d C\_Sums\_iweights\_dsubset
@{
double C_Sums_iweights_dsubset
(
    @<Sums Inputs@>
    int *weights,
    @<C real subset Input@>
) {

    double *s;
    int *w; 

  @<Sums Body@>
}
@}

@d C\_Sums\_iweights\_isubset
@{
double C_Sums_iweights_isubset
(
    @<Sums Inputs@>
    int *weights,
    @<C integer subset Input@>
) {

    int *s, *w;

  @<Sums Body@>
}
@}

@d C\_Sums\_dweights\_isubset
@{
double C_Sums_dweights_isubset
(
    @<Sums Inputs@>
    double *weights,
    @<C integer subset Input@>
) {

    int *s; 
    double *w;

  @<Sums Body@>
}
@}

@d RC\_Sums
@{

double RC_Sums
(
    R_xlen_t N,
    SEXP weights,
    SEXP subset,
    R_xlen_t offset,
    R_xlen_t Nsubset
) {
    
    if (XLENGTH(weights) == 0) {
        if (XLENGTH(subset) == 0) {
            return((double) N);
        } else {
            return((double) Nsubset - offset);
        }
    } 
    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            return(C_Sums_iweights_isubset(N, INTEGER(weights), INTEGER(subset), offset, Nsubset));
        } else {
            return(C_Sums_iweights_dsubset(N, INTEGER(weights), REAL(subset), offset, Nsubset));
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            return(C_Sums_dweights_isubset(N, REAL(weights), INTEGER(subset), offset, Nsubset));
        } else {
            return(C_Sums_dweights_dsubset(N, REAL(weights), REAL(subset), offset, Nsubset));
        }
    }
}

@}

@d R\_Sums
@{
SEXP R_Sums
(
    SEXP x,
    SEXP weights,
    SEXP subset
) {

    SEXP ans;
    int P;
    R_xlen_t N, Nsubset;

    P = NCOL(x);
    N = XLENGTH(x) / P;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = RC_Sums(N, weights, subset, 0, Nsubset);
    UNPROTECT(1);
    return(ans);
}
@}

\subsection{Kronecker Sums}

<<regression-test-KronSum>>=
### replace with library("libcoin")
dyn.load("Sums.so")
set.seed(29)
N <- 20L
P <- 3L
Q <- 2L
x <- matrix(runif(N * P), nrow = N)
y <- matrix(runif(N * Q), nrow = N)
weights <- sample(0:5, size = N, replace = TRUE)
subset <- sort(sample(1:N, floor(N/2)))
r1 <- rep(1:ncol(x), ncol(y))
r2 <- rep(1:ncol(y), each = ncol(x))

a0 <- colSums(x[subset,r1] * y[subset,r2] * weights[subset])

a1 <- .Call("R_KronSums", x, y, weights, subset - 1L)

a2 <- .Call("R_KronSums", x, y, as.double(weights), as.double(subset - 1L))

max(a0 - a1)
max(a1 - a2)
@@


@d KronSums Body
@{
    R_xlen_t diff;
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
            s = subset;
            w = weights;
            if (Nsubset > 0) {
                diff = (R_xlen_t) s[offset];
            } else {
                diff = 0;
            }
            for (R_xlen_t i = offset; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) {
                xx = xx + diff;
                yy = yy + diff;
                if (HAS_WEIGHTS) {
                    w = w + diff;
                    PQ_ans[0] += (xx[0] - cx) * (yy[0] - cy) * w[0];
                } else {
                    PQ_ans[0] += (xx[0] - cx) * (yy[0] - cy);
                }
                if (Nsubset > 0) {
                    diff = (R_xlen_t) s[1] - s[0];
                    s++;
                } else {
                    diff = 1;
                }
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


@d C x Input
@{
    double *x,
    const R_xlen_t N,
    const int P,
@}

@d C y Input
@{
    double *y,
    const int Q,
@}

@d KronSums Inputs
@{
    @<C x Input@>
    @<C y Input@>
    const int SYMMETRIC,
    double *centerx,
    double *centery,
    const int CENTER,
@}

@d KronSums Answer
@{
    double *PQ_ans
@}



@d C\_KronSums\_dweights\_dsubset
@{
void C_KronSums_dweights_dsubset
(
    @<KronSums Inputs@>
    double *weights,
    const int HAS_WEIGHTS,
    @<C real subset Input@>,
    @<KronSums Answer@>
) {

    double *s, *w; 

  @<KronSums Body@>
}
@}

@d C\_KronSums\_iweights\_dsubset
@{
void C_KronSums_iweights_dsubset
(
    @<KronSums Inputs@>
    int *weights,
    const int HAS_WEIGHTS,
    @<C real subset Input@>,
    @<KronSums Answer@>
) {

    double *s;
    int *w; 

  @<KronSums Body@>
}
@}

@d C\_KronSums\_iweights\_isubset
@{
void C_KronSums_iweights_isubset
(
    @<KronSums Inputs@>
    int *weights,
    const int HAS_WEIGHTS,
    @<C integer subset Input@>,
    @<KronSums Answer@>
) {

    int *s, *w;

  @<KronSums Body@>
}
@}

@d C\_KronSums\_dweights\_isubset
@{
void C_KronSums_dweights_isubset
(
    @<KronSums Inputs@>
    double *weights,
    const int HAS_WEIGHTS,
    @<C integer subset Input@>,
    @<KronSums Answer@>
) {

    int *s; 
    double *w;

  @<KronSums Body@>
}
@}

@d RC\_KronSums
@{

void RC_KronSums
(
    @<KronSums Inputs@>
    SEXP weights,
    SEXP subset,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<KronSums Answer@>
) {

    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_KronSums_iweights_isubset(x, N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_KronSums_iweights_dsubset(x, N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_KronSums_dweights_isubset(x, N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_KronSums_dweights_dsubset(x, N, P, y, Q, SYMMETRIC, centerx, centery, CENTER,
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    }
}

@}

@d R\_KronSums
@{
SEXP R_KronSums
(
    SEXP x,
    SEXP y,
    SEXP weights,
    SEXP subset
) {

    SEXP ans;
    int P, Q;
    R_xlen_t N, Nsubset;
    double *center;

    P = NCOL(x);
    Q = NCOL(y);
    N = XLENGTH(x) / P;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, P * Q));
    RC_KronSums(REAL(x), N, P, REAL(y), Q, 0, center, center, 0, weights, subset, 0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

<<regression-test-KronSums-Permutation>>=

<<regression-Permutation>>=
a0 <- colSums(x[subset,r1] * y[subsety, r2])
a1 <- .Call("R_KronSums_Permutation", x, y, subset - 1L, subsety -1L)
max(a0 - a1)
@@

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
            PQ_ans[0] += xx[0] * y[ (R_xlen_t) sy[0] + q *N];
        }
        PQ_ans++;
    }
@}


@d C\_KronSums\_Permutation\_dsubset
@{
void C_KronSums_Permutation_dsubset
(
    @<C x Input@>
    @<C y Input@>
    @<C real subset Input@>,
    double *subsety,
    @<KronSums Answer@>
) {

    double *sx, *sy;

  @<KronSums Permutation Body@>
}
@}

@d C\_KronSums\_Permutation\_isubset
@{
void C_KronSums_Permutation_isubset
(
    @<C x Input@>
    @<C y Input@>
    @<C integer subset Input@>,
    int *subsety,
    @<KronSums Answer@>
) {

    int *sx, *sy;

  @<KronSums Permutation Body@>
}
@}

@d RC\_KronSums\_Permutation
@{

void RC_KronSums_Permutation
(
    @<C x Input@>
    @<C y Input@>
    SEXP subset,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    SEXP subsety,
    @<KronSums Answer@>
) {

    if (TYPEOF(subset) == INTSXP) {
        C_KronSums_Permutation_isubset(x, N, P, y, Q, 
                                       INTEGER(subset), offset, Nsubset, INTEGER(subsety), PQ_ans);
        } else {
        C_KronSums_Permutation_dsubset(x, N, P, y, Q, 
                                       REAL(subset), offset, Nsubset, REAL(subsety), PQ_ans);
    }
}

@}

@d R\_KronSums\_Permutation
@{
SEXP R_KronSums_Permutation
(
    SEXP x,
    SEXP y,
    SEXP subset,
    SEXP subsety
) {

    SEXP ans;
    int P, Q;
    R_xlen_t N, Nsubset;
    double *center;

    P = NCOL(x);
    Q = NCOL(y);
    N = XLENGTH(x) / P;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, P * Q));
    RC_KronSums_Permutation(REAL(x), N, P, REAL(y), Q, subset, 0, Nsubset, subsety, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}


\subsection{Column Sums}

<<regression-test-colSums>>=
a0 <- colSums(x[subset,r1] * weights[subset])
a1 <- .Call("R_colSums", x, weights, subset - 1L)

a2 <- .Call("R_colSums", x, as.double(weights), as.double(subset - 1L))

max(a0 - a1)
max(a1 - a2)
@@


@d colSums Body
@{
    R_xlen_t diff;
    double *xx, cx = 0.0;

    for (int p = 0; p < P; p++) {
        P_ans[0] = 0.0;
        xx = x + N * p;
        if (CENTER) {
            cx = centerx[p];
        }
        s = subset;
        w = weights;
        if (Nsubset > 0) {
            diff = (R_xlen_t) s[offset];
        } else {
            diff = 0;
        }
        for (R_xlen_t i = offset; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) {
            xx = xx + diff;
            if (HAS_WEIGHTS) {
                w = w + diff;
                P_ans[0] += pow(xx[0] - cx, power) * w[0];
            } else {
                P_ans[0] += pow(xx[0] - cx, power);
            }
            if (Nsubset > 0) {
                diff = (R_xlen_t) s[1] - s[0];
                s++;
            } else {
                diff = 1;
            }
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

@d colSums Inputs
@{
    @<C x Input@>
    const int power,
    double *centerx,
    const int CENTER,
@}

@d colSums Answer
@{
    double *P_ans
@}

@d C\_colSums\_dweights\_dsubset
@{
void C_colSums_dweights_dsubset
(
    @<colSums Inputs@>
    double *weights,
    const int HAS_WEIGHTS,
    @<C real subset Input@>,
    @<colSums Answer@>
) {

    double *s, *w; 

  @<colSums Body@>
}
@}

@d C\_colSums\_iweights\_dsubset
@{
void C_colSums_iweights_dsubset
(
    @<colSums Inputs@>
    int *weights,
    const int HAS_WEIGHTS,
    @<C real subset Input@>,
    @<colSums Answer@>
) {

    double *s;
    int *w; 

  @<colSums Body@>
}
@}

@d C\_colSums\_iweights\_isubset
@{
void C_colSums_iweights_isubset
(
    @<colSums Inputs@>
    int *weights,
    const int HAS_WEIGHTS,
    @<C integer subset Input@>,
    @<colSums Answer@>
) {

    int *s, *w;

  @<colSums Body@>
}
@}

@d C\_colSums\_dweights\_isubset
@{
void C_colSums_dweights_isubset
(
    @<colSums Inputs@>
    double *weights,
    const int HAS_WEIGHTS,
    @<C integer subset Input@>,
    @<colSums Answer@>
) {

    int *s; 
    double *w;

  @<colSums Body@>
}
@}

@d RC\_colSums
@{

void RC_colSums
(
    @<colSums Inputs@>
    SEXP weights,
    SEXP subset,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<colSums Answer@>
) {

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

@}

@d R\_colSums
@{
SEXP R_colSums
(
    SEXP x,
    SEXP weights,
    SEXP subset
) {

    SEXP ans;
    int P;
    R_xlen_t N, Nsubset;
    double *center;

    P = NCOL(x);
    N = XLENGTH(x) / P;
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, P));
    RC_colSums(REAL(x), N, P, 1, center, 0, weights, subset, 0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
}
@}

\subsection{Table Sums}

<<regression-test-TableSum>>=
### replace with library("libcoin")
dyn.load("Sums.so")
set.seed(29)
N <- 20L
P <- 3L
Q <- 2L
x <- sample(1:P, size = N, replace = TRUE)
y <- sample(1:Q, size = N, replace = TRUE)
weights <- sample(0:5, size = N, replace = TRUE)
subset <- sort(sample(1:N, floor(N/2)))

a0 <- xtabs(weights ~ x + y, data = data.frame(x = x, y = y, weights =
weights)[subset,])

a1 <- .Call("R_TableSums", x, P + 1L, y, Q + 1L, weights, subset - 1L)
a1 <- matrix(a1, nrow = P + 1, ncol = Q + 1)[-1,-1]

a2 <- .Call("R_TableSums", x, P + 1L, y, Q + 1L, as.double(weights), as.double(subset - 1L))
a2 <- matrix(a2, nrow = P + 1, ncol = Q + 1)[-1,-1]

max(a0 - a1)
max(a1 - a2)
@@


@d TableSums Body
@{
    R_xlen_t diff;
    int *xx, *yy;

    for (int p = 0; p < Q * P; p++) PQ_ans[p] = 0.0;

    yy = y;
    xx = x;
    s = subset;
    w = weights;
    if (Nsubset > 0) {
        diff = (R_xlen_t) s[offset];
    } else {
        diff = 0;
    }
    for (R_xlen_t i = offset; i < (Nsubset == 0 ? N : Nsubset) - 1; i++) {
        xx = xx + diff;
        yy = yy + diff;
        if (HAS_WEIGHTS) {
            w = w + diff;
            PQ_ans[yy[0] * P + xx[0]] += (double) w[0];
        } else {
            PQ_ans[yy[0] * P + xx[0]]++;
        }
        if (Nsubset > 0) {
            diff = (R_xlen_t) s[1] - s[0];
            s++;
        } else {
            diff = 1;
        }
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


@d C integer x Input
@{
    int *x,
    const R_xlen_t N,
    const int P,
@}

@d C integer y Input
@{
    int *y,
    const int Q,
@}

@d TableSums Inputs
@{
    @<C integer x Input@>
    @<C integer y Input@>
@}

@d TableSums Answer
@{
    double *PQ_ans
@}

@d C\_TableSums\_dweights\_dsubset
@{
void C_TableSums_dweights_dsubset
(
    @<TableSums Inputs@>
    double *weights,
    const int HAS_WEIGHTS,
    @<C real subset Input@>,
    @<TableSums Answer@>
) {

    double *s, *w; 

  @<TableSums Body@>
}
@}

@d C\_TableSums\_iweights\_dsubset
@{
void C_TableSums_iweights_dsubset
(
    @<TableSums Inputs@>
    int *weights,
    const int HAS_WEIGHTS,
    @<C real subset Input@>,
    @<TableSums Answer@>
) {

    double *s;
    int *w; 

  @<TableSums Body@>
}
@}

@d C\_TableSums\_iweights\_isubset
@{
void C_TableSums_iweights_isubset
(
    @<TableSums Inputs@>
    int *weights,
    const int HAS_WEIGHTS,
    @<C integer subset Input@>,
    @<TableSums Answer@>
) {

    int *s, *w;

  @<TableSums Body@>
}
@}

@d C\_TableSums\_dweights\_isubset
@{
void C_TableSums_dweights_isubset
(
    @<TableSums Inputs@>
    double *weights,
    const int HAS_WEIGHTS,
    @<C integer subset Input@>,
    @<TableSums Answer@>
) {

    int *s; 
    double *w;

  @<TableSums Body@>
}
@}

@d RC\_TableSums
@{

void RC_TableSums
(
    @<TableSums Inputs@>
    SEXP weights,
    SEXP subset,
    R_xlen_t offset,
    R_xlen_t Nsubset,
    @<TableSums Answer@>
) {

    if (TYPEOF(weights) == INTSXP) {
        if (TYPEOF(subset) == INTSXP) {
            C_TableSums_iweights_isubset(x, N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_TableSums_iweights_dsubset(x, N, P, y, Q, 
                                        INTEGER(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    } else {
        if (TYPEOF(subset) == INTSXP) {
            C_TableSums_dweights_isubset(x, N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, INTEGER(subset), 
                                        offset, Nsubset, PQ_ans);
        } else {
            C_TableSums_dweights_dsubset(x, N, P, y, Q, 
                                        REAL(weights), XLENGTH(weights) > 0, REAL(subset), 
                                        offset, Nsubset, PQ_ans);
        }
    }
}

@}

@d R\_TableSums
@{
SEXP R_TableSums
(
    SEXP x,
    SEXP P,
    SEXP y,
    SEXP Q,
    SEXP weights,
    SEXP subset
) {

    SEXP ans;
    R_xlen_t N, Nsubset;
    double *center;

    N = XLENGTH(x);
    Nsubset = XLENGTH(subset);
    
    PROTECT(ans = allocVector(REALSXP, INTEGER(P)[0] * INTEGER(Q)[0]));
    RC_TableSums(INTEGER(x), N, INTEGER(P)[0], INTEGER(y), INTEGER(Q)[0], 
                 weights, subset, 0, Nsubset, REAL(ans));
    UNPROTECT(1);
    return(ans);
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
