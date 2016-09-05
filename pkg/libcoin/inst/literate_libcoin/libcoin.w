
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
of a specific linear form. Let $A \subseteq \{1, \dots, N\}$ denote some subset of the
observation numbers and consider the linear statistic
\begin{eqnarray} \label{linstat}
\T(A) = \vec\left(\sum_{i \in A} w_i g(\X_i) h(\Y_i, \{\Y_i \mid i \in A\})^\top\right)
\in \R^{pq}.  
\end{eqnarray}
Here, $g: \mathcal{X} \rightarrow \R^{p}$ is a transformation of
$\X$ known as the \emph{regression function} and $h: \mathcal{Y} \times 
\mathcal{Y}^n \rightarrow \R^q$ is a transformation of $\Y$ known as the
\emph{influence function}, where the latter may depend on $\Y_i$ for $i \in A$
in a permutation symmetric way.  We will give specific examples on how to choose
$g$ and $h$ later on.

With $\x_i = g(\X_i)$ and $\y_i = h(\Y_i, \{\Y_i, i \in A\})$ we write
\begin{eqnarray} \label{linstat}
\T(A) = \vec\left(\sum_{i \in A} w_i \x_i \y_i^\top\right)
\in \R^{pq}.  
\end{eqnarray}
The \pkg{libcoin} package doesn't handle neither $g$ nor $h$, this is the job
of \pkg{coin} and we therefore continue with $\x_i$ and $\y_i$.

The distribution of $\T$  depends on the joint
distribution of $\Y$ and $\X$, which is unknown under almost all practical  
circumstances. At least under the null hypothesis one can dispose of this
dependency by fixing $\X_i, i \in A$ and conditioning on all possible
permutations $S(A)$ of the responses $\Y_i, i \in A$.
This principle leads to test procedures known
as \textit{permutation tests}.
The conditional expectation $\mu(A) \in \R^{pq}$ and covariance
$\Sigma(A) \in \R^{pq \times pq}$
of $\T$ under $H_0$ given
all permutations $\sigma \in S(A)$ of the responses are derived by
\cite{strasserweber1999}:
\begin{eqnarray}
\mu(A) & = & \E(\T(A) \mid S(A)) = \vec \left( \left( \sum_{i \in A} w_i \x_i \right) \E(h \mid S(A))^\top
\right), \nonumber \\
\Sigma(A) & = & \V(\T(A) \mid S(A)) \nonumber \\
& = &
    \frac{\ws}{\ws(A) - 1}  \V(h \mid S(A)) \otimes
        \left(\sum_{i \in A} w_i  \x_i \otimes w_i \x_i^\top \right)
\label{expectcovar}
\\
& - & \frac{1}{\ws(A) - 1}  \V(h \mid S(A))  \otimes \left(
        \sum_{i \in A} w_i \x_i \right)
\otimes \left( \sum_{i \in A} w_i \x_i\right)^\top
\nonumber
\end{eqnarray}
where $\ws(A) = \sum_{i \in A} w_i$ denotes the sum of the case weights,
and $\otimes$ is the Kronecker product. The conditional expectation of the
influence function is
\begin{eqnarray*}
\E(h \mid S(A)) = \ws(A)^{-1} \sum_{i \in A} w_i \y_i \in
\R^q
\end{eqnarray*}
with corresponding $q \times q$ covariance matrix
\begin{eqnarray*}
\V(h \mid S(A)) = \ws(A)^{-1} \sum_{i \in A} w_i \left(\y_i - \E(h \mid S(A)) \right) \left(\y_i  - \E(h \mid S(A))\right)^\top.
\end{eqnarray*}

With $A_b = \{i \mid b_i = b\}$ we get $\T = \sum_{b = 1}^B T(A_b)$, 
$\mu = \sum_{b = 1}^B \mu(A_b)$ and $\Sigma = \sum_{b = 1}^B \Sigma(A_b)$.

Having the conditional expectation and covariance at hand we are able to
standardize a linear statistic $\T \in \R^{pq}$ of the form
(\ref{linstat}). Univariate test statistics~$c$ mapping an observed linear
statistic $\mathbf{t} \in
\R^{pq}$ into the real line can be of arbitrary form.  An obvious choice is
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
the \pkg{coin} package is
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

\section{Includes}

@o libcoin_internal.h -cc
@{
@<Include R API@>
@<Define macros@>
@<Define input constants@>
@<Define output constants@>
@}

We include the relevant header files from R.

@d Include R API
@{
/* R API */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Lapack.h> /* for dgesdd */
@}

\verb|S[i, j, n]| computes the index of element $(i, j)$ of
some $n \times n$ symmetric matrix stored in lower packaged
form, allowing for $i < j$. \verb|LE| and \verb|GE| 
implement $\le$ and $\ge$ with given tolerance \verb|tol|.

@d Define macros
@{
/* S[i, j] for n x n symmetric matrix in lower packed 
   storage allowing for i < j */
#define S(i, j, n) ((i) >= (j) ? (n) * (j) + (i) - (j) * ((j) + 1) / 2 : (n) * (i) + (j) - (i) * ((i) + 1) / 2)
#define LE(x, y, tol)  ((x) < (y)) || (fabs((x) - (y)) < (tol))
#define GE(x, y, tol)  ((x) > (y)) || (fabs((x) - (y)) < (tol))
@| S LE GE
@}

@d Define input constants
@{
#define ALTERNATIVE_twosided		1
#define ALTERNATIVE_less                2
#define ALTERNATIVE_greater             3
#define TESTSTAT_maximum                1
#define TESTSTAT_quadratic              2
@| ALTERNATIVE_twosided ALTERNATIVE_less ALTERNATIVE_greater          
TESTSTAT_maximum TESTSTAT_quadratic 
@}

@d Define output constants
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
@| LinearStatistic_SLOT Expectation_SLOT
Covariance_SLOT Variance_SLOT MPinv_SLOT ExpectationX_SLOT 
varonly_SLOT dim_SLOT ExpectationInfluence_SLOT CovarianceInfluence_SLOT 
VarianceInfluence_SLOT Xfactor_SLOT Work_SLOT tol_SLOT 
PermutedLinearStatistic_SLOT TableBlock_SLOT Sumweights_SLOT
Table_SLOT
@}

\section{Sums}

@s

The \proglang{C} file \verb|Sums.c| defines the \proglang{C}
functions and a corresponding \proglang{R} interface (via \verb|.C()|)

@o Sums.c -cc
@{
#include "libcoin_internal.h"
@<Sums of weights@>
@<Row sums@>
@<Column sums@>
@}

The corresponding header file contains definitions of
functions to be used outside \verb|Sums.c|

@o Sums.h -cc
@{
@<Function prototypes@>
@}

The \proglang{R} interfaces are used to implement 
regression tests to be called from within \proglang{R}

<<regression-test>>=
### replace with library("libcoin")
dyn.load("Sums.so")
set.seed(29)
N <- 20L
P <- 3L
Q <- 2L
x <- matrix(runif(N * P), nrow = N)
y <- matrix(runif(N * Q), nrow = N)
weights <- sample(0:3, size = N, replace = TRUE)
subset <- sample(1:N, size = N, replace = TRUE)
centerx <- as.double(runif(P))
P_ans <- as.double(numeric(P))
@@

\subsection{Sums of weights}

@d Function prototypes
@{
/* sum(weights) */
extern int C_sum_weights(const *double, const int);
extern int C_sum_weights_subset(const *double, const *int, const int);
@}

For an $N$ vector of integer \verb|weights| $w = (w_0, \dots, w_{N - 1})$

@d Input weights
@{
const int *weights
@| weights
@}

we compute the sum of certain elements of $w$ by the scrap

@d Addition of @'N@' summands @'summands@'
@{
  long int ans = 0;
  for (int i = 0; i < @1; i++) ans += @2;
@}

Because the sum can be larger than \verb|INT_MAX|, the sum is
stored as a \verb|long int| first and a cast to \verb|int| performed 
if this is possible, otherwise an error is reported.

@d Check Integer Overlow and Return
@{
/* integer overflow can only happen here */
if (ans > INT_MAX) 
  error("sum of weights is larger than INT_MAX");
int ret = (int) ans;
return(ret);
@}

We compute the total sum $\sum_{i = 0}^{N - 1} w_i$ 

@d Sums of weights
@{
/* sum(weights) */
int C_sum_weights
(
  @<Input weights@>,
  int N
) {
  @<Addition of @'N@' summands @'weights[i]@'@>
  @<Check Integer Overlow and Return@>
}
void RC_sum_weights
(
  @<Input weights@>,
  int *N,
  int *ans
) {
  ans[0] = C_sum_weights(weights, N[0]);
}
@| C_sum_weights
@}

and regression test

<<>>=
ans <- .C("RC_sum_weights", weights = as.integer(weights),
                            N = as.integer(N),
                            ans = integer(1))$ans
stopifnot(all.equal(ans, sum(weights)))
@@

and the sum over a subset $\sum_{i \in S} w_i$ 
where the subset $S$ is represented by an integer vector 
\verb|subsetx| of length \verb|Nsubset| 

@d Input subsetx
@{
const int *subsetx,
@| subsetx
@}

@d Sums of weights
@{
/* sum(weights[subsetx]) */
int C_sum_weights_subset
(
  @<Input weights@>,
  const int N,
  @<Input subsetx@>
  const int Nsubset
) {
  @<Addition of @'Nsubset@' summands @'weights[subsetx[i]]@'@>
  @<Check Integer Overlow and Return@>
}
void RC_sum_weights_subset
(
  @<Input weights@>,
  const int *N,
  @<Input subsetx@>
  const int *Nsubset,
  int *ans
) {
  ans[0] = C_sum_weights_subset(weights, N[0], subsetx, Nsubset[0]);
}
@| C_sum_weights_subset
@}

and regression test

<<>>=
ans <- .C("RC_sum_weights_subset", 
          weights = as.integer(weights),
          N = as.integer(N),
          subset = as.integer(subset - 1L),
          Nsubset = as.integer(length(subset)),
          ans = integer(1))$ans
stopifnot(all.equal(ans, sum(weights[subset])))
@@

\subsection{Row sums}

For a double matrix $X \sim (N, P)$ 
@d Input xNP 
@{
const double *x,
const int N,
const int P,
@}

@d R Input xNP 
@{
const double *x,
const int *N,
const int *P,
@}

@d Function prototypes
@{
external void C_rowSums_(@<Input xNP@>, double *N_ans);
external void C_rowSums_i(const int *x, const int N,  
                          const int P, int *N_ans);
@}

The $N$ vector \verb|N_ans| of row sums 
$\sum_{p = 0}^{P - 1} X_{ip}, i = 0, \dots, N - 1$ for a double
matrix $X \sim (N, P)$ is computed by  

@d Row sum scrap
@{
for (int i = 0; i < N; i++) {
    N_ans[i] = 0.0;
    for (int p = 0; p < P; p++)
        N_ans[i] += x[p * N + i];
}
@}

@d Row sums
@{
/* rowSums(x) */
void C_rowSums_
(
  @<Input xNP@>
  double *N_ans
) {
  @<Row sum scrap@>
}
void RC_rowSums_
(
  @<R Input xNP@>
  double *N_ans
) {
  C_rowSums_(x, N[0], P[0], N_ans);
}
@| C_rowSums_
@}

<<>>=
N_ans <- .C("RC_rowSums_", x = as.double(x), N = as.integer(N),
                           P = as.integer(P), N_ans = double(N))$N_ans
stopifnot(all.equal(N_ans, rowSums(x)))
@@

A version for integer $X$ is

@d Row sums
@{
/* rowSums(x), integer version */
void C_rowSums_i
(
  const int *x,
  const int N,
  const int P,
  int *N_ans
) {
  @<Row sum scrap@>
}
@| C_rowSums_i
@}

\subsection{Column sums}


@d Input xNP weights 
@{
@<Input xNP@>
@<Input weights@>,
@}

@d Input xNP subset
@{
@<Input xNP@>
@<Input subsetx@>
const int Nsubset,
@}

@d Input xNP weights subset 
@{
@<Input xNP weights@>
@<Input subsetx@>
const int Nsubset,
@}


@d R Input xNP weights 
@{
@<R Input xNP@>
@<Input weights@>,
@}

@d R Input xNP subset
@{
@<R Input xNP@>
@<Input subsetx@>
const int *Nsubset,
@}

@d R Input xNP weights subset 
@{
@<R Input xNP weights@>
@<Input subsetx@>
const int *Nsubset,
@}

\subsubsection{Simple columns sums}

@d Function prototypes
@{
external void C_colSums_(@<Input xNP@>, double *P_ans);
external void C_colSums_i(const int *x, const int N,  
                          const int P, int *P_ans);
external void C_colSums_weights(@<Input xNP weights@>, double *P_ans);
external void C_colSums_subset(@<Input xNP subset@>, double *P_ans);
external void C_colSums_weights_subset(@<Input xNP weights subset@>, double *P_ans);
@}

@d Column Sum of @'sum@' with @'index@'
@{
for (int p = 0; p < P; p++) {
    P_ans[p] = 0; /* initialise with 0 */
    int pN = p * N;
    for (int i = 0; i < @2; i++)
        P_ans[p] += @1;
}
@}

the $P$ vector of column sums $\sum_{i = 0}^{N - 1} X_{ip}, p = 0, \dots, P - 1$
is computed by \verb|C_colSums_| and written to a $P$ vector pointed to by
\verb|P_ans|

@d Column sums
@{
/* colSums(x) */
void C_colSums_
(
  @<Input xNP@>
  double *P_ans  /* Return P_ans */
) {
  @<Column Sum of @'x[pN + i]@' with @'N@'@>
}
void RC_colSums_ 
(
  @<R Input xNP@>
  double *P_ans  /* Return P_ans */
) {
  C_colSums_(x, N[0], P[0], P_ans);
}
@| C_colSums_
@}

<<>>=
P_ans <- .C("RC_colSums_", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(x)))
@@

A version for an integer matrix $X$ 

@d Column sums
@{
/* colSums(x), integer version */
void C_colSums_i
(
  const int *x, const int N,
  const int P, double *P_ans
) {
  @<Column Sum of @'x[pN + i]@' with @'N@'@>
}
@| C_colSums_i
@}

@d Column sums
@{
/* colSums(x * weights) */
void C_colSums_weights
(
  @<Input xNP weights@>
  double *P_ans
){
  @<Column Sum of @'weights[i] * x[pN + i]@' with @'N@'@>
}
void RC_colSums_weights
(
  @<R Input xNP weights@>
  double *P_ans
){
  C_colSums_weights(x, N[0], P[0], weights, P_ans);
}
@| C_colSums_weights
@}

<<>>=
P_ans <- .C("RC_colSums_weights", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           weights = as.integer(weights),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(x * weights)))
@@

@d Column sums
@{
/* colSums(x[subsetx,]) */
void C_colSums_subset
(
  @<Input xNP subset@>
  double *P_ans
){
  @<Column Sum of @'x[pN + subsetx[i]]@' with @'Nsubset@'@>
}
void RC_colSums_subset
(
  @<R Input xNP subset@>
  double *P_ans
){
  C_colSums_subset(x, N[0], P[0], subsetx, Nsubset[0], P_ans);
}
@| C_colSums_subset
@}

<<>>=
P_ans <- .C("RC_colSums_subset", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           subset = as.integer(subset - 1L),
                           Nsubset = as.integer(length(subset)),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(x[subset,])))
@@


@d Column sums
@{
/* colSums(x[subsetx,] * weights[subsetx]) */
void C_colSums_weights_subset
(
  @<Input xNP weights subset@>
  double *P_ans
){
  @<Column Sum of @'weights[subsetx[i]] * x[pN + subsetx[i]]@' with @'Nsubset@'@>
}
void RC_colSums_weights_subset
(
  @<R Input xNP weights subset@>
  double *P_ans
){
  C_colSums_weights_subset(x, N[0], P[0], weights, subsetx, Nsubset[0], P_ans);
}
@| C_colSums_weights_subset
@}

<<>>=
P_ans <- .C("RC_colSums_weights_subset", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           weights = as.integer(weights),
                           subset = as.integer(subset - 1L),
                           Nsubset = as.integer(length(subset)),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(x[subset,] * weights[subset])))
@@


\subsubsection{Squared columns sums}

@d Function prototypes
@{
external void C_colSums2_(@<Input xNP@>, double *P_ans);
external void C_colSums2_weights(@<Input xNP weights@>, double *P_ans);
external void C_colSums2_subset(@<Input xNP subset@>, double *P_ans);
external void C_colSums2_weights_subset(@<Input xNP weights subset@>, double *P_ans);
@}

the $P$ vector of column sums $\sum_{i = 0}^{N - 1} X^2_{ip}, p = 0, \dots, P - 1$
is computed by \verb|C_colSums_| and written to a $P$ vector pointed to by
\verb|P_ans|

@d Column sums
@{
/* colSums(x^2) */
void C_colSums2_
(
  @<Input xNP@>
  double *P_ans  /* Return P_ans */
) {
  @<Column Sum of @'pow(x[pN + i], 2)@' with @'N@'@>
}
void RC_colSums2_ 
(
  @<R Input xNP@>
  double *P_ans  /* Return P_ans */
) {
  C_colSums2_(x, N[0], P[0], P_ans);
}
@| C_colSums2_
@}

<<>>=
P_ans <- .C("RC_colSums2_", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(x^2)))
@@

@d Column sums
@{
/* colSums(x^2 * weights) */
void C_colSums2_weights
(
  @<Input xNP weights@>
  double *P_ans
){
  @<Column Sum of @'weights[i] * pow(x[pN + i], 2)@' with @'N@'@>
}
void RC_colSums2_weights
(
  @<R Input xNP weights@>
  double *P_ans
){
  C_colSums2_weights(x, N[0], P[0], weights, P_ans);
}
@| C_colSums2_weights
@}

<<>>=
P_ans <- .C("RC_colSums2_weights", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           weights = as.integer(weights),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(x^2 * weights)))
@@

@d Column sums
@{
/* colSums(x[subsetx,]^2) */
void C_colSums2_subset
(
  @<Input xNP subset@>
  double *P_ans
){
  @<Column Sum of @'pow(x[pN + subsetx[i]], 2)@' with @'Nsubset@'@>
}
void RC_colSums2_subset
(
  @<R Input xNP subset@>
  double *P_ans
){
  C_colSums2_subset(x, N[0], P[0], subsetx, Nsubset[0], P_ans);
}
@| C_colSums2_subset
@}

<<>>=
P_ans <- .C("RC_colSums2_subset", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           subset = as.integer(subset - 1L),
                           Nsubset = as.integer(length(subset)),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(x[subset,]^2)))
@@


@d Column sums
@{
/* colSums(x[subsetx,]^2 * weights[subsetx]) */
void C_colSums2_weights_subset
(
  @<Input xNP weights subset@>
  double *P_ans
){
  @<Column Sum of @'weights[subsetx[i]] * pow(x[pN + subsetx[i]], 2)@' with @'Nsubset@'@>
}
void RC_colSums2_weights_subset
(
  @<R Input xNP weights subset@>
  double *P_ans
){
  C_colSums2_weights_subset(x, N[0], P[0], weights, subsetx, Nsubset[0], P_ans);
}
@| C_colSums2_weights_subset
@}

<<>>=
P_ans <- .C("RC_colSums2_weights_subset", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           weights = as.integer(weights),
                           subset = as.integer(subset - 1L),
                           Nsubset = as.integer(length(subset)),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(x[subset,]^2 * weights[subset])))
@@


\subsubsection{Squared centered columns sums}

@d Function prototypes
@{
external void C_colSums2_center(@<Input xNP@>, const double *centerx, double *P_ans);
external void C_colSums2_center_weights(@<Input xNP weights@>, const double *centerx, double *P_ans);
external void C_colSums2_center_subset(@<Input xNP subset@>, const double *centerx, double *P_ans);
external void C_colSums2_center_weights_subset(@<Input xNP weights subset@>, const double *centerx, double *P_ans);
@}

the $P$ vector of column sums $\sum_{i = 0}^{N - 1} X^2_{ip}, p = 0, \dots, P - 1$
is computed by \verb|C_colSums_| and written to a $P$ vector pointed to by
\verb|P_ans|

@d Column sums
@{
/* colSums(x^2) */
void C_colSums2_center
(
  @<Input xNP@>
  const double *centerx,
  double *P_ans  /* Return P_ans */
) {
  @<Column Sum of @'pow(x[pN + i] - centerx[p], 2)@' with @'N@'@>
}
void RC_colSums2_center 
(
  @<R Input xNP@>
  const double *centerx,
  double *P_ans  /* Return P_ans */
) {
  C_colSums2_center(x, N[0], P[0], centerx, P_ans);
}
@| C_colSums2_center
@}

<<>>=
P_ans <- .C("RC_colSums2_center", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           centerx = as.double(centerx),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(t((t(x) - centerx))^2)))
@@

@d Column sums
@{
/* colSums(x^2 * weights) */
void C_colSums2_center_weights
(
  @<Input xNP weights@>
  const double *centerx,
  double *P_ans
){
  @<Column Sum of @'weights[i] * pow(x[pN + i] - centerx[p], 2)@' with @'N@'@>
}
void RC_colSums2_center_weights
(
  @<R Input xNP weights@>
  const double *centerx,
  double *P_ans
){
  C_colSums2_center_weights(x, N[0], P[0], weights, centerx, P_ans);
}
@| C_colSums2_center_weights
@}

<<>>=
P_ans <- .C("RC_colSums2_center_weights", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           weights = as.integer(weights),
                           centerx = as.double(centerx),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(t((t(x) - centerx))^2 * weights)))
@@

@d Column sums
@{
/* colSums(x[subsetx,]^2) */
void C_colSums2_center_subset
(
  @<Input xNP subset@>
  const double *centerx,
  double *P_ans
){
  @<Column Sum of @'pow(x[pN + subsetx[i]] - centerx[p], 2)@' with @'Nsubset@'@>
}
void RC_colSums2_center_subset
(
  @<R Input xNP subset@>
  const double *centerx,
  double *P_ans
){
  C_colSums2_center_subset(x, N[0], P[0], subsetx, Nsubset[0], centerx, P_ans);
}
@| C_colSums2_center_subset
@}

<<>>=
P_ans <- .C("RC_colSums2_center_subset", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           subset = as.integer(subset - 1L),
                           Nsubset = as.integer(length(subset)),
                           centerx = as.double(centerx),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(t((t(x[subset,]) - centerx))^2)))
@@


@d Column sums
@{
/* colSums(x[subsetx,]^2 * weights[subsetx]) */
void C_colSums2_center_weights_subset
(
  @<Input xNP weights subset@>
  const double *centerx,
  double *P_ans
){
  @<Column Sum of @'weights[subsetx[i]] * pow(x[pN + subsetx[i]] - centerx[p], 2)@' with @'Nsubset@'@>
}
void RC_colSums2_center_weights_subset
(
  @<R Input xNP weights subset@>
  const double *centerx,
  double *P_ans
){
  C_colSums2_center_weights_subset(x, N[0], P[0], weights, subsetx, Nsubset[0], centerx, P_ans);
}
@| C_colSums2_center_weights_subset
@}

<<>>=
P_ans <- .C("RC_colSums2_center_weights_subset", x = as.double(x),
                           N = as.integer(N),
                           P = as.integer(P),
                           weights = as.integer(weights),
                           subset = as.integer(subset - 1L),
                           Nsubset = as.integer(length(subset)),
                           centerx = as.double(centerx),
                           P_ans = as.double(P_ans))$P_ans
stopifnot(all.equal(P_ans, colSums(t((t(x[subset,]) - centerx))^2 * weights[subset])))
@@


@S



\section{R Code}

@s

@o ctabs.R -cp
@{
@<Cross Tabulation@>
@}

@d Cross Tabulation
@{
ctabs <- function(ix, iy = integer(0), weights = integer(0),
                  subset = integer(0), block = integer(0))
{

    if (is.null(attr(ix, "levels")))
            attr(ix, "levels") <- 1:max(ix)

    if (length(iy) > 0) {
        if (is.null(attr(iy, "levels")))
            attr(iy, "levels") <- 1:max(iy)
    }

    if (length(subset) > 0) subset <- subset - 1L

    ret <- .Call("R_tables", ix, iy, weights, subset, block,
                 PACKAGE = "libcoin")
    if (length(block) > 0) {
        if (length(iy) == 0)
            ret <- ret[,,-1, drop = FALSE]
        else
            ret <- ret[,,-dim(ret)[3], drop = FALSE]
    } else {
        ret <- ret[,,,drop = TRUE]
    }
    ret
}
@| ctabs
@}

<<>>=
source("ctabs.R")
library("libcoin")

x <- 1:6
y <- rep(1:3, length.out = 6)
w <- rep(2L, length(x))
b <- gl(2, 3)
s <- 1:4

ctabs(x)
ctabs(x, weights = w)
ctabs(x, subset = s)
ctabs(x, weights = w, subset = s)
ctabs(x, block = b)
ctabs(x, weights = w, block = b)
ctabs(x, subset = s, block = b)
ctabs(x, weights = w, subset = s, block = b)

ctabs(x, y)
ctabs(x, y, weights = w)
ctabs(x, y, subset = s)
ctabs(x, y, weights = w, subset = s)
ctabs(x, y, block = b)
ctabs(x, y, weights = w, block = b)
ctabs(x, y, subset = s, block = b)
ctabs(x, y, weights = w, subset = s, block = b)

w <- round(1:length(x) / length(x), 2)
ctabs(x, weights = w)
ctabs(x, weights = w, subset = s)
ctabs(x, weights = w, block = b)
ctabs(x, weights = w, subset = s, block = b)

ctabs(x, y, weights = w)
ctabs(x, y, weights = w, subset = s)
ctabs(x, y, weights = w, block = b)
ctabs(x, y, weights = w, subset = s, block = b)
@@


\section{Example}

<<ex>>=
summary(1)
@@

@S

\end{document}
