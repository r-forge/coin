
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


\author{Torsten Hothorn \\ Universit\"at Z\"urich}

\title{The \pkg{libcoin} Package}

\begin{document}

\pagenumbering{roman}
\maketitle
\tableofcontents

\chapter{Introduction}
\pagenumbering{arabic}

\verb|code|

\chapter{The libcoin package}

\section{C Code}

\subsection{Includes}

@o libcoin_internal.h -cc
@{
@<Include R API@>
@<Define Macros@>
@<Define Input Constants@>
@<Define Output Constants@>
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

@d Define Macros
@{
/* S[i, j] for n x n symmetric matrix in lower packed 
   storage allowing for i < j */
#define S(i, j, n) ((i) >= (j) ? (n) * (j) + (i) - (j) * ((j) + 1) / 2 : (n) * (i) + (j) - (i) * ((i) + 1) / 2)
#define LE(x, y, tol)  ((x) < (y)) || (fabs((x) - (y)) < (tol))
#define GE(x, y, tol)  ((x) > (y)) || (fabs((x) - (y)) < (tol))
@| S LE GE
@}

@d Define Input Constants
@{
#define ALTERNATIVE_twosided		1
#define ALTERNATIVE_less                2
#define ALTERNATIVE_greater             3
#define TESTSTAT_maximum                1
#define TESTSTAT_quadratic              2
@| ALTERNATIVE_twosided ALTERNATIVE_less ALTERNATIVE_greater          
TESTSTAT_maximum TESTSTAT_quadratic 
@}

@d Define Output Constants
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

\subsection{Sums}

@s

@o Sums.c -cc
@{
#include "libcoin_internal.h"
@<Sums of Weights@>
@<Row and Column Sums@>
@}

\subsubsection{Sums of Weights}

For an $N$ vector of integer \verb|weights| $(w_0, \dots, w_{N - 1})$

@d Input weights, N
@{
@<Input weights@>,
@<Input N@>
@}


@d Input weights
@{
const int *weights
@| weights
@}

@d Input N
@{
const int N
@| N
@}


we compute the total sum $\sum_{i = 0}^{N - 1} w_i$ (\verb|C_sum_weights|) 
or the sum over a subset $\sum_{i \in S} w_i$ (\verb|C_sum_weights_subset|), 
where the subset $S$ is represented by an integer vector \verb|subsetx| of length 
\verb|Nsubset| 

@d Input subsetx, Nsubset
@{
const int *subsetx,
const int Nsubset
@| subsetx Nsubset
@}

@d Sum over @'subset@' of length @'N@'
@{
  long int ans = 0;
  for (int i = 0; i < @2; i++) ans += weights[@1];
@}

@d Sums of Weights
@{
/* sum(weights) */
int C_sum_weights
(
  @<Input weights, N@>
) {
  @<Sum over @'i@' of length @'N@'@>
  @<Check Integer Overlow and Return@>
}

/* sum(weights[subsetx]) */
int C_sum_weights_subset
(
  @<Input weights, N@>,
  @<Input subsetx, Nsubset@>
) {
  @<Sum over @'subsetx[i]@' of length @'Nsubset@'@>
  @<Check Integer Overlow and Return@>
}
@| C_sum_weights C_sum_weights_subset
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

\subsubsection{Row and Column Sums}

@d Row and Column Sums
@{
@<Double Column Sum@>
@<Integer Column Sum@>
@<Double Row Sum@>
@<Integer Row Sum@>
@<Double Column Sum with Weights@>
@<Double Column Sum with Subset@>
@<Double Column Sum with Weights and Subset@>
@<Double Squared Column Sum@>
@<Double Squared Column Sum with Weights@>
@<Double Squared Column Sum with Subset@>
@<Double Squared Column Sum with Weights and Subset@>
@| C_colSums_ C_rowSums_i C_colSums_ C_rowSums_i C_colSums_weights
@}

For a double matrix $X \sim (N, P)$ 
@d Input double x, N, P
@{
@<Input double x@>,
@<Input N@>,
@<Input P@>
@}

@d Input double x
@{
const double *x
@}

@d Input P
@{
const int P
@}


the $P$ vector of column sums $\sum_{i = 0}^{N - 1} X_{ip}, p = 0, \dots, P - 1$
is computed by \verb|C_colSums_| and written to a $P$ vector pointed to by
\verb|P_ans|

@d Double Column Sum
@{
/* colSums(x) */
void C_colSums_
(
  @<Input double x, N, P@>,
  double *P_ans  /* Return P_ans */
) {
  @<Column Sum of @'x[pN + i]@' with @'N@'@>
}
@| C_colSums_
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

A version for an integer matrix $X$ 
@d Input int x, N, P
@{
const int *x,
@<Input N@>,
@<Input P@>
@}

@d Integer Column Sum
@{
/* colSums(x), integer version */
void C_colSums_i
(
  @<Input int x, N, P@>,
  int *P_ans
) {
  @<Column Sum of @'x[pN + i]@' with @'N@'@>
}
@| C_colSums_i
@}

The $N$ vector \verb|N_ans| of row sums 
$\sum_{p = 0}^{P - 1} X_{ip}, i = 0, \dots, N - 1$ for a double
matrix $X \sim (N, P)$ is computed by  \verb|C_rowSums_|

@d Double Row Sum
@{
/* rowSums(x) */
void C_rowSums_
(
  @<Input double x, N, P@>,
  double *N_ans
) {
  @<Row Sum@>
}
@| C_rowSums_
@}

@d Row Sum
@{
for (int i = 0; i < N; i++) {
    N_ans[i] = 0.0;
    for (int p = 0; p < P; p++)
        N_ans[i] += x[p * N + i];
}
@}

and a version for integer $X$ is

@d Integer Row Sum
@{
/* rowSums(x), integer version */
void C_rowSums_i
(
  @<Input int x, N, P@>,
  int *N_ans
) {
  @<Row Sum@>
}
@| C_rowSums_i
@}

@d Double Column Sum with Weights
@{
/* colSums(x * weights) */
void C_colSums_weights
(
  @<Input double x, N, P@>,
  @<Input weights@>,
  double *P_ans
){
  @<Column Sum of @'weights[i] * x[pN + i]@' with @'N@'@>
}
@| C_colSums_weights
@}

@d Double Column Sum with Subset
@{
/* colSums(x[subsetx,]) */
void C_colSums_subset
(
  @<Input double x, N, P@>,
  @<Input subsetx, Nsubset@>,
  double *P_ans
){
  @<Column Sum of @'x[pN + subsetx[i]]@' with @'Nsubset@'@>
}
@| C_colSums_subset
@}

@d Double Column Sum with Weights and Subset
@{
/* colSums(x[subsetx,]) */
void C_colSums_weights_subset
(
  @<Input double x, N, P@>,
  @<Input weights@>,
  @<Input subsetx, Nsubset@>,
  double *P_ans
){
  @<Column Sum of @'weights[subsetx[i]] * x[pN + subsetx[i]]@' with @'Nsubset@'@>
}
@| C_colSums_weights_subset
@}

@d Double Squared Column Sum
@{
/* colSums(x^2) */
void C_colSums2_
(
  @<Input double x, N, P@>,
  @<Input subsetx, Nsubset@>,
  double *P_ans
){
  @<Column Sum of @'pow(x[pN + i], 2)@' with @'N@'@>
}
@| C_colSums2_
@}

@d Double Squared Column Sum with Weights
@{
/* colSums(x^2 * weights) */
void C_colSums2_weights
(
  @<Input double x, N, P@>,
  @<Input weights@>,
  double *P_ans
){
  @<Column Sum of @'weights[i] * pow(x[pN + i], 2)@' with @'N@'@>
}
@| C_colSums2_weigths
@}

@d Double Squared Column Sum with Subset
@{
/* colSums(x^2 * weights) */
void C_colSums2_subset
(
  @<Input double x, N, P@>,
  @<Input subsetx, Nsubset@>,
  double *P_ans
){
  @<Column Sum of @'pow(x[pN + subsetx[i]], 2)@' with @'Nsubset@'@>
}
@| C_colSums2_subset
@}

@d Double Squared Column Sum with Weights and Subset
@{
/* colSums(x^2 * weights) */
void C_colSums2_weights_subset
(
  @<Input double x, N, P@>,
  @<Input weights@>,
  @<Input subsetx, Nsubset@>,
  double *P_ans
){
  @<Column Sum of @'weights[subsetx[i]] * pow(x[pN + subsetx[i]], 2)@' with @'Nsubset@'@>
}
@| C_colSums2_weights_subset
@}


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

@o regtest_ctabs.R -cp
@{
@<Cross Tabulation Tests@>
@}

@d Cross Tabulation Tests
@{
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
@}


@S

\end{document}
