\name{LocationTests}
\alias{perm_test}
\alias{wilcox_test}
\alias{normal_test}
\alias{median_test}
\alias{kruskal_test}
\title{ Independent Two- and K-Sample Location Tests }
\description{
    Testing the equality of the distributions of a numeric response in
    two or more independent groups against shift alternatives.
}
\usage{

\method{perm_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{perm_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), ...)

\method{wilcox_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{wilcox_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),
    conf.int = FALSE, conf.level = 0.95, ...)

\method{normal_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{normal_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...)

\method{median_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{median_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),
    conf.int = FALSE, conf.level = 0.95, ...)

\method{kruskal_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{kruskal_test}{IndependenceProblem}(x, 
    distribution = c("asympt", "approx"), ...)

}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    is a numeric variable giving the data values and \code{x} a factor
    with two or more levels giving the corresponding groups. \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
      integer valued weights for the observations.}
  \item{x}{an object of class \code{IndependenceProblem}.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or    
    \code{"less"}.  You can specify just the initial letter.}   
  \item{distribution}{the null distribution of the test statistic
    can be computed (\code{exact}) or can be approximated by its
    asymptotic distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{ties.method}{two methods are available to adjust scores for ties.
      Either the score generating function is applied to \code{mid-ranks}
      or scores, based on random ranks, are averaged \code{average-scores}.}
  \item{conf.int}{a logical indicating whether a confidence interval
    for the difference in location should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of \code{y} in
  the groups given by \code{x} is tested. In particular, the methods
  documented here are designed to detect shift alternatives. 

  The test procedures apply a rank transformation to the response values 
  (\code{y}), except of \code{perm_test} which computes a test statistic
  using the raw response values.

  The asymptotic null distribution is computed by default for all
  procedures. Exact p-values may be computed for the two-sample problems and
  can be approximated via Monte-Carlo for all procedures. Exact p-values
  are computed either by the shift algorithm (Streitberg & Roehmel, 1986,
  1987) or by the split-up algorithm (van de Wiel, 2001).

  Confidence intervals for the difference in location are available for the 
  rank-based procedures and are computed according to Bauer (1972). 

  In case \code{x} is an ordered factor, \code{kruskal_test} computes the
  linear-by-linear association test for ordered alternatives.

  For the adjustment of scores for tied values see Hajek, Sidak and Sen
  (1999), page 131ff.

}
\value{

  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue}, \code{statistic} and (optionally)
  \code{confint}.

}
\references{

        Bernd Streitberg & Joachim R\"ohmel (1986),
        Exact distributions for permutations and rank tests: 
        An introduction to some recently published algorithms.
        \emph{Statistical Software Newsletter} \bold{12}(1), 10--17.

        Bernd Streitberg & Joachim R\"ohmel (1987),
        Exakte Verteilungen f\"ur Rang- und Randomisierungstests
        im allgemeinen $c$-Stichprobenfall.
        \emph{EDV in Medizin und Biologie} \bold{18}(1), 12--19.

        Mark A. van de Wiel (2001), The split-up algorithm: a fast 
        symbolic method for computing p-values of rank statistics.
        \emph{Computational Statistics} \bold{16}, 519--538.

        David F. Bauer (1972), Constructing confidence sets using 
        rank statistics. \emph{Journal of the American Statistical Association}
        \bold{67}, 687--690.

        Jaroslav Hajek, Zbynek Sidak & Pranab K. Sen (1999), 
        \emph{Theory of Rank Tests}. San Diego, London: Academic Press.

}
\examples{

data(bloodp)

# Permutation test

# asymptotic p-value
perm_test(bp ~ group, data = bloodp)
# approximate p-value (with 99% confidence interval)
pvalue(perm_test(bp ~ group, data = bloodp, 
                 distribution = "approx", B = 9999))
# exact p-value
pt <- perm_test(bp ~ group, data = bloodp, distribution = "exact")
pvalue(pt)

# plot density and distribution of the standardised 
# test statistic
layout(matrix(1:2, nrow = 2))
s <- support(pt)
d <- sapply(s, function(x) dperm(pt, x))
p <- sapply(s, function(x) pperm(pt, x))
plot(s, d, type = "S", xlab = "Teststatistic", ylab = "Density")
plot(s, p, type = "S", xlab = "Teststatistic", ylab = "Cumm. Probability")


# Wilcoxon-Mann-Whitney test

# exact p-value and confidence interval for the difference in location
wt <- wilcox_test(bp ~ group, data = bloodp, 
            distribution = "exact", conf.int = TRUE)
pvalue(wt)
confint(wt)


# Kruskal-Wallis test (data in Hollander and Wolfe, page 200)
YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44, 
                             42, 60, 32, 42, 45, 58, 27, 51, 42, 52, 
                             38, 33, 26, 25, 28, 28, 26, 27, 27, 27, 
                             31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                  site = factor(c(rep("I", 10), rep("II", 10),
                                  rep("III", 10), rep("IV", 10))))

# approximate exact p-value
kw <- kruskal_test(length ~ site, data = YOY, 
             distribution = "approx", B = 9999)

kw
pvalue(kw)


# Nemenyi-Damico-Wolfe-Dunn test (joint ranking)
# Hollander and Wolfe, 1999, pages 200 and 244 
# (where Steel-Dwass results are given)
if (require(multcomp)) {

    NDWD <- perm_test(length ~ site, data = YOY,
        ytrafo = function(data) trafo(data, numeric_trafo = rank),
        xtrafo = function(data) trafo(data, factor_trafo = function(x)
            model.matrix(~x - 1) \%*\% t(contrMat(table(x), "Tukey"))),
        teststat = "maxtype", distribution = "approx", B = 90000)

    # global p-value
    print(pvalue(NDWD))

    # sites (I = II) != (III = IV) at alpha = 0.01 (page 244)
    print(pvalue(NDWD, adjusted = TRUE))
}

}
\keyword{htest}
