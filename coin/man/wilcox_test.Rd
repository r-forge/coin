\name{LocationTests}
\alias{perm_test}
\alias{wilcox_test}
\alias{normal_test}
\alias{median_test}
\title{ Independent Two- and K-Sample Location Tests }
\description{
    Testing the equality of the distributions of a numeric response in 
    two or more independent groups against shift alternatives.
}
\usage{
\method{perm_test}{formula}(formula, data, subset = NULL, \dots)
\method{perm_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"), ...)
\method{wilcox_test}{formula}(formula, data, subset = NULL, \dots)
\method{wilcox_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),
    conf.int = FALSE, conf.level = 0.95, ...)
\method{kruskal_test}{formula}(formula, data, subset = NULL, \dots)
\method{kruskal_test}{IndependenceProblem}(x, 
    distribution = c("asympt", "approx"), ...)
\method{normal_test}{formula}(formula, data, subset = NULL, \dots)
\method{normal_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),
    conf.int = FALSE, conf.level = 0.95, ...)
\method{median_test}{formula}(formula, data, subset = NULL, \dots)
\method{median_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),
    conf.int = FALSE, conf.level = 0.95, ...)
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
  \item{x}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or    
    \code{"less"}.  You can specify just the initial letter.}   
  \item{distribution}{the null distribution of the test statistic
    can be computed (\code{exact}) or can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{conf.int}{a logical indicating whether a confidence interval
    for the difference in location should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of \code{y} in
  the groups given by \code{x} is tested. In particular, the methods
  documented here are designed to detect shift alternatives. 

  

}
\value{
  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue}, \code{statistic} and (optionally)
  \code{confint}.

}
\examples{

data(bloodp)

# Permutation test

perm_test(bp ~ group, data = bloodp)
perm_test(bp ~ group, data = bloodp, distribution = "exact")
perm_test(bp ~ group, data = bloodp, distribution = "approx", B = 9999)

# Wilcoxon-Mann-Whitney test

wilcox_test(bp ~ group, data = bloodp)
wilcox_test(bp ~ group, data = bloodp, distribution = "exact")
wilcox_test(bp ~ group, data = bloodp, distribution = "approx", B = 9999)

# the v.d. Waerden test

normal_test(bp ~ group, data = bloodp)
normal_test(bp ~ group, data = bloodp, distribution = "exact")
normal_test(bp ~ group, data = bloodp, distribution = "approx", B = 9999)

# Median

median_test(bp ~ group, data = bloodp)
median_test(bp ~ group, data = bloodp, distribution = "exact")
median_test(bp ~ group, data = bloodp, distribution = "approx", B = 9999)

# Nemenyi-Damico-Wolfe-Dunn test (joint ranking)
# Hollander and Wolfe, 1999, pages 200 and 244 
# (where Steel-Dwass results are given)
if (require(multcomp)) {
    YOY <- data.frame(length = c(46, 28, 46, 37, 32, 41, 42, 45, 38, 44, 
                                 42, 60, 32, 42, 45, 58, 27, 51, 42, 52, 
                                 38, 33, 26, 25, 28, 28, 26, 27, 27, 27, 
                                 31, 30, 27, 29, 30, 25, 25, 24, 27, 30),
                      site = factor(c(rep("I", 10), rep("II", 10),
                                      rep("III", 10), rep("IV", 10))))

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
